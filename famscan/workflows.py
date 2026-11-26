"""Use-case workflows built from reusable stages.

Each workflow is an ordered list of stages operating on a shared `Ctx`. The
interactive menu (`famscan start`) and the non-interactive command
(`famscan workflow <name> --set key=value`) both run the same definitions.
"""
from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd

from . import analyze, features, motifs, msa, search, seqio, structure


# ======================================================================= core
@dataclass
class Param:
    name: str
    prompt: str
    kind: str = "file"          # file | files | str | int | float | bool | choice
    required: bool = False
    default: object = None
    choices: tuple = ()
    help: str = ""


@dataclass
class Workflow:
    key: str
    title: str
    question: str               # what the user is trying to do, shown in the menu
    params: list[Param]
    stages: list[str]
    notes: str = ""


class StopWorkflow(Exception):
    """Raised by a stage when there is nothing left to analyse."""


@dataclass
class Ctx:
    params: dict
    out: Path
    log: Callable = print
    proteins: dict = field(default_factory=dict)     # every protein loaded
    focus: list = field(default_factory=list)        # ids still in play
    table: pd.DataFrame = field(default_factory=pd.DataFrame)
    proteomes: dict = field(default_factory=dict)    # label -> {id: seq}
    sections: list = field(default_factory=list)     # (title, markdown)
    files: list = field(default_factory=list)
    timings: dict = field(default_factory=dict)
    aln: dict = field(default_factory=dict)
    ask: Callable | None = None          # set by `famscan start` to prompt the user
    outages: list = field(default_factory=list)

    def p(self, key, default=None):
        v = self.params.get(key)
        return default if v in (None, "", []) else v

    def seqs(self, ids=None) -> dict:
        return {i: self.proteins[i] for i in (self.focus if ids is None else ids) if i in self.proteins}

    def add(self, rows: dict[str, dict]):
        """Merge {protein_id: {col: value}} into the main table."""
        if not rows:
            return
        df = pd.DataFrame.from_dict(rows, orient="index")
        df.index.name = "id"
        overlap = [c for c in df.columns if c in self.table.columns]
        if overlap:
            self.table.update(df[overlap])
            df = df.drop(columns=overlap)
        self.table = self.table.join(df, how="left") if len(self.table) else df

    def section(self, title: str, md: str):
        self.sections.append((title, md))

    def service_down(self, service: str, stage: str, error: Exception, without: str,
                     fatal: bool = False) -> bool:
        """A web service used by this stage is unreachable. Record it, tell the
        user, and ask whether to carry on without it. Returns True to continue.

        `on_service_down` decides when nobody can be asked: "ask" (the default in
        `famscan start`), "skip" (carry on, the default for scripts) or "stop".
        """
        policy = str(self.p("on_service_down", "ask" if self.ask else "skip"))
        error = _plain_error(error)
        msg = f"{service} is not available: {str(error).strip()}"
        self.log(f"      {msg}")
        self.outages.append({"stage": stage, "service": service, "detail": str(error).strip(),
                             "continued_without": without})
        if policy == "stop":
            raise StopWorkflow(msg)
        if policy == "ask" and self.ask and not fatal:
            if not self.ask(f"\n  {service} could not be reached, so this run cannot {without}.\n"
                            f"  {str(error).strip()}\n  Continue without it?"):
                raise StopWorkflow(f"stopped at your request: {service} unavailable")
        tail = ("This stage was the whole search, so there is nothing to continue with."
                if fatal else
                f"The run continued without this step, so **{without}** is missing from the results "
                f"below. Everything else is unaffected.")
        self.section(f"{service} unavailable",
                     f"{msg}\n\n{tail} Rerun when the service is back, or avoid it altogether: "
                     f"`famscan fetchdb swissprot` for searches, `pdb_dir` for your own structures.")
        return not fatal

    def save(self, name: str, obj) -> Path:
        path = self.out / name
        if isinstance(obj, pd.DataFrame):
            obj.to_csv(path, sep="\t", index=False)
        elif isinstance(obj, (list, dict)) and name.endswith(".tsv"):
            pd.DataFrame(obj).to_csv(path, sep="\t", index=False)
        else:
            path.write_text(obj if isinstance(obj, str) else json.dumps(obj, indent=2, default=str))
        self.files.append(name)
        return path


def md_table(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df is None or not len(df):
        return "_none_"
    df = df.head(max_rows)
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "---|" * len(df.columns)

    def fmt(v):
        if isinstance(v, float):
            if pd.isna(v):
                return ""
            return f"{v:.1e}" if (v != 0 and abs(v) < 1e-3) else f"{v:.3g}"
        return str(v).replace("|", "/")[:60]

    return "\n".join([head, sep] + ["| " + " | ".join(fmt(v) for v in r) + " |"
                                    for r in df.itertuples(index=False)])


# ===================================================================== stages
STAGES: dict[str, tuple[Callable, str]] = {}


def stage(name: str, desc: str):
    def deco(fn):
        STAGES[name] = (fn, desc)
        return fn
    return deco


@stage("load", "Read input; call genes with Prodigal if it is a genome")
def st_load(ctx: Ctx):
    path = ctx.p("input")
    ctx.proteins, kind = seqio.to_proteins(path, log=ctx.log)
    ctx.focus = list(ctx.proteins)
    ctx.table = pd.DataFrame(index=pd.Index(ctx.focus, name="id"))
    ctx.table["length"] = [len(ctx.proteins[i]) for i in ctx.focus]
    if kind == "genome":
        seqio.write_fasta(ctx.proteins, ctx.out / "predicted_proteins.faa")
        ctx.files.append("predicted_proteins.faa")
    ctx.section("Input", f"`{path}` read as **{kind}**: {len(ctx.proteins)} proteins.")


@stage("load_many", "Read several genomes/proteomes (one per species or strain)")
def st_load_many(ctx: Ctx):
    paths = ctx.p("proteomes")
    paths = paths if isinstance(paths, list) else str(paths).split()
    expanded = []
    for p in paths:
        pp = Path(p)
        expanded += sorted(x for x in pp.iterdir() if x.suffix in (".faa", ".fa", ".fasta", ".fna", ".gbk")) \
            if pp.is_dir() else [pp]
    rows = []
    for p in expanded:
        seqs, kind = seqio.to_proteins(p, log=ctx.log)
        label = p.stem
        ctx.proteomes[label] = seqs
        for k, v in seqs.items():
            ctx.proteins[f"{label}|{k}"] = v
        rows.append({"proteome": label, "kind": kind, "proteins": len(seqs)})
    if len(ctx.proteomes) < 2:
        raise StopWorkflow("need at least two proteomes/genomes")
    ctx.focus = list(ctx.proteins)
    ctx.table = pd.DataFrame(index=pd.Index(ctx.focus, name="id"))
    ctx.table["proteome"] = [i.split("|", 1)[0] for i in ctx.focus]
    ctx.table["length"] = [len(ctx.proteins[i]) for i in ctx.focus]
    ctx.section("Inputs", md_table(pd.DataFrame(rows)))


@stage("family_search", "Profile/homology search for the family (HMMER)")
def st_family_search(ctx: Ctx):
    query = ctx.p("query")
    hits = search.dedupe(search.hmmer_search(query, ctx.seqs(), evalue=float(ctx.p("evalue", 1e-5)),
                                             workdir=ctx.out))
    if not hits:
        raise StopWorkflow(f"no hits for {query}")
    ctx.focus = [h.target for h in hits]
    ctx.add({h.target: {"evalue": h.evalue, "score": round(h.score, 1), "method": h.method,
                        "dom_from": h.domain_from, "dom_to": h.domain_to} for h in hits})
    ctx.table = ctx.table.loc[ctx.focus]
    msg = f"{len(hits)} proteins match `{query}` at E <= {ctx.p('evalue', 1e-5)}."
    if "proteome" in ctx.table:
        counts = ctx.table.groupby("proteome").size().rename("copies").reset_index()
        missing = [p for p in ctx.proteomes if p not in set(counts.proteome)]
        msg += "\n\nCopies per proteome:\n\n" + md_table(counts, 200)
        if missing:
            msg += f"\n\nNo member found in: {', '.join(missing)}"
    ctx.section("Family search", msg)


@stage("homology", "Homology search in any BLAST direction (blastp/n/x, tblastn, tblastx)")
def st_homology(ctx: Ctx):
    q, t = ctx.p("query"), ctx.p("target")
    prog, eng = ctx.p("program", "auto"), ctx.p("engine", "builtin")
    from .remote import SearchTimeout, ServiceUnavailable
    try:
        hits = search.homology_search(q, t, program=prog, engine=eng, evalue=float(ctx.p("evalue", 1e-5)),
                                      workdir=ctx.out, email=ctx.p("email"),
                                      timeout=float(ctx.p("timeout", 900)), log=ctx.log)
    except (ServiceUnavailable, SearchTimeout) as e:
        ctx.service_down({"ebi": "EBI", "mmseqs": "MMseqs2", "ncbi": "NCBI"}.get(eng, eng),
                         "homology", e, f"search {t} over the network", fatal=True)
        raise StopWorkflow(f"no hits: the {eng} service is unavailable. Download a database with "
                           f"`famscan fetchdb swissprot` and search it locally instead.") from None
    df = pd.DataFrame([h.row() for h in hits])
    ctx.save("homology_hits.tsv", df)
    if len(df):
        ctx.table = df.set_index(pd.Index([f"{r.query}->{r.target}" for r in df.itertuples()], name="id"))
        ctx.proteins.update({f"{h.query}->{h.target}": h.sequence for h in hits if h.sequence})
        ctx.focus = list(dict.fromkeys(ctx.table.index))
    ctx.section("Homology search", f"{len(df)} hits ({prog}, engine {eng}).\n\n" +
                md_table(df[[c for c in ["query", "target", "evalue", "score", "method", "domain_from",
                                         "domain_to", "frame", "nt_start", "nt_end", "strand", "pident"]
                             if c in df]] if len(df) else df))


@stage("name_by_homology", "Name proteins by best hit in a reference set (e.g. Swiss-Prot)")
def st_name(ctx: Ctx):
    ref = ctx.p("annotation_db")
    if not ref:
        ctx.section("Naming by homology", "Skipped: no `annotation_db` (e.g. uniprot_sprot.fasta) given.")
        return
    desc = seqio.read_descriptions(ref)
    q = seqio.write_fasta(ctx.seqs(), ctx.out / "_query.faa")
    hits = search.best_hit_per_query(search.builtin_search(str(q), ref, "blastp", evalue=1e-5))
    ctx.add({qid: {"best_hit": h.target, "best_hit_evalue": h.evalue,
                   "product": desc.get(h.target, "")[:80]} for qid, h in hits.items()})
    (ctx.out / "_query.faa").unlink()
    ctx.section("Naming by homology", f"{len(hits)}/{len(ctx.focus)} proteins have a hit in `{ref}`.")


@stage("domains", "Domain architecture (hmmscan-style against Pfam or your HMMs)")
def st_domains(ctx: Ctx):
    db = ctx.p("hmm_db")
    if not db:
        ctx.section("Domain architecture", "Skipped: give `hmm_db` (Pfam-A.hmm from "
                    "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) to annotate domains.")
        return
    rows = search.hmmscan_domains(ctx.seqs(), db)
    ctx.save("domains.tsv", rows)
    arch: dict[str, list] = {}
    for r in rows:
        arch.setdefault(r["protein"], []).append(f"{r['domain']}({r['start']}-{r['end']})")
    ctx.add({k: {"domains": " | ".join(v), "n_domains": len(v)} for k, v in arch.items()})
    top = pd.DataFrame(rows).domain.value_counts().head(20).rename_axis("domain").reset_index() if rows else None
    ctx.section("Domain architecture", f"{len(arch)}/{len(ctx.focus)} proteins have >= 1 domain.\n\n" + md_table(top))


@stage("primary", "Primary structure: ProtParam, composition, solubility")
def st_primary(ctx: Ctx):
    rows = {}
    for i, s in ctx.seqs().items():
        r = analyze.protparam(s)
        comp = features.composition(s)
        r.update({k: v for k, v in comp.items() if not k.startswith("pct_") or len(k) > 6})
        r.update(features.solubility_wh(s))
        rows[i] = r
    ctx.add(rows)
    t = ctx.table.loc[ctx.focus]
    ctx.section("Primary structure", md_table(t[["length", "mw_kda", "pI", "instability", "gravy"]]
                                              .describe().loc[["mean", "min", "max"]].round(2).reset_index()))


@stage("features", "Secondary structure, super-secondary motifs, TM helices, signal peptides, disorder, coiled coils")
def st_features(ctx: Ctx):
    rows, ss_lines = {}, []
    for i, s in ctx.seqs().items():
        f = features.profile(s)
        ss_lines.append(f">{i}\n{s}\n{f.pop('ss_pred')}")
        rows[i] = f
    ctx.add(rows)
    ctx.save("secondary_structure.txt", "\n".join(ss_lines) + "\n")
    t = ctx.table.loc[ctx.focus]
    loc = t.localisation_guess.value_counts().rename_axis("localisation").reset_index()
    ctx.section("Sequence features (fast predictors)",
                f"Signal peptide: {int(t.signal_peptide.sum())}, >=1 TM helix: {int((t.n_tm > 0).sum())}, "
                f"coiled coil: {int((t.coiled_coils != '').sum())}, >30 % disordered: {int((t.disorder_frac > 0.3).sum())}.\n\n"
                + md_table(loc) + "\n\nPer-residue H/E/C strings: `secondary_structure.txt`. These are quick "
                "heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.")


def _user_motifs(ctx) -> dict:
    raw = ctx.p("motifs")
    if not raw:
        return {}
    if isinstance(raw, dict):
        return raw
    out = {}
    for part in str(raw).split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            v = v.strip()
            is_prosite = re.fullmatch(r"<?[A-Zx\[\]{}()\d,]+(-[A-Zx\[\]{}()\d,]+)+>?", v)
            out[k.strip()] = motifs.prosite_to_regex(v) if is_prosite else v
    return out


@stage("motifs", "Motif scan (PROSITE patterns + your motifs) and de novo motif discovery")
def st_motifs(ctx: Ctx):
    pats = motifs.load_prosite_dat(ctx.p("prosite_dat")) if ctx.p("prosite_dat") else None
    user = _user_motifs(ctx)
    all_hits, rows = [], {}
    for i, s in ctx.seqs().items():
        h = motifs.scan(s, pats, user)
        for x in h:
            x["protein"] = i
        all_hits += h
        rows[i] = {"motifs": motifs.summarise(h), "n_motifs": len(h)}
    ctx.add(rows)
    ctx.save("motif_hits.tsv", all_hits)
    msg = (pd.DataFrame(all_hits).motif.value_counts().rename_axis("motif").reset_index()
           if all_hits else None)
    text = f"{len(all_hits)} motif matches.\n\n" + md_table(msg, 25)
    seqs = ctx.seqs()
    if 3 <= len(seqs) <= 300:
        dn = motifs.denovo_kmers(seqs)
        ctx.save("denovo_motifs.tsv", dn)
        text += "\n\n**De novo motifs shared by the set** (x = any residue):\n\n" + md_table(pd.DataFrame(dn), 10)
    ctx.section("Motifs", text)


@stage("compare_refs", "BLOSUM62 alignment to known family members (or the HMM consensus)")
def st_compare(ctx: Ctx):
    refs = seqio.read_seqs(ctx.p("references")) if ctx.p("references") else {}
    q = ctx.p("query")
    if not refs and q and search._is_hmm(q):
        refs = {f"{h.name}_consensus": h.consensus.upper() for h in search.load_hmms(q)}
    if not refs:
        return
    ctx.add({i: analyze.best_reference_match(s, refs) for i, s in ctx.seqs().items()})
    ctx.section("Similarity to references", f"Compared with {len(refs)} reference sequence(s).")


@stage("filter", "Apply your property filters")
def st_filter(ctx: Ctx):
    rules = ctx.p("rules")
    rules = json.loads(rules) if isinstance(rules, str) else (rules or {"length": [50, 3000]})
    rows, keep = {}, []
    for i in ctx.focus:
        ok, why = analyze.passes(ctx.table.loc[i].to_dict(), rules)
        rows[i] = {"pass_filter": ok, "fail_reason": why}
        if ok:
            keep.append(i)
    ctx.add(rows)
    ctx.section("Filters", f"`{rules}`: {len(keep)}/{len(ctx.focus)} pass.")
    ctx.focus = keep
    if not keep:
        raise StopWorkflow("no proteins pass the filters")


@stage("cluster", "Remove redundancy by clustering at % identity (CD-HIT style)")
def st_cluster(ctx: Ctx):
    ident = float(ctx.p("cluster_identity", 90))
    seqs = ctx.seqs()
    if len(seqs) > 3000:
        ctx.section("Clustering", "Skipped (>3,000 sequences): use MMseqs2 `easy-cluster`.")
        return
    cl = msa.cluster(seqs, min_identity=ident)
    rows = {m: {"cluster_rep": rep, "cluster_size": len(mem)} for rep, mem in cl.items() for m in mem}
    ctx.add(rows)
    ctx.focus = [i for i in ctx.focus if i in cl]
    ctx.section("Clustering", f"{len(seqs)} sequences -> {len(cl)} clusters at {ident:.0f} % identity; "
                "continuing with one representative per cluster.")


@stage("msa", "Multiple sequence alignment (FAMSA) + per-column conservation")
def st_msa(ctx: Ctx):
    seqs = ctx.seqs()
    ref_id = _reference_id(ctx)
    if ref_id and ref_id not in seqs:
        seqs = {ref_id: ctx.proteins[ref_id], **seqs}
    if len(seqs) < 3:
        ctx.section("Alignment", "Skipped: fewer than 3 sequences.")
        return
    if len(seqs) > 500:
        top = ctx.table.loc[ctx.focus].sort_values("score", ascending=False).index[:500] \
            if "score" in ctx.table else ctx.focus[:500]
        seqs = {i: ctx.proteins[i] for i in top}
    ctx.aln = msa.align(seqs)
    msa.write_aln(ctx.aln, ctx.out / "alignment.afa")
    ctx.files.append("alignment.afa")
    cons = msa.conservation(ctx.aln)
    ctx.save("conservation.tsv", cons)
    ref = ref_id or next(iter(ctx.aln))
    hot = msa.conserved_positions(ctx.aln, ref, 0.85)
    ctx.save("conserved_residues.tsv", hot)
    ctx.section("Alignment and conservation",
                f"{len(ctx.aln)} sequences, {len(cons)} columns. {len(hot)} residues of `{ref}` are "
                f"in columns with conservation >= 0.85 (candidate active/binding-site residues):\n\n"
                + ", ".join(h["residue"] for h in hot[:80]))


@stage("tree", "Phylogenetic tree (neighbour joining, BLOSUM62 distances)")
def st_tree(ctx: Ctx):
    if len(ctx.aln) < 3:
        return
    t = msa.tree(ctx.aln, method=ctx.p("tree_method", "nj"))
    ctx.save("tree.nwk", msa.tree_to_newick(t) + "\n")
    ctx.section("Phylogeny", "Newick: `tree.nwk` (open in iTOL or FigTree). For a maximum-likelihood "
                "tree run IQ-TREE on `alignment.afa`.\n\n```\n" + msa.tree_ascii(t)[:6000] + "```")


def _reference_id(ctx) -> str | None:
    r = ctx.p("reference_id")
    if r and r in ctx.proteins:
        return r
    if ctx.p("reference"):
        seqs = seqio.read_seqs(ctx.p("reference"))
        rid, rseq = next(iter(seqs.items()))
        ctx.proteins[rid] = rseq
        ctx.params["reference_id"] = rid
        return rid
    return None


@stage("variants", "Mine natural variation for engineering (back-to-consensus suggestions)")
def st_variants(ctx: Ctx):
    ref = _reference_id(ctx)
    if not ctx.aln or not ref or ref not in ctx.aln:
        ctx.section("Variant mining", "Skipped: needs an alignment that includes the reference.")
        return
    rows = msa.variant_table(ctx.aln, ref, min_freq=float(ctx.p("min_variant_freq", 0.15)))
    ctx.save("variants.tsv", rows)
    df = pd.DataFrame(rows)
    sugg = df[df.back_to_consensus] if len(df) else df
    ctx.section("Variant mining",
                f"{len(df)} positions of `{ref}` vary among {len(ctx.aln) - 1} homologs. "
                f"{len(sugg)} back-to-consensus mutations suggested (a common way to gain stability):\n\n"
                + md_table(sugg[["suggested", "consensus_freq", "conservation", "alternatives"]] if len(sugg) else sugg, 40))


@stage("subtract_host", "Remove proteins with homologs in the host (e.g. human) proteome")
def st_subtract(ctx: Ctx):
    host = ctx.p("host")
    if not host:
        raise StopWorkflow("drug-target screening needs `host` (e.g. the human proteome FASTA)")
    q = seqio.write_fasta(ctx.seqs(), ctx.out / "_q.faa")
    best = search.best_hit_per_query(search.builtin_search(str(q), host, "blastp",
                                                           evalue=float(ctx.p("host_evalue", 1e-4))))
    (ctx.out / "_q.faa").unlink()
    rows = {i: {"host_homolog": i in best, "host_best_hit": best[i].target if i in best else "",
                "host_evalue": best[i].evalue if i in best else None} for i in ctx.focus}
    ctx.add(rows)
    keep = [i for i in ctx.focus if i not in best]
    ctx.section("Host subtraction", f"{len(best)} proteins have a homolog in `{host}` and are dropped; "
                f"{len(keep)} pathogen-specific proteins remain.")
    ctx.focus = keep
    ess = ctx.p("essential_db")
    if ess and keep:
        q = seqio.write_fasta(ctx.seqs(), ctx.out / "_q.faa")
        eh = search.best_hit_per_query(search.builtin_search(str(q), ess, "blastp", evalue=1e-10))
        (ctx.out / "_q.faa").unlink()
        ctx.add({i: {"essential_hit": i in eh} for i in keep})
        ctx.section("Essentiality", f"{len(eh)} of them match the essential-gene set `{ess}` (e.g. DEG).")


def _screen(ctx, db: str, label: str, strict: float) -> dict:
    """Best hit of each focus protein in a screening db, with identity/coverage."""
    if search._is_hmm(db):
        out: dict = {}
        for h in search.hmmer_search(db, ctx.seqs(), evalue=1e-10):
            if h.target not in out or h.evalue < out[h.target][f"{label}_evalue"]:
                out[h.target] = {f"{label}_hit": h.query, f"{label}_evalue": h.evalue,
                                 f"{label}_class": "profile"}
        return out
    dbseqs = seqio.read_seqs(db)
    desc = seqio.read_descriptions(db)
    q = seqio.write_fasta(ctx.seqs(), ctx.out / "_q.faa")
    best = search.best_hit_per_query(search.builtin_search(str(q), db, "blastp", evalue=1e-5))
    (ctx.out / "_q.faa").unlink()
    rows = {}
    for i, h in best.items():
        a = analyze.blosum_align(ctx.proteins[i], dbseqs[h.target])
        cls = ("perfect" if a.identity >= 99.5 and a.coverage >= 95 else
               "strict" if a.identity >= strict and a.coverage >= 60 else "loose")
        rows[i] = {f"{label}_hit": h.target, f"{label}_desc": desc.get(h.target, "")[:70],
                   f"{label}_identity": a.identity, f"{label}_coverage": a.coverage,
                   f"{label}_evalue": h.evalue, f"{label}_class": cls}
    return rows


@stage("resistance", "Screen for antibiotic-resistance genes (CARD / AMRFinder db)")
def st_resistance(ctx: Ctx):
    db = ctx.p("resistance_db")
    if not db:
        raise StopWorkflow("give `resistance_db`: CARD protein_fasta_protein_homolog_model.fasta "
                           "(https://card.mcmaster.ca/download) or NCBI AMRFinder HMMs")
    rows = _screen(ctx, db, "amr", float(ctx.p("strict_identity", 80)))
    ctx.add(rows)
    ctx.focus = [i for i in ctx.focus if i in rows]
    df = pd.DataFrame.from_dict(rows, orient="index").rename_axis("protein").reset_index()
    ctx.section("Antibiotic resistance",
                f"{len(rows)} proteins match `{db}`. Classes follow CARD RGI logic: perfect (identical), "
                f"strict (>= {ctx.p('strict_identity', 80)} % identity, >= 60 % coverage), loose (weaker homologs, "
                "possibly novel variants).\n\n" + md_table(df, 50))


@stage("toxin_allergen", "Screen for toxins and allergens (FAO/WHO Codex rules)")
def st_toxin(ctx: Ctx):
    text = []
    if ctx.p("toxin_db"):
        rows = _screen(ctx, ctx.p("toxin_db"), "toxin", float(ctx.p("strict_identity", 50)))
        ctx.add(rows)
        text.append(f"**Toxins:** {len(rows)} proteins match `{ctx.p('toxin_db')}`.")
    if ctx.p("allergen_db"):
        alg = seqio.read_seqs(ctx.p("allergen_db"))
        k8 = {}
        for aid, s in alg.items():
            for j in range(len(s) - 7):
                k8.setdefault(s[j:j + 8], aid)
        q = seqio.write_fasta(ctx.seqs(), ctx.out / "_q.faa")
        best = search.best_hit_per_query(search.builtin_search(str(q), ctx.p("allergen_db"), "blastp", evalue=1.0))
        (ctx.out / "_q.faa").unlink()
        rows = {}
        for i, s in ctx.seqs().items():
            exact = next((k8[s[j:j + 8]] for j in range(len(s) - 7) if s[j:j + 8] in k8), None)
            r = {"allergen_8mer": exact or ""}
            if i in best:
                a = analyze.blosum_align(s, alg[best[i].target])
                r.update(allergen_hit=best[i].target, allergen_identity=a.identity, allergen_aln_len=a.aln_len)
                r["allergen_flag"] = bool(exact) or (a.aln_len >= 80 and a.identity > 35)
            else:
                r["allergen_flag"] = bool(exact)
            rows[i] = r
        ctx.add(rows)
        n = sum(r["allergen_flag"] for r in rows.values())
        text.append(f"**Allergens:** {n} proteins flagged by Codex criteria (> 35 % identity over >= 80 aa, "
                    "or an identical 8-residue stretch) against `" + ctx.p("allergen_db") + "`.")
    t = ctx.table.loc[ctx.focus]
    if {"signal_peptide", "n_cys"} <= set(t.columns):
        scr = t[(t.length <= 150) & (t.n_cys >= 4) & t.signal_peptide.astype(bool)]
        ctx.add({i: {"small_secreted_cys_rich": True} for i in scr.index})
        text.append(f"**Toxin-like by sequence:** {len(scr)} small (<= 150 aa), secreted, Cys-rich (>= 4 Cys) "
                    "proteins - the profile of many peptide toxins.")
    flag_cols = [c for c in ("toxin_hit", "allergen_flag", "small_secreted_cys_rich") if c in ctx.table]
    if flag_cols:
        t = ctx.table.loc[ctx.focus, flag_cols]
        hit = t.apply(lambda r: any(isinstance(v, str) and v or v is True for v in r), axis=1)
        ctx.focus = list(t.index[hit])
        text.append(f"**{len(ctx.focus)} proteins flagged in total** (listed under Final candidates).")
    if not text:
        raise StopWorkflow("give `toxin_db` (e.g. UniProt keyword KW-0800 FASTA) and/or `allergen_db` "
                           "(e.g. COMPARE or AllergenOnline FASTA)")
    ctx.section("Toxins and allergens", "\n\n".join(text))


def _plain_error(e: Exception) -> Exception:
    """Swap requests' long connection-pool messages for one plain sentence."""
    import requests
    from urllib.parse import urlparse
    url = getattr(getattr(e, "request", None), "url", "") or ""
    host = urlparse(url).hostname or "the site"
    if isinstance(e, requests.ConnectionError):
        return RuntimeError(f"could not connect to {host} (the site is down, or no internet connection)")
    if isinstance(e, requests.Timeout):
        return RuntimeError(f"{host} did not answer in time")
    if isinstance(e, ValueError) and not isinstance(e, requests.RequestException):
        return RuntimeError("the site sent back an error page instead of results")
    return e


def _is_outage(e: Exception) -> bool:
    """Connection failure, time-out, 5xx/429, or an HTML error page instead of JSON."""
    import requests
    if isinstance(e, (requests.ConnectionError, requests.Timeout, ValueError)):  # JSONDecodeError is a ValueError
        return True
    resp = getattr(e, "response", None)
    return resp is not None and (resp.status_code >= 500 or resp.status_code == 429)


@stage("structure", "Tertiary structure for top hits (AlphaFold DB / ESMFold / your PDBs), pLDDT, TM-score, disulfides")
def st_structure(ctx: Ctx):
    method = ctx.p("structure", "auto")
    if method == "none":
        ctx.section("Tertiary structure", "Skipped (structure=none).")
        return
    n = int(ctx.p("top_n", 10))
    t = ctx.table.loc[ctx.focus]
    key = next((k for k in ("rank_score", "score") if k in t), None)
    top = list(t.sort_values(key, ascending=False).index[:n]) if key else ctx.focus[:n]
    ref_pdb = ctx.p("reference_pdb")
    rows, ok = {}, 0
    foldseek_down = False
    for i in top:
        pdb, src = structure.get_structure(i, ctx.proteins[i], ctx.out / "structures",
                                           Path(ctx.p("pdb_dir")) if ctx.p("pdb_dir") else None, method,
                                           log=ctx.log)
        if src == "unreachable":
            ctx.service_down("The structure service (ESMFold / AlphaFold DB)", "structure",
                             RuntimeError("the prediction servers could not be reached"),
                             "confirm these hits with a 3D structure")
            ctx.add(rows)
            return
        r = {"structure_source": src}
        if pdb:
            ok += 1
            ss = features.ss_from_structure(pdb)
            r.update(pdb=str(pdb), plddt=structure.plddt(pdb), ss_3d_helix=features.ss_summary(ss)["ss_helix"],
                     ss_3d_strand=features.ss_summary(ss)["ss_strand"],
                     supersecondary_3d=",".join(sorted({m["motif"] for m in features.supersecondary(ss)})),
                     disulfides=";".join(f"C{a}-C{b}" for a, b in features.disulfides_from_structure(pdb)))
            if ref_pdb:
                r.update(structure.tm_score(pdb, Path(ref_pdb)))
            if ctx.p("foldseek_web") and not foldseek_down:
                try:
                    fs = structure.foldseek_web(pdb)
                    r["foldseek_top"] = fs[0]["target"] if fs else ""
                except Exception as e:  # noqa: BLE001
                    if _is_outage(e):  # the site, not this protein: report once, ask, stop calling it
                        foldseek_down = True
                        ctx.service_down("Foldseek", "structure", e, "compare structures with known folds "
                                         "(foldseek_top)")
                    else:
                        ctx.log(f"      foldseek failed: {str(e)[:100]}")
        rows[i] = r
    ctx.add(rows)
    plddt_min, tm_min = float(ctx.p("plddt_min", 70)), float(ctx.p("tm_min", 0.5))

    def verdict(r):
        if pd.isna(r.get("plddt", float("nan"))):
            return "no structure"
        if r["plddt"] < plddt_min:
            return "low-confidence fold"
        if "tm_score" in r and not pd.isna(r.get("tm_score")):
            return "confirmed same fold" if r["tm_score"] >= tm_min else "different fold"
        return "well folded"
    ctx.add({i: {"structure_verdict": verdict(ctx.table.loc[i].to_dict())} for i in top})
    msg = f"{ok}/{len(top)} structures obtained."
    if ok < len(top):
        msg += (" Missing ones need network access to AlphaFold DB / ESMFold, or fold them with "
                "ColabFold and pass the folder as `pdb_dir`.")
    ctx.section("Tertiary structure", msg)


@stage("rank", "Rank candidates for the lab (stability, solubility, fold confidence)")
def st_rank(ctx: Ctx):
    t = ctx.table.loc[ctx.focus]
    lo, hi = (int(x) for x in str(ctx.p("length_range", "100-600")).split("-"))

    def score(r):
        s = 0.0
        s += 1.0 if r.get("instability", 99) < 40 else 0
        s += float(r.get("wh_probability", 0.5)) if r.get("wh_soluble", False) else -float(r.get("wh_probability", 0.5))
        s -= 1.0 * float(r.get("disorder_frac", 0) or 0)
        s -= 0.75 if (r.get("n_tm", 0) or 0) > 0 else 0
        s -= 0.5 if not (lo <= r.get("length", 0) <= hi) else 0
        s += 0.5 if r.get("host_homolog") is False else 0
        s += 0.5 if r.get("essential_hit") else 0
        if not pd.isna(r.get("plddt", float("nan"))):
            s += (float(r["plddt"]) - 70) / 30
        if not pd.isna(r.get("score", float("nan"))):
            s += min(float(r["score"]) / 500, 1)
        return round(s, 2)
    ctx.add({i: {"rank_score": score(t.loc[i].to_dict())} for i in ctx.focus})
    t = ctx.table.loc[ctx.focus].sort_values("rank_score", ascending=False)
    cols = [c for c in ["rank_score", "length", "pI", "instability", "wh_soluble", "disorder_frac", "n_tm",
                        "localisation_guess", "plddt", "host_homolog", "essential_hit", "domains"] if c in t]
    ctx.section("Ranking", "Higher is better: stable (instability < 40), predicted soluble, ordered, no TM "
                f"helices, length {lo}-{hi} aa, confident fold.\n\n" + md_table(t[cols].reset_index(), 25))


@stage("strain_compare", "Compare strains: reciprocal best hits, core / accessory / unique genes")
def st_strains(ctx: Ctx):
    labels = list(ctx.proteomes)
    ref = ctx.p("reference_strain", labels[0])
    tmp = ctx.out / "_strains"
    tmp.mkdir(exist_ok=True)
    paths = {l: seqio.write_fasta(s, tmp / f"{l}.faa") for l, s in ctx.proteomes.items()}
    presence = pd.DataFrame(index=pd.Index(list(ctx.proteomes[ref]), name=f"{ref}_gene"))
    unique = {}
    for l in labels:
        if l == ref:
            presence[l] = True
            continue
        ctx.log(f"      RBH {ref} <-> {l}")
        fwd = search.best_hit_per_query(search.builtin_search(str(paths[ref]), str(paths[l]), "blastp", 1e-10))
        rev = search.best_hit_per_query(search.builtin_search(str(paths[l]), str(paths[ref]), "blastp", 1e-10))
        rbh = {q: h.target for q, h in fwd.items() if h.target in rev and rev[h.target].target == q}
        presence[l] = [g in rbh for g in presence.index]
        unique[l] = [g for g in ctx.proteomes[l] if g not in rev]
    n = len(labels)
    k = presence.sum(axis=1)
    summary = pd.DataFrame([{"category": "core (all strains)", "genes": int((k == n).sum())},
                            {"category": "accessory", "genes": int(((k > 1) & (k < n)).sum())},
                            {"category": f"only in {ref}", "genes": int((k == 1).sum())}]
                           + [{"category": f"only in {l} (no hit in {ref})", "genes": len(u)} for l, u in unique.items()])
    presence.reset_index().to_csv(ctx.out / "presence_absence.tsv", sep="\t", index=False)
    ctx.files.append("presence_absence.tsv")
    ctx.save("strain_unique_genes.json", unique)
    for p in paths.values():
        p.unlink()
    tmp.rmdir()
    ctx.section("Strain comparison", f"Reference strain `{ref}`; orthologs by reciprocal best hit (E <= 1e-10).\n\n"
                + md_table(summary, 50))
    ctx.focus = [f"{ref}|{g}" for g in presence.index[k < n]][:5000]


# ================================================================== workflows
P_INPUT = Param("input", "Genome (DNA FASTA/GenBank) or proteome (protein FASTA)", required=True)
P_QUERY = Param("query", "Family: Pfam accession, .hmm file, or FASTA of example proteins", required=True)
P_EVAL = Param("evalue", "E-value cut-off", "float", default=1e-5)
P_TOPN = Param("top_n", "How many top hits get a 3D structure", "int", default=10)
P_STRUCT = Param("structure", "Structure source", "choice", default="auto", choices=("auto", "afdb", "esmfold", "none"))
P_HMMDB = Param("hmm_db", "Domain library for architecture (Pfam-A.hmm), optional")
P_REFPDB = Param("reference_pdb", "Known structure (PDB) to compare folds with, optional")
P_RULES = Param("rules", 'Filters as JSON, e.g. {"length":[200,600],"pI":[null,7]}', "str",
                default='{"length":[50,3000]}')
P_DOWN = Param("on_service_down", "If a web service is down", "choice", default="ask",
               choices=("ask", "skip", "stop"),
               help="ask (famscan start), skip (carry on, default for scripts) or stop")
P_MOTIFS = Param("motifs", "Your own motifs, name=PATTERN;... (regex or PROSITE syntax), optional", "str")

WORKFLOWS: dict[str, Workflow] = {w.key: w for w in [
    Workflow("family", "Find members of a protein family",
             "I want every member of a protein family in my genome/proteome, checked by structure",
             [P_INPUT, P_QUERY, Param("references", "FASTA of known members, optional"), P_REFPDB, P_RULES,
              P_EVAL, P_HMMDB, P_MOTIFS, P_TOPN, P_STRUCT, P_DOWN],
             ["load", "family_search", "domains", "primary", "features", "motifs", "compare_refs", "filter",
              "structure", "report"]),
    Workflow("enzyme", "Discover new enzymes",
             "I'm looking for new enzymes of a known type (e.g. novel lipases, PETases)",
             [P_INPUT, P_QUERY, Param("references", "FASTA of characterised enzymes, optional"), P_REFPDB,
              P_RULES, P_MOTIFS, P_HMMDB, Param("cluster_identity", "Cluster identity % to remove redundancy", "float",
                                                default=90), P_EVAL, P_TOPN, P_STRUCT, P_DOWN],
             ["load", "family_search", "domains", "primary", "features", "motifs", "compare_refs", "filter",
              "cluster", "msa", "tree", "structure", "rank", "report"],
             "Use `motifs` for catalytic signatures, e.g. lipase=G-x-S-x-G."),
    Workflow("annotate", "Annotate a new genome",
             "I have a newly sequenced genome and want to know what its proteins are",
             [P_INPUT, Param("annotation_db", "Reference proteins for naming (e.g. uniprot_sprot.fasta), optional"),
              P_HMMDB, P_MOTIFS],
             ["load", "name_by_homology", "domains", "primary", "features", "motifs", "report"],
             "For a fully curated annotation also run Bakta or Prokka; eukaryotes need BRAKER first."),
    Workflow("drug_target", "Find drug targets",
             "I want proteins in a pathogen that humans don't have (potential drug or vaccine targets)",
             [P_INPUT, Param("host", "Host proteome FASTA (e.g. UniProt human UP000005640)", required=True),
              Param("essential_db", "Essential genes FASTA (e.g. DEG), optional"), P_HMMDB, P_TOPN,
              P_STRUCT, P_DOWN],
             ["load", "subtract_host", "domains", "primary", "features", "structure", "rank", "report"],
             "Secreted/membrane proteins in the output are vaccine candidates; essential, cytoplasmic "
             "enzymes with good folds are classic drug targets."),
    Workflow("resistance", "Screen for antibiotic resistance",
             "I want to know which antibiotic-resistance genes a bacterium carries",
             [P_INPUT, Param("resistance_db", "CARD protein FASTA or AMRFinder HMMs", required=True),
              Param("strict_identity", "% identity for a 'strict' match", "float", default=80)],
             ["load", "resistance", "primary", "features", "report"]),
    Workflow("toxin_allergen", "Screen for toxins and allergens",
             "I need to check proteins for toxin or allergen risk (safety, food, biotech)",
             [P_INPUT, Param("toxin_db", "Toxin protein FASTA (UniProt KW-0800), optional"),
              Param("allergen_db", "Allergen FASTA (COMPARE / AllergenOnline), optional")],
             ["load", "primary", "features", "toxin_allergen", "motifs", "report"]),
    Workflow("evolution", "Study a family's evolution across species",
             "I want to compare a protein family across several species (copy number, conservation, tree)",
             [Param("proteomes", "Genomes/proteomes, one per species (paths or a folder)", "files", required=True),
              P_QUERY, P_EVAL, P_MOTIFS, Param("tree_method", "Tree method", "choice", default="nj", choices=("nj", "upgma"))],
             ["load_many", "family_search", "primary", "features", "motifs", "msa", "tree", "report"]),
    Workflow("unknown", "Characterise unknown proteins",
             "I have proteins with no or weak annotation and want to know what they might be",
             [P_INPUT, P_HMMDB, Param("annotation_db", "Reference proteins for naming, optional"), P_MOTIFS, P_TOPN,
              P_STRUCT, Param("foldseek_web", "Search predicted structures with Foldseek web", "bool",
                              default=False), P_DOWN],
             ["load", "name_by_homology", "domains", "primary", "features", "motifs", "structure", "report"],
             "Structure search (Foldseek) often finds function when sequence search fails."),
    Workflow("lab_candidates", "Pick candidates for the lab",
             "I have candidate proteins and need to choose which to clone and express",
             [P_INPUT, Param("length_range", "Preferred length range (aa)", "str", default="100-600"), P_TOPN,
              P_STRUCT, P_RULES],
             ["load", "primary", "features", "filter", "rank", "structure", "rank", "report"]),
    Workflow("engineering", "Mine variants for protein engineering",
             "I want to improve a protein (stability, activity) using natural sequence variation",
             [Param("reference", "Your protein (FASTA, first record)", required=True),
              P_INPUT, Param("query", "Family HMM/Pfam to collect homologs, optional (default: your protein)"),
              Param("cluster_identity", "Cluster identity % to reduce redundancy", "float", default=95),
              Param("min_variant_freq", "Minimum frequency of a substitution", "float", default=0.15), P_EVAL],
             ["load", "family_search", "cluster", "msa", "variants", "tree", "report"],
             "`input` is the sequence database to mine homologs from (e.g. UniRef50 subset)."),
    Workflow("strains", "Compare strains",
             "I have several strains/isolates and want shared and unique genes",
             [Param("proteomes", "Genomes/proteomes, one per strain (paths or a folder)", "files", required=True),
              Param("reference_strain", "Reference strain label (file name without extension), optional", "str"),
              P_HMMDB],
             ["load_many", "strain_compare", "domains", "features", "report"]),
    Workflow("search", "Run a BLAST-style search (any direction)",
             "I just want a homology search: protein or DNA, translated or not",
             [Param("query", "Query FASTA (protein or DNA) or HMM", required=True),
              Param("target", "Target FASTA (protein or DNA), or a db name for engine ncbi/blast+", required=True),
              Param("program", "Program", "choice", default="auto",
                    choices=("auto", "blastp", "blastn", "blastx", "tblastn", "tblastx", "hmmsearch")),
              Param("engine", "Engine", "choice", default="builtin",
                    choices=("builtin", "blast+", "diamond", "ncbi", "ebi", "mmseqs")),
              Param("email", "Contact email (required by the EBI service only)", "str"),
              Param("timeout", "Give up on a web service after this many seconds", "int", default=900),
              P_DOWN, P_EVAL],
             ["homology", "report"],
             "builtin needs no install and no network. ncbi/ebi/mmseqs are web services, fine for a few "
             "queries; for big jobs run `famscan fetchdb swissprot` and search the file locally."),
]}


# ===================================================================== runner
def run_workflow(key: str, params: dict, out: str | Path = "famscan_out", log=print,
                 ask: Callable | None = None) -> Ctx:
    """Run workflow `key`. params["formats"] = "md", "xlsx", "pdf", a comma list or "all"."""
    wf = WORKFLOWS[key]
    from .export import parse_formats
    parse_formats(params.get("formats"))  # fail fast on a typo
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    missing = [p.name for p in wf.params if p.required and not params.get(p.name)]
    if missing:
        raise ValueError(f"workflow '{key}' needs: {', '.join(missing)}")
    defaults = {p.name: p.default for p in wf.params if p.default is not None}
    if key == "engineering" and not params.get("query"):
        params["query"] = params["reference"]
    from . import config
    merged = {**defaults, **params}
    if merged.get("engine") == "ebi":
        merged["email"] = config.email(merged.get("email"))
    ctx = Ctx(params=merged, out=out, log=log, ask=ask)
    t0 = time.time()
    stopped = None
    for i, name in enumerate(wf.stages, 1):
        if name == "report":
            break
        fn, desc = STAGES[name]
        log(f"[{i}/{len(wf.stages) - 1}] {desc}")
        t = time.time()
        try:
            fn(ctx)
        except StopWorkflow as e:
            stopped = str(e)
            log(f"      stopped: {e}")
            break
        ctx.timings[name] = round(time.time() - t, 2)
    write_report(ctx, wf, stopped, time.time() - t0)
    from .export import export
    export(out, ctx.p("formats", "md"), log=log)
    log(f"done in {time.time() - t0:.1f}s -> {out / 'report.md'}")
    return ctx


PRIVATE_PARAMS = ("email",)  # personal details: used for the run, never written to its outputs


def public_params(params: dict) -> dict:
    """Parameters safe to write into run.json / report.md / report.xlsx / report.pdf,
    which people share and sometimes commit."""
    return {k: ("(set)" if k in PRIVATE_PARAMS and v else v) for k, v in params.items()}


def write_report(ctx: Ctx, wf: Workflow, stopped: str | None, elapsed: float):
    if len(ctx.table):
        ctx.table.reset_index().to_csv(ctx.out / "results.tsv", sep="\t", index=False)
        ctx.files.insert(0, "results.tsv")
        seqio.write_fasta(ctx.seqs(), ctx.out / "final_candidates.faa")
        ctx.files.insert(1, "final_candidates.faa")
    shown = public_params(ctx.params)
    (ctx.out / "run.json").write_text(json.dumps({"workflow": wf.key, "params": shown,
                                                  "timings": ctx.timings, "outages": ctx.outages},
                                                 indent=2, default=str))
    lines = [f"# famscan: {wf.title}", "", f"_{wf.question}_", "",
             f"- Stages: {' -> '.join(wf.stages[:-1])}",
             f"- Parameters: `{ {k: v for k, v in shown.items() if v not in (None, '')} }`",
             f"- Run time: {elapsed:.1f} s; final candidates: {len(ctx.focus)}"]
    if stopped:
        lines += ["", f"> **Stopped early:** {stopped}"]
    if ctx.outages:
        lines += ["", "> **Steps missing because a web service was down:** "
                  + "; ".join(f"{o['service']} ({o['stage']}) - no {o['continued_without']}"
                              for o in ctx.outages)]
    if wf.notes:
        lines += ["", f"> {wf.notes}"]
    for title, md in ctx.sections:
        lines += ["", f"## {title}", "", md]
    if len(ctx.table) and ctx.focus:
        t = ctx.table.loc[[i for i in ctx.focus if i in ctx.table.index]]
        sort = next((k for k in ("rank_score", "score") if k in t), None)
        if sort:
            t = t.sort_values(sort, ascending=False)
        cols = [c for c in ["rank_score", "evalue", "score", "length", "pI", "instability", "localisation_guess",
                            "n_tm", "signal_peptide", "motifs", "domains", "product", "ref_identity", "plddt",
                            "tm_score", "structure_verdict", "amr_class", "amr_desc", "allergen_flag"] if c in t]
        lines += ["", "## Final candidates", "", md_table(t[cols].reset_index(), 40)]
    lines += ["", "## Files", ""] + [f"- `{f}`" for f in ctx.files]
    (ctx.out / "report.md").write_text("\n".join(lines) + "\n")
