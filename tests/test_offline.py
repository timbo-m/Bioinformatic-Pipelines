"""Offline checks: python -m pytest tests/"""
from pathlib import Path

import numpy as np
import pytest

from famscan import analyze, structure

AA = "ACDEFGHIKLMNPQRSTVWY"
THREE = dict(zip(AA, "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()))


def _helix_pdb(path: Path, seq: str, b: float) -> Path:
    lines = []
    for i, aa in enumerate(seq):
        t = np.radians(100 * i)
        x, y, z = 2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * i
        lines.append(f"ATOM  {i+1:5d}  CA  {THREE[aa]} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}           C")
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_protparam():
    p = analyze.protparam("MKTAYIAKQRQISFVKSHFSRQ")
    assert p["length"] == 22 and 9 < p["pI"] < 12


def test_blosum_identity():
    r = analyze.blosum_align("MKTAYIAKQRQISFVK", "MKTAYIAKQRQISFVK")
    assert r.identity == 100 and r.coverage == 100


def test_structure_helpers(tmp_path):
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLG"
    a = _helix_pdb(tmp_path / "a.pdb", seq, 0.85)   # ESMFold-style 0-1 pLDDT
    assert structure.plddt(a) == 85.0
    assert structure.tm_score(a, a)["tm_score"] > 0.99


def test_uniprot_parse():
    assert structure.uniprot_accession("sp|P0A9T0|SERA_ECOLI") == "P0A9T0"
    assert structure.uniprot_accession("938293.PRJEB85.HG003685_328") is None


# ----------------------------------------------------------- v0.2 additions
DATA = Path(__file__).resolve().parents[1] / "demo" / "data"


def test_six_frame_coordinates():
    from famscan import seqio
    nt = "ATG" + "GCT" * 40 + "TAA"
    pep = seqio.six_frame({"x": nt}, min_len=30)
    k = next(k for k in pep if "|f+1|" in k)
    assert pep[k].startswith("MAAA") and k.endswith("|1-123")


def test_all_blast_directions():
    from famscan import search
    cases = {
        "blastp": (DATA / "luxC_reference.faa", DATA / "luxC_homologs.faa"),
        "blastn": (DATA / "bmyD_gene.fna", DATA / "bacillus_contig.fna"),
        "blastx": (DATA / "bmyD_gene.fna", DATA / "bacillus_proteins.faa"),
        "tblastn": (DATA / "bmyD_protein.faa", DATA / "bacillus_contig.fna"),
        "tblastx": (DATA / "bmyD_gene.fna", DATA / "bacillus_contig.fna"),
    }
    for prog, (q, t) in cases.items():
        hits = search.homology_search(str(q), str(t), program=prog)
        assert hits, prog
    # translated hits map back to the same genome locus
    tn = search.homology_search(str(cases["tblastn"][0]), str(cases["tblastn"][1]), program="tblastn")[0]
    bn = search.homology_search(str(cases["blastn"][0]), str(cases["blastn"][1]), program="blastn")[0]
    assert abs(tn.extra["nt_start"] - bn.domain_from) < 30


def test_prosite_conversion_and_scan():
    from famscan import motifs
    assert motifs.prosite_to_regex("C-x(2,4)-[DE]-{P}>") == "C.{2,4}[DE][^P]$"
    hits = motifs.scan("MAGQSGGKSTLLK")
    assert any(h["id"] == "PS00017" for h in hits)


def test_features_known_proteins():
    from famscan import features
    ompa = "MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKG"
    assert features.signal_peptide(ompa)["sp_cleavage"] == 21
    gcn4 = "MKVKQLEDKVEELLSKNYHLENEVARLKKLVGERMKQLEDKVEELLSKNYHLENEVARLKKLVGER"
    assert features.coiled_coils(gcn4)
    assert "H" in features.secondary_structure("AEELLKKAEELLKKAEELLKK")


def test_msa_tree_variants():
    from famscan import msa, seqio
    s = seqio.read_seqs(DATA / "luxC_homologs.faa")
    aln = msa.align(s)
    assert len({len(v) for v in aln.values()}) == 1
    assert msa.tree_to_newick(msa.tree(aln)).endswith(";")
    assert msa.variant_table(aln, next(iter(aln)))


def test_workflow_family(tmp_path):
    from famscan.workflows import run_workflow
    ctx = run_workflow("family", {"input": str(DATA / "proteome_938293.faa"), "query": str(DATA / "PF02826.hmm"),
                                  "structure": "none"}, tmp_path, log=lambda *_: None)
    assert len(ctx.focus) == 7 and (tmp_path / "report.md").exists()


def test_workflow_genome_input(tmp_path):
    from famscan.workflows import run_workflow
    ctx = run_workflow("enzyme", {"input": str(DATA / "bacillus_contig.fna"),
                                  "query": str(DATA / "Thioesterase.hmm"), "structure": "none"},
                       tmp_path, log=lambda *_: None)
    assert ctx.focus and len(ctx.proteins) > 200  # genes were called from DNA


def test_export_formats(tmp_path):
    from famscan.export import parse_formats
    from famscan.workflows import run_workflow
    assert parse_formats("excel") == ["md", "xlsx"]
    assert parse_formats("all") == ["md", "xlsx", "pdf"]
    run_workflow("family", {"input": str(DATA / "proteome_938293.faa"), "query": str(DATA / "PF02826.hmm"),
                            "structure": "none", "formats": "all"}, tmp_path, log=lambda *_: None)
    for name in ("report.md", "report.xlsx", "report.pdf"):
        assert (tmp_path / name).stat().st_size > 1000


def test_database_catalog_and_dir(tmp_path, monkeypatch):
    from famscan import databases
    monkeypatch.setenv("FAMSCAN_DB", str(tmp_path / "db"))
    assert databases.db_dir().is_dir()
    assert "swissprot" in databases.describe()
    (tmp_path / "db" / "uniprot_sprot.fasta").write_text(">a\nMKV\n")
    assert databases.installed()[0]["name"] == "uniprot_sprot.fasta"
    try:
        databases.fetch("not_a_db")
    except ValueError as e:
        assert "unknown database" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_remote_timeout_is_bounded():
    from famscan.remote import SearchTimeout, _poll
    try:
        _poll(lambda: "WAITING", {"READY"}, {"FAILED"}, timeout=1.0, interval=0.2,
              service="test", job_id="J1", log=lambda *_: None)
    except SearchTimeout as e:
        assert e.job_id == "J1" and "local database" in str(e)
    else:
        raise AssertionError("expected SearchTimeout")


def test_ebi_requires_email():
    from famscan.remote import ebi_blast
    try:
        ebi_blast("MKV", email=None)
    except ValueError as e:
        assert "email" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_builtin_hits_carry_descriptions():
    from famscan import search
    hits = search.builtin_search(str(DATA / "luxC_reference.faa"), str(DATA / "luxC_homologs.faa"), "blastp")
    assert hits and all(h.description for h in hits)


def test_remote_outage_is_a_clear_error(monkeypatch):
    from famscan import remote

    class Down:
        status_code, headers = 503, {"content-type": "text/html"}
        text = "<html><title>MMseqs2 search will be back!</title></html>"

    monkeypatch.setattr(remote.requests, "request", lambda *a, **k: Down())
    with pytest.raises(RuntimeError, match=r"unavailable \(HTTP 503: MMseqs2 search will be back!\)"):
        remote.mmseqs_search(DATA / "luxC_reference.faa", log=lambda *_: None)


def test_unreachable_host_is_a_clear_error(monkeypatch):
    import requests

    from famscan import remote

    def boom(*a, **k):
        raise requests.ConnectionError("no route")

    monkeypatch.setattr(remote.requests, "request", boom)
    for call in (lambda: remote.mmseqs_search(__file__, timeout=1),
                 lambda: remote.ncbi_blast("MKV", timeout=1, log=lambda *_: None)):
        try:
            call()
        except remote.ServiceUnavailable as e:
            assert "could not reach" in str(e) and "fetchdb" in str(e)
        else:
            raise AssertionError("expected ServiceUnavailable")


def _down(*_a, **_k):
    raise RuntimeError("boom")


def test_service_outage_continues_and_is_reported(tmp_path, monkeypatch):
    """A structure-server outage must not sink the whole run."""
    from famscan import structure
    from famscan.workflows import run_workflow

    monkeypatch.setattr(structure, "get_structure", lambda *a, **k: (None, "unreachable"))
    ctx = run_workflow("family", {"input": str(DATA / "proteome_938293.faa"),
                                  "query": str(DATA / "PF02826.hmm"), "structure": "esmfold",
                                  "top_n": 2}, tmp_path, log=lambda *_: None)
    assert len(ctx.focus) == 7                      # every earlier stage survived
    assert ctx.outages and ctx.outages[0]["stage"] == "structure"
    assert "unavailable" in (tmp_path / "report.md").read_text()


def test_service_outage_can_ask_and_stop(tmp_path, monkeypatch):
    from famscan import structure
    from famscan.workflows import run_workflow

    monkeypatch.setattr(structure, "get_structure", lambda *a, **k: (None, "unreachable"))
    asked = []
    ctx = run_workflow("family", {"input": str(DATA / "proteome_938293.faa"),
                                  "query": str(DATA / "PF02826.hmm"), "structure": "esmfold",
                                  "top_n": 2, "on_service_down": "ask"}, tmp_path,
                       log=lambda *_: None, ask=lambda q: (asked.append(q), False)[1])
    assert asked and "Continue without it?" in asked[0]
    assert "Stopped early" in (tmp_path / "report.md").read_text()
    assert ctx.outages


def test_email_config(tmp_path, monkeypatch):
    from famscan import config
    monkeypatch.setenv("FAMSCAN_CONFIG", str(tmp_path / "c.json"))
    monkeypatch.delenv("FAMSCAN_EMAIL", raising=False)
    assert config.email() is None
    config.save(email="a@b.org")
    assert config.email() == "a@b.org"
    assert config.email("other@c.org") == "other@c.org"      # explicit wins
    assert not config.valid_email("nope")
