#!/usr/bin/env bash
# Runs every workflow on the bundled demo data (offline; structure steps need internet).
# Stand-ins: data/mini_domain_library.hmm replaces Pfam-A.hmm, the Bacillus proteins
# replace the human proteome, and small FASTAs replace CARD / toxin / allergen dbs,
# so the screening demos show the mechanics, not real biological findings.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONWARNINGS=ignore
D=data; O=out
S="--set structure=none"

famscan workflow family --out $O/01_family --set input=$D/proteome_938293.faa \
  --set query=$D/PF02826.hmm --set hmm_db=$D/mini_domain_library.hmm --set top_n=5
famscan workflow enzyme --out $O/02_enzyme --set input=$D/bacillus_contig.fna \
  --set query=$D/Thioesterase.hmm --set hmm_db=$D/mini_domain_library.hmm \
  --set "motifs=thioesterase_nucleophile=G-x-S-x-G" $S
famscan workflow annotate --out $O/03_annotate --set input=$D/bacillus_contig.fna \
  --set hmm_db=$D/mini_domain_library.hmm --set annotation_db=$D/reference_proteins.faa
famscan workflow drug_target --out $O/04_drug_target --set input=$D/proteome_938293.faa \
  --set host=$D/bacillus_proteins.faa --set hmm_db=$D/mini_domain_library.hmm $S
famscan workflow resistance --out $O/05_resistance --set input=$D/bacillus_contig.fna \
  --set resistance_db=$D/standin_screen_db.faa
famscan workflow toxin_allergen --out $O/06_toxin_allergen --set input=$D/proteome_938293.faa \
  --set toxin_db=$D/standin_screen_db.faa --set allergen_db=$D/standin_screen_db.faa
famscan workflow evolution --out $O/07_evolution --set query=$D/PF02826.hmm \
  --set "proteomes=$D/proteome_938293.faa $D/bacillus_contig.fna $D/cp040672_cds.fna"
famscan workflow unknown --out $O/08_unknown --set input=$D/unknown_proteins.faa \
  --set hmm_db=$D/mini_domain_library.hmm $S
famscan workflow lab_candidates --out $O/09_lab_candidates --set input=$O/01_family/final_candidates.faa $S
famscan workflow engineering --out $O/10_engineering --set reference=$D/luxC_reference.faa \
  --set input=$D/luxC_homologs.faa
famscan workflow strains --out $O/11_strains \
  --set "proteomes=$D/bacillus_contig.fna $D/bacillus_contig_variant.fna" --set hmm_db=$D/mini_domain_library.hmm
for p in blastp blastn blastx tblastn tblastx; do
  case $p in
    blastp)  q=$D/luxC_reference.faa;  t=$D/luxC_homologs.faa ;;
    blastn)  q=$D/bmyD_gene.fna;       t=$D/bacillus_contig.fna ;;
    blastx)  q=$D/bmyD_gene.fna;       t=$D/bacillus_proteins.faa ;;
    tblastn) q=$D/bmyD_protein.faa;    t=$D/bacillus_contig.fna ;;
    tblastx) q=$D/bmyD_gene.fna;       t=$D/bacillus_contig.fna ;;
  esac
  famscan search -q $q -t $t -p $p --out $O/12_search_$p
done
