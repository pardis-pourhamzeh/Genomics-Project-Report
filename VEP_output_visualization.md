# Variant Effect Predictor Results

Links to view the visual results from the online Ensembl Variant Effect Predictor
(<https://www.ensembl.org/Tools/VEP>, GRCh38 assembly) for each candidate variant
identified in this project. Each link opens the per-variant report that was used
to validate the gene, the consequence, the ClinVar phenotype and the gnomAD
allele frequency — the step the project guidance recommended for all indel
candidates.

- [trio_1 (healthy — no candidate)](#)
- [trio_2 — CYP24A1 frameshift (Hypercalcemia, infantile, 1)](https://www.ensembl.org/Homo_sapiens/Tools/VEP/Results?tl=QqZSTWfEAKOVdo5G-11886343)
- [trio_3 — SAMHD1 frameshift (Aicardi-Goutières syndrome 5)](https://www.ensembl.org/Homo_sapiens/Tools/VEP/Results?tl=00AeZ2HEDyAjDzX9-11886365)
- [trio_4 — MKKS frameshift (Bardet-Biedl / McKusick-Kaufman, BBS6)](https://www.ensembl.org/Homo_sapiens/Tools/VEP/Results?tl=bdnK7mP9Vzy4o5fS-11886369)
- [trio_5 — ADA frameshift (ADA-SCID)](https://www.ensembl.org/Homo_sapiens/Tools/VEP/Results?tl=ICjlJido3L6VZESp-11886377)

All four affected-trio candidates have been re-checked on the online VEP and the
gene, consequence, ClinVar phenotype and gnomAD allele frequency match the
diagnosis reported in Table 1 of the report. The Ensembl VEP result pages are
persistent for several days and can be re-generated from the one-line VCFs
below.

## One-line VCFs (input to online VEP)

```
chr20  54157411  .  G    GA    # trio_2  CYP24A1  Hypercalcemia, infantile, 1
chr20  36898463  .  CT   C     # trio_3  SAMHD1   Aicardi-Goutières syndrome 5
chr20  10407697  .  CA   C     # trio_4  MKKS     Bardet-Biedl / MKKS (BBS6)
chr20  44620357  .  CTT  C     # trio_5  ADA      ADA-SCID
```

Trio 1 is reported as healthy (no qualifying candidate after `filter_vep` on
either the HIGH-rare or MODERATE-rare short-lists) and therefore has no online
VEP entry.
