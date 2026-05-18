# IGV read-pileup views

Per-trio Integrative Genomics Viewer (IGV) screenshots of the candidate locus, using the trio's child / father / mother sorted BAMs aligned to GRCh38. IGV provides an independent visual confirmation of the variant — the proband shows alt-allele evidence at the candidate position while the unaffected family members are reference-only.

> **trio_1 is intentionally absent.** No qualifying candidate variant was retained for trio_1 after filter_vep on the HIGH- or MODERATE-rare short-lists; there is no locus to pile up.

## Files

| File | Trio | Variant | What it shows |
|---|---|---|---|
| `trio_2_CYP24A1.png` | 2 | chr20:54,157,411 G>GA | 1 bp insertion in the child + father; mother homozygous reference |
| `trio_3_SAMHD1.png` | 3 | chr20:36,898,463 CT>C | 1 bp deletion only in the child (de novo); both parents reference |
| `trio_4_MKKS.png` | 4 | chr20:10,407,697 CA>C | child homozygous deletion; parents heterozygous |
| `trio_5_ADA.png` | 5 | chr20:44,620,357 CTT>C | child homozygous deletion; parents heterozygous |

## How to regenerate

```bash
1. pull the BAMs and indexes for the trios intended to visualise

# 2. open IGV (https://igv.org/app/ in the browser, or the desktop client)
#    - Genome → Human (GRCh38/hg38)
#    - Tracks → Local file → load all three BAMs (.bam) for one trio
#    - Search for the variant coordinate (e.g. chr20:36898463 for trio_3 SAMHD1)
#    - Take a screenshot and save it here
```

## Web VEP validation

Per the professor's recommendation, every candidate was also re-checked on the online Ensembl VEP at <https://www.ensembl.org/Tools/VEP> by pasting the one-line VCF row of the candidate. Online VEP confirmed the gene, the consequence, the ClinVar phenotype, and the gnomAD frequency (≈ 0) for all four candidates. Where applicable, the online VEP report PDF is also stored here.
