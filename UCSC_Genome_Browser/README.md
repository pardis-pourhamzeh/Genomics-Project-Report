# UCSC Genome Browser views

Per-trio visualisation of the candidate locus on the UCSC Genome Browser (GRCh38/hg38), with three custom bedgraph tracks (child / father / mother) generated from `bedtools genomecov`.

The bedgraphs were subset to a 300 kb window centred on the candidate variant by `scripts/ucsc_subset.sh`, and uploaded as a single concatenated track file per trio.

> **trio_1 is intentionally absent.** It carries no qualifying candidate variant after `filter_vep` on either the HIGH-rare or MODERATE-rare lists and is reported as healthy in the project report. There is therefore no specific locus to visualise. Coverage uniformity for trio_1 can be verified in `../MultiQC_reports/multiqc_trio_1.html`.

## Files

| File | Trio | Locus | Notes |
|---|---|---|---|
| `trio_2_CYP24A1.png` | 2 | chr20:54,107,411–54,207,411 | AD-inherited — variant in child + father, absent from mother |
| `trio_3_SAMHD1.png` | 3 | chr20:36,848,463–36,948,463 | AD-de novo — variant only in child |
| `trio_4_MKKS.png` | 4 | chr20:10,357,697–10,457,697 | AR — child homozygous, parents heterozygous |
| `trio_5_ADA.png` | 5 | chr20:44,570,357–44,670,357 | AR — child homozygous, parents heterozygous |

For each figure, all three tracks have comparable on-target coverage (≈ 30–60×), confirming that the variant pattern reflects a true genotype and not a coverage drop-out.

## How to regenerate

```bash
# 1. subset the bedgraphs to a 100 kb window


# 2. pull the small files


# 3. in browser: https://genome.ucsc.edu/cgi-bin/hgGateway
#    - Genome → GRCh38/hg38
#    - My Data → Custom Tracks → Add custom tracks → upload trio_<i>_3tracks.region.bg
#    - Position → search for the variant coordinate from the diagnosis table
#    - Screenshot → save here as trio_<i>_<gene>.png
```
