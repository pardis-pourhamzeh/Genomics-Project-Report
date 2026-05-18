# Clinical Diagnosis of Rare Mendelian Disorders by Trio-Based WES on Chromosome 20

Final project for the **Genomics course** at the Università degli Studi di Milano Statale.
A trio-based whole-exome sequencing (WES) analysis is performed on five simulated trios (chr20, GRCh38) to identify candidate disease-causing variants under prescribed inheritance models, validate them against ClinVar and online VEP, and report a final diagnosis for each trio.

> **Author:** Pardis Pourhamzeh
> **Course:** BCG 2026 — Genomics, Università degli Studi di Milano Statale
> **Final report (PDF):** [`Pourhamzeh_BCG2026_Genomics_report.pdf`](Pourhamzeh_BCG2026_Genomics_report.pdf)

---

## 1. Final diagnosis

| Trio | Mode | Gene | Position (chr20, GRCh38) | REF → ALT | Impact | Phenotype | Status |
|:----:|:----|:-----|:------------------------:|:---------:|:------:|:----------|:-------|
| 1 | AR | — | — | — | — | — | **Healthy** |
| 2 | AD-inh (father) | **CYP24A1** | 54,157,411 | G → GA | frameshift (HIGH) | Hypercalcemia, infantile, type 1 | Affected |
| 3 | AD-de novo | **SAMHD1** | 36,898,463 | CT → C | frameshift (HIGH) | Aicardi-Goutières syndrome 5 | Affected |
| 4 | AR | **MKKS** | 10,407,697 | CA → C | frameshift (HIGH) | Bardet-Biedl / McKusick-Kaufman (BBS6) | Affected |
| 5 | AR | **ADA** | 44,620,357 | CTT → C | frameshift (HIGH) | ADA-SCID | Affected |

All four candidate variants are flagged Pathogenic in ClinVar, are within an exome target region, and segregate with the inheritance model prescribed for the trio.

---

## 2. Repository layout

```
.
├── README.md                                  # this file
├── Pourhamzeh_BCG2026_Genomics_report.pdf     # final 6-page report
├── genome_pipeline.sh                         # single-file pipeline (read this first)
├── VEP_output_visualization.md                # how to read the VEP / filter_vep output
├── VCF_trios/                                 # candidate VCFs per trio (HIGH-rare short-list)
├── UCSC_Genome_Browser/                       # UCSC screenshots
├── IGV_genomes_view/                          # IGV screenshots
└── MultiQC_reports/                           # MultiQC HTML
```

**Two equivalent ways to run the analysis** are provided. `genome_pipeline.sh` at the repository root is a single-file consolidation of the entire workflow.

---

## 3. Data and environment

- **Input:** simulated paired-end FASTQ files (R1 + R2; 150 bp) for five trios in `/home/BCG2026_exam/.../trio_{1..5}/`.
- **Reference:** chromosome 20 of GRCh38 (`/home/BCG2026_exam/chr20.fa`) with a Bowtie2 index (`chr20.{1..4}.bt2` and `chr20.rev.{1,2}.bt2`).
- **Target panel:** `/home/BCG2026_exam/chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed`.
- **Disorder list:** `/home/BCG2026_exam/list_disorders.txt` (AD and AR partitions).
- **Trio individuals:** trio_1–4 share genetic background HG00427 (child) / HG00428 (father) / HG00429 (mother); trio_5 is HG00421 / HG00422 / HG00423.
- **Inheritance models** (per project instructions): AR for trios 1, 4, 5; AD inherited for trio 2; AD de novo for trio 3.


### Software

| Tool | Version | Use |
|---|---|---|
| FastQC | 0.12 | per-read quality control |
| MultiQC | latest | report aggregation |
| Bowtie2 | 2.5 | paired-end alignment |
| Samtools | 1.x | SAM → BAM, sort, index |
| Qualimap | 2.x | per-BAM target QC |
| FreeBayes | 1.x | joint variant calling per trio |
| bcftools | 1.x | normalisation + QUAL/missing-GT filters |
| BEDTools | 2.x | exome-target intersect, coverage tracks |
| Ensembl VEP | 115 | functional annotation (offline cache + ClinVar) |
| filter_vep | 115 | IMPACT + MAX_AF short-list |
| UCSC Genome Browser | — | per-trio coverage visualisation |
| IGV | 2.x | per-trio read pileup validation |

---

## 4. Reproducible pipeline

The full pipeline takes ≈ 2–4 hours, mostly Bowtie2. You can run it either way.



### Key methodological notes

- **Inheritance grep order.** After `bcftools view -S samples.txt` (alphabetical sort) the sample column order is **child / father / mother**. Pattern `1/1.*0/1.*0/1` is therefore AR; `0/1.*0/0.*0/0` is AD-de novo; `0/1.*0/1.*0/0` and `0/1.*0/0.*0/1` are the two AD-inherited possibilities.
- **Missing genotypes.** `./.` indicates low coverage, *not* homozygous reference. The pipeline drops any record carrying `./.` in any sample (`bcftools filter -e 'GT[*]~"\."'`) before applying the inheritance pattern.
- **VEP `--pick_allele` was deliberately removed.** With `--pick_allele`, the trio_3 SAMHD1 frameshift was overwritten by an overlapping TLDC2 downstream-MODIFIER consequence and the variant fell out at the IMPACT filter. Removing `--pick_allele` keeps consequences across all transcripts and the SAMHD1 HIGH call survives.
- **Indel representation.** The trio_5 ADA candidate was initially called by FreeBayes as a longer compound indel (`CTTTTCA → CTTCA`) which `bcftools norm` collapsed to `CTT → C`. Online VEP confirmed the resulting frameshift and the ADA-SCID phenotype. Indel candidates were always re-validated on online VEP.
- **VEP cache FASTA workaround.** The bundled FASTA in `/home/data/vep_cache/...` is read-only, so VEP cannot write its `.vep.lock` file there. We pass our own writable copy of `chr20.fa` via `--fasta`.

---

## 5. Visual evidence

For each affected trio we provide two independent visualisations of the candidate locus:

1. **UCSC Genome Browser** screenshots (`UCSC_Genome_Browser/`) — coverage tracks for child / father / mother around the variant, generated from `bedtools genomecov`-derived bedgraphs.
2. **IGV** screenshots (`IGV_genomes_view/`) — read pileups around the variant, with alt-allele evidence visible in the proband and reference-only reads in the unaffected family members.

The pair of views is what allows us to conclude that the variant is real, fully covered, and segregates with the prescribed inheritance — the same logic applied in the reference report from the previous edition of the course.

---

## 6. Candidate VCFs

`VCF_trios/` contains the `candidates_<i>_HIGH_rare.vcf` short-list for every trio after `filter_vep`. Each file has the full VEP `CSQ` field, MAX_AF, SIFT, PolyPhen, and ClinVar `CLNSIG` / `CLNDN` columns. Every candidate listed in the diagnosis table comes from these files unchanged.

For trio 2 we keep both `candidates_2_fatherAff_HIGH_rare.vcf` and `candidates_2_motherAff_HIGH_rare.vcf` as evidence that only the father-affected hypothesis yields a qualifying variant.

---

## 7. Validation strategy

- **ClinVar:** every candidate is annotated as Pathogenic for the matched phenotype.
- **list_disorders.txt:** every candidate gene maps to a disorder in the project's disorder list under the correct inheritance partition.
- **Online VEP:** every candidate was re-checked at <https://www.ensembl.org/Tools/VEP> to validate the consequence and phenotype, especially for indels.
- **Allele frequency:** every candidate has `MAX_AF` close to or equal to 0 in gnomAD/1000G — consistent with rare-disease alleles.
- **Inheritance segregation:** every candidate's genotype in child / father / mother matches the prescribed model.
- **Coverage:** every candidate is supported by ≥10 reads in all three family members at the variant site (Qualimap mean ≈ 40× per sample; ≥95% of the target panel ≥30×).

---

## 8. Limitations

The analysis is restricted to a single chromosome and a fixed exome target panel; X-linked, mitochondrial and compound-heterozygous patterns are not searched; the proband phenotype is not given so candidate-gene prioritisation is purely genotype-driven; and the simulated nature of the data means certain quirks (e.g. the trio_5 indel representation) reflect properties of the simulator rather than of a clinical sequencing assay.

