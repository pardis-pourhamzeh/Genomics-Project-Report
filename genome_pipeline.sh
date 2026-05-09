#!/bin/bash
# ============================================================================
# BCG 2026 Genomics — final project pipeline (single-file equivalent)
# Author : Babi Pourhamzeh (BCG2026_Pourhamzeh_P)
#
# Paired-end exome trio analysis on chromosome 20 (GRCh38). Five trios are
# analysed under three inheritance models: AR (trio_1, trio_4, trio_5),
# AD inherited (trio_2 — affected parent inferred empirically), AD de novo
# (trio_3). The pipeline mirrors the structure of the BCG 2025 reference
# script (Caccia & Gautieri, single-end chr16) and is adapted for paired-end
# data, GRCh38, an updated VEP cache and a stricter filter_vep step.
#
# Run this from ~/exam_project_scripts/ on the leon server. All output goes
# to ~/exam_project/, organised per trio.
# ============================================================================
set -euo pipefail

# ---------- TRIOS / INHERITANCE -------------------------------------------
trios=("trio_1" "trio_2" "trio_3" "trio_4" "trio_5")
declare -A MODE
MODE[trio_1]="AR"
MODE[trio_2]="AD_inh"
MODE[trio_3]="AD_dn"
MODE[trio_4]="AR"
MODE[trio_5]="AR"

# Per-trio individual IDs (taken from FASTQ filenames in the personal folder)
declare -A CHILD FATHER MOTHER
for t in trio_1 trio_2 trio_3 trio_4; do
  CHILD[$t]="HG00427"; FATHER[$t]="HG00428"; MOTHER[$t]="HG00429"
done
CHILD[trio_5]="HG00421"; FATHER[trio_5]="HG00422"; MOTHER[trio_5]="HG00423"

# ---------- PATHS ---------------------------------------------------------
USER_FOLDER="BCG2026_Pourhamzeh_P"
EXAM_DIR="/home/BCG2026_exam"
SRC="${EXAM_DIR}/${USER_FOLDER}"
WORKDIR="${HOME}/exam_project"

REFERENCE="${EXAM_DIR}/chr20.fa"
BOWTIE_INDEX="${EXAM_DIR}/chr20"
TARGET_BED="${EXAM_DIR}/chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed"
DISORDERS="${EXAM_DIR}/list_disorders.txt"

VEP="/opt/vep/ensembl-vep/vep"
FILTER_VEP="/opt/vep/ensembl-vep/filter_vep"
VEP_CACHE="/home/data/vep_cache"
CLINVAR="/data/vep_annotations/clinvar.vcf.gz"

# Local writable copy of chr20.fa — VEP needs to write a .vep.lock file next
# to its FASTA, and the cache directory is read-only for us.
LOCAL_FASTA="${WORKDIR}/chr20.fa"

THREADS=6

# ---------- 0.  WORK DIR + SYMBOLIC LINKS ---------------------------------
echo "[0/9] Setting up work directory ${WORKDIR}"
mkdir -p "${WORKDIR}"
[ -e "${LOCAL_FASTA}" ]      || cp "${REFERENCE}"        "${LOCAL_FASTA}"
[ -e "${LOCAL_FASTA}.fai" ]  || cp "${REFERENCE}.fai"    "${LOCAL_FASTA}.fai" 2>/dev/null \
                              || samtools faidx "${LOCAL_FASTA}"

for t in "${trios[@]}"; do
  src="${SRC}/${t}"
  dst="${WORKDIR}/${t}"
  mkdir -p "${dst}"
  for fq in "${src}"/*.fq.gz; do
    bn=$(basename "${fq}")
    [ -e "${dst}/${bn}" ] || ln -s "${fq}" "${dst}/${bn}"
  done
done

# ---------- HELPER --------------------------------------------------------
filter_by_inheritance () {
    # awk-based filter, independent of VCF column order
    local vcf_in="$1" vcf_out="$2" mode="$3" affected="${4:-}"
    awk -v mode="${mode}" -v affected="${affected}" '
      BEGIN { OFS="\t" }
      /^##/ { print; next }
      /^#CHROM/ {
        for (i=10; i<=NF; i++) {
          if      ($i=="child")  c=i
          else if ($i=="father") f=i
          else if ($i=="mother") m=i
        }
        print; next
      }
      {
        split($c, ca, ":"); cg=ca[1]
        split($f, fa, ":"); fg=fa[1]
        split($m, ma, ":"); mg=ma[1]
        keep=0
        if      (mode=="AR"     && cg=="1/1" && fg=="0/1" && mg=="0/1") keep=1
        else if (mode=="AD_dn"  && cg=="0/1" && fg=="0/0" && mg=="0/0") keep=1
        else if (mode=="AD_inh" && cg=="0/1") {
          if (affected=="father" && fg=="0/1" && mg=="0/0") keep=1
          if (affected=="mother" && fg=="0/0" && mg=="0/1") keep=1
        }
        if (keep) print
      }
    ' "${vcf_in}" > "${vcf_out}"
}

run_vep () {
    local in_vcf="$1" out_vcf="$2" log="$3"
    [ -s "${in_vcf}" ] || { echo "    skipping ${in_vcf} (empty)"; return; }
    "${VEP}" \
        -i "${in_vcf}" \
        -o "${out_vcf}" \
        --vcf --cache --offline \
        --assembly GRCh38 \
        --dir_cache "${VEP_CACHE}" \
        --fasta "${LOCAL_FASTA}" \
        --use_given_ref \
        --everything \
        --fork "${THREADS}" \
        --custom "file=${CLINVAR},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNDN" \
        --stats_file "${out_vcf%.vcf}_summary.html" \
        --force_overwrite \
        2> "${log}"
}

apply_filter_vep () {
    local in_vcf="$1" out_high="$2" out_mod="$3"
    [ -s "${in_vcf}" ] || { echo "    skipping ${in_vcf} (empty)"; return; }
    "${FILTER_VEP}" -i "${in_vcf}" -o "${out_high}" --force_overwrite \
        --filter "IMPACT is HIGH and (not MAX_AF or MAX_AF < 0.0001)"
    "${FILTER_VEP}" -i "${in_vcf}" -o "${out_mod}"  --force_overwrite \
        --filter "IMPACT is MODERATE and (not MAX_AF or MAX_AF < 0.0001)"
}

# ---------- MAIN PER-TRIO LOOP --------------------------------------------
for t in "${trios[@]}"; do
  inh="${MODE[$t]}"
  C="${CHILD[$t]}"; F="${FATHER[$t]}"; M="${MOTHER[$t]}"

  echo
  echo "============================================================"
  echo "  ${t}  (${inh})  —  child=${C}  father=${F}  mother=${M}"
  echo "============================================================"

  cd "${WORKDIR}/${t}"
  mkdir -p fastqc fastqc/qualimap ucsc_tracks logs

  # ---- 1. FastQC ----------------------------------------------------------
  echo "[1/9] FastQC"
  fastqc -t ${THREADS} -o fastqc/ \
         ${C}.targets_R1.fq.gz ${C}.targets_R2.fq.gz \
         ${F}.targets_R1.fq.gz ${F}.targets_R2.fq.gz \
         ${M}.targets_R1.fq.gz ${M}.targets_R2.fq.gz \
         2> logs/fastqc.log

  # ---- 2. Bowtie2 paired-end alignment ------------------------------------
  echo "[2/9] Bowtie2 (paired-end) + samtools sort + index"
  for r in child:${C} father:${F} mother:${M}; do
    role="${r%%:*}"; id="${r##*:}"
    bowtie2 -x "${BOWTIE_INDEX}" \
            -1 "${id}.targets_R1.fq.gz" \
            -2 "${id}.targets_R2.fq.gz" \
            --rg-id "${role}" --rg "SM:${role}" \
            -p ${THREADS} \
            2> logs/bowtie2_${role}.log \
      | samtools view -Sb - \
      | samtools sort -@ ${THREADS} -o ${role}_sorted.bam -
    samtools index ${role}_sorted.bam
  done

  # ---- 3. Qualimap --------------------------------------------------------
  echo "[3/9] Qualimap"
  for role in child father mother; do
    qualimap bamqc \
        -bam ${role}_sorted.bam \
        --feature-file "${TARGET_BED}" \
        -outdir fastqc/qualimap/${role} \
        --java-mem-size=2G -nt ${THREADS} \
        2>/dev/null
  done

  # ---- 4. FreeBayes joint trio calling ------------------------------------
  echo "[4/9] FreeBayes joint variant calling"
  freebayes -f "${REFERENCE}" \
            -m 20 -C 5 -Q 10 -q 10 --min-coverage 20 \
            child_sorted.bam father_sorted.bam mother_sorted.bam \
        > ${t}.vcf 2> logs/freebayes.log

  # ---- 5. Normalise + GT/QUAL filter --------------------------------------
  echo "[5/9] bcftools norm + sample sort + QUAL>30 + drop missing GTs"
  bcftools norm -m -any -f "${REFERENCE}" ${t}.vcf > ${t}_norm.vcf
  bcftools query -l ${t}_norm.vcf | sort > samples.txt   # alphabetical: child, father, mother
  bcftools view -S samples.txt ${t}_norm.vcf \
    | bcftools filter -i 'QUAL>30' \
    | bcftools filter -e 'GT[*]~"\."' \
    > ${t}_sorted_filtered.vcf

  # ---- 6. Inheritance-pattern filter + exome-target intersect -------------
  echo "[6/9] Inheritance-based GT filter + exome-target intersect"
  if [ "${inh}" = "AD_inh" ]; then
    # mode_inheritance.tsv was missing — produce both candidate sets and
    # let downstream interpretation decide which parent is affected.
    filter_by_inheritance ${t}_sorted_filtered.vcf candidates_${t}_fatherAff.vcf AD_inh father
    filter_by_inheritance ${t}_sorted_filtered.vcf candidates_${t}_motherAff.vcf AD_inh mother
    for tag in fatherAff motherAff; do
      grep "^#" candidates_${t}_${tag}.vcf > candidates_${t}_${tag}_final.vcf
      bedtools intersect -a candidates_${t}_${tag}.vcf -b "${TARGET_BED}" -u \
        | grep -v '^#' >> candidates_${t}_${tag}_final.vcf || true
    done
  else
    filter_by_inheritance ${t}_sorted_filtered.vcf candidates_${t}.vcf "${inh}"
    grep "^#" candidates_${t}.vcf > candidates_${t}_final.vcf
    bedtools intersect -a candidates_${t}.vcf -b "${TARGET_BED}" -u \
      | grep -v '^#' >> candidates_${t}_final.vcf || true
  fi

  # ---- 7. VEP annotation + filter_vep -------------------------------------
  # NOTE: --pick_allele was deliberately omitted. With it, the trio_3 SAMHD1
  # frameshift was rewritten as a TLDC2 downstream MODIFIER and lost at the
  # IMPACT step; --everything keeps all transcripts.
  echo "[7/9] VEP + filter_vep (IMPACT==HIGH and MAX_AF<1e-4)"
  if [ "${inh}" = "AD_inh" ]; then
    for tag in fatherAff motherAff; do
      run_vep candidates_${t}_${tag}_final.vcf candidates_${t}_${tag}_vep.vcf logs/vep_${t}_${tag}.log
      apply_filter_vep candidates_${t}_${tag}_vep.vcf \
                       candidates_${t}_${tag}_HIGH_rare.vcf \
                       candidates_${t}_${tag}_MODERATE_rare.vcf
    done
  else
    run_vep candidates_${t}_final.vcf candidates_${t}_vep.vcf logs/vep_${t}.log
    apply_filter_vep candidates_${t}_vep.vcf \
                     candidates_${t}_HIGH_rare.vcf \
                     candidates_${t}_MODERATE_rare.vcf
  fi

  # ---- 8. UCSC coverage tracks --------------------------------------------
  echo "[8/9] Bedgraph coverage tracks for UCSC"
  for role in child father mother; do
    bedtools genomecov -ibam ${role}_sorted.bam \
        -bg -trackline -trackopts "name=\"${t}_${role}\"" -max 100 \
      > ucsc_tracks/${role}.bg
  done

  echo "Completed ${t}"
done

# ---------- 9. AGGREGATED MULTIQC -----------------------------------------
echo
echo "[9/9] MultiQC aggregation across all trios"
mkdir -p "${WORKDIR}/multiqc"
cd "${WORKDIR}"
multiqc -f -o multiqc trio_*/fastqc trio_*/fastqc/qualimap 2> multiqc/multiqc.log

echo
echo "============================================================"
echo "Pipeline COMPLETED for all five trios."
echo "  candidate VCFs : ${WORKDIR}/trio_*/candidates_*_HIGH_rare.vcf"
echo "  MultiQC report : ${WORKDIR}/multiqc/multiqc_report.html"
echo "  bedgraphs      : ${WORKDIR}/trio_*/ucsc_tracks/*.bg"
echo "============================================================"
