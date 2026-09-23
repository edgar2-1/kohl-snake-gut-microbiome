#!/usr/bin/env bash
# =============================================================================
# Kohl / N. rhombifer neonate gut microbiome -- 16S rRNA pipeline (QIIME 2)
# =============================================================================
# Activate the QIIME 2 amplicon env, then:
#   bash run_16s.sh 2>&1 | tee run_16s.log
#
# Marker:    16S V4 (515F/806R, Earth Microbiome Project)
# Primers:   515F = GTGYCAGCMGCCGCGGTAA   (in-read, MUST be trimmed)
#            806R = GGACTACNVGGGTWTCTAAT  (in-read, MUST be trimmed)
# Reads:     2 x 251 bp PE MiSeq
# Samples:   32 N. rhombifer LI swabs (Newborn=3, Sterile=16, Muck=10, Fetus=3)
#            + 4 Blank PCR negative controls (used by decontam)
# =============================================================================

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/root/kohl_v1}"
META="$PROJECT_DIR/meta/metadata.tsv"
MANIFEST="$PROJECT_DIR/meta/manifest_16s.tsv"
OUT="$PROJECT_DIR/qiime/16s"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"

FWD_PRIMER="GTGYCAGCMGCCGCGGTAA"
REV_PRIMER="GGACTACNVGGGTWTCTAAT"

# Truncation lengths for DADA2.
# V4 amplicon ~252 bp; after primer trimming reads are ~230 bp.
# F=220, R=180 -> 220+180-252 = 148 bp overlap (>= 12 needed).
TRUNC_LEN_F=220
TRUNC_LEN_R=180

# Sampling depth for core-diversity. Tune after viewing alpha-rarefaction.qzv.
SAMPLING_DEPTH=5000

# SILVA 138.2 (symlinked from inornatus_v1)
SILVA_SEQS="$PROJECT_DIR/qiime/reference/silva-138-seqs-derep.qza"
SILVA_TAXA="$PROJECT_DIR/qiime/reference/silva-138-taxa-derep.qza"

mkdir -p "$OUT" "$OUT/figures"
cd "$OUT"

echo "==> 16S pipeline starting at $(date)"
echo "    Threads: $THREADS  Project: $PROJECT_DIR"

# 1. Import paired-end reads
if [[ ! -f demux.qza ]]; then
  echo "==> [1/12] Importing PE fastqs"
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --output-path demux.qza \
    --input-format PairedEndFastqManifestPhred33V2
  qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
fi

# 2. Cutadapt: trim 515F/806R primers (present in-read)
if [[ ! -f demux-trimmed.qza ]]; then
  echo "==> [2/12] Cutadapt: trimming 515F/806R primers"
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux.qza \
    --p-front-f "$FWD_PRIMER" \
    --p-front-r "$REV_PRIMER" \
    --p-discard-untrimmed \
    --p-cores "$THREADS" \
    --o-trimmed-sequences demux-trimmed.qza \
    --verbose
  qiime demux summarize --i-data demux-trimmed.qza --o-visualization demux-trimmed.qzv
fi

# 3. DADA2 denoise -> ASV table + rep-seqs
if [[ ! -f table.qza ]]; then
  echo "==> [3/12] DADA2 denoise-paired (trunc-f=$TRUNC_LEN_F trunc-r=$TRUNC_LEN_R)"
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux-trimmed.qza \
    --p-trunc-len-f "$TRUNC_LEN_F" \
    --p-trunc-len-r "$TRUNC_LEN_R" \
    --p-n-threads "$THREADS" \
    --o-representative-sequences rep-seqs.qza \
    --o-table table.qza \
    --o-denoising-stats denoising-stats.qza \
    --o-base-transition-stats base-transition-stats.qza \
    --verbose
  qiime feature-table summarize \
    --i-table table.qza --o-summary table.qzv \
    --o-feature-frequencies table-feature-frequencies.qza \
    --o-sample-frequencies table-sample-frequencies.qza \
    --m-metadata-file "$META"
  qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
  qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv
fi

# 4. Phylogeny (for UniFrac / Faith's PD)
if [[ ! -f rooted-tree.qza ]]; then
  echo "==> [4/12] Align + mask + FastTree"
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza --p-n-threads "$THREADS" \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza
fi

# 5. Taxonomy (SILVA 138.2 via vsearch consensus)
if [[ ! -f taxonomy.qza ]]; then
  echo "==> [5/12] Classifying with SILVA (classify-consensus-vsearch)"
  qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seqs.qza \
    --i-reference-reads "$SILVA_SEQS" \
    --i-reference-taxonomy "$SILVA_TAXA" \
    --p-threads "$THREADS" \
    --p-perc-identity 0.80 \
    --p-maxaccepts 10 \
    --o-classification taxonomy.qza \
    --o-search-results silva-search-results.qza
  qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
fi

# 6. Filter host/non-bacterial, then decontam vs Blank negatives
if [[ ! -f table-clean.qza ]]; then
  echo "==> [6a/12] Filter mitochondria/chloroplast/eukaryota, require p__"
  qiime taxa filter-table \
    --i-table table.qza --i-taxonomy taxonomy.qza \
    --p-include p__ --p-exclude mitochondria,chloroplast,eukaryota \
    --o-filtered-table table-noHost.qza
  qiime taxa filter-seqs \
    --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza \
    --p-include p__ --p-exclude mitochondria,chloroplast,eukaryota \
    --o-filtered-sequences rep-seqs-noHost.qza

  echo "==> [6b/12] decontam prevalence filter vs Blank1-4"
  bash "$PROJECT_DIR/qiime/decontam_filter.sh" \
       "$OUT/table-noHost.qza" "$OUT/rep-seqs-noHost.qza" "$META" "$OUT"

  mv table-decontam.qza    table-clean.qza
  mv rep-seqs-decontam.qza rep-seqs-clean.qza
  qiime feature-table summarize \
    --i-table table-clean.qza --o-summary table-clean.qzv \
    --o-feature-frequencies table-clean-feature-frequencies.qza \
    --o-sample-frequencies table-clean-sample-frequencies.qza \
    --m-metadata-file "$META"
fi

# 6c. Drop the Blank samples from the clean table -- they were only needed for
#     the decontam step. Downstream analyses operate on biological samples only.
if [[ ! -f table-bio.qza ]]; then
  echo "==> [6c/12] Dropping Blank samples from biological table"
  qiime feature-table filter-samples \
    --i-table table-clean.qza \
    --m-metadata-file "$META" \
    --p-where "[is-negative-control]='false'" \
    --o-filtered-table table-bio.qza
  qiime feature-table summarize \
    --i-table table-bio.qza --o-summary table-bio.qzv \
    --m-metadata-file "$META" \
    --o-feature-frequencies table-bio-feature-frequencies.qza \
    --o-sample-frequencies table-bio-sample-frequencies.qza
fi

# 7. Alpha rarefaction (choose sampling depth from this)
if [[ ! -f alpha-rarefaction.qzv ]]; then
  echo "==> [7/12] Alpha rarefaction"
  qiime diversity alpha-rarefaction \
    --i-table table-bio.qza --i-phylogeny rooted-tree.qza \
    --p-max-depth 20000 --m-metadata-file "$META" \
    --o-visualization alpha-rarefaction.qzv
fi

# 8. Core diversity metrics
if [[ ! -d core-metrics ]]; then
  echo "==> [8/12] core-metrics-phylogenetic depth=$SAMPLING_DEPTH"
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza --i-table table-bio.qza \
    --p-sampling-depth "$SAMPLING_DEPTH" \
    --m-metadata-file "$META" \
    --output-dir core-metrics \
    --p-n-jobs-or-threads "$THREADS"
fi

# 9. Alpha-diversity group significance (Kruskal-Wallis)
echo "==> [9/12] Alpha group significance"
for metric in shannon_vector faith_pd_vector evenness_vector observed_features_vector; do
  [[ -f core-metrics/${metric}.qza ]] && \
  qiime diversity alpha-group-significance \
    --i-alpha-diversity "core-metrics/${metric}.qza" \
    --m-metadata-file "$META" \
    --o-visualization "core-metrics/${metric}-group-significance.qzv" || true
done

# 10. Beta-diversity PERMANOVA on the focal grouping variables
echo "==> [10/12] Beta PERMANOVA: primary-comparison, environment, mother, age-class, fed"
for matrix in weighted_unifrac_distance_matrix unweighted_unifrac_distance_matrix \
              bray_curtis_distance_matrix jaccard_distance_matrix; do
  for col in primary-comparison environment mother age-class fed; do
    [[ -f core-metrics/${matrix}.qza ]] && \
    qiime diversity beta-group-significance \
      --i-distance-matrix "core-metrics/${matrix}.qza" \
      --m-metadata-file "$META" \
      --m-metadata-column "$col" \
      --p-pairwise \
      --o-visualization "core-metrics/${matrix}-${col}-permanova.qzv" 2>/dev/null || true
  done
done

# 11. Taxa bar plots
if [[ ! -f taxa-bar-plots.qzv ]]; then
  echo "==> [11/12] Taxa bar plots"
  qiime taxa barplot \
    --i-table table-bio.qza --i-taxonomy taxonomy.qza \
    --m-metadata-file "$META" \
    --o-visualization taxa-bar-plots.qzv
fi

# 12. ANCOM-BC differential abundance for primary-comparison and environment
echo "==> [12/12] ANCOM-BC differential abundance"
mkdir -p ancombc
for col in primary-comparison environment; do
  for level in 6 5; do  # genus, family
    # Collapse to the chosen taxonomic level first
    qiime taxa collapse \
      --i-table table-bio.qza --i-taxonomy taxonomy.qza \
      --p-level "$level" \
      --o-collapsed-table "ancombc/table-L${level}.qza" 2>/dev/null || true
    qiime composition ancombc \
      --i-table "ancombc/table-L${level}.qza" \
      --m-metadata-file "$META" \
      --p-formula "$col" \
      --o-differentials "ancombc/ancombc-${col}-L${level}.qza" 2>/dev/null || true
    qiime composition da-barplot \
      --i-data "ancombc/ancombc-${col}-L${level}.qza" \
      --p-significance-threshold 0.05 \
      --o-visualization "ancombc/ancombc-${col}-L${level}.qzv" 2>/dev/null || true
  done
done

echo "==> 16S pipeline complete at $(date)"
echo "    Outputs under: $OUT"
echo "    View any .qzv at https://view.qiime2.org/"
