#!/usr/bin/env bash
# Export QIIME 2 artifacts to TSV / Newick / FASTA for downstream Python analysis.
set -euo pipefail
PROJECT_DIR="${PROJECT_DIR:-/root/kohl_v1}"
OUT="$PROJECT_DIR/qiime/16s"
EXP="$PROJECT_DIR/qiime/exports"
mkdir -p "$EXP"
cd "$OUT"

echo "==> Exporting feature tables, taxonomy, tree, sequences, diversity"

# Feature tables
for tag in table table-noHost table-clean table-bio; do
  if [[ -f "${tag}.qza" ]]; then
    qiime tools export --input-path "${tag}.qza" --output-path "$EXP/${tag}"
    if [[ -f "$EXP/${tag}/feature-table.biom" ]]; then
      biom convert -i "$EXP/${tag}/feature-table.biom" -o "$EXP/${tag}/feature-table.tsv" --to-tsv
    fi
  fi
done

# Taxonomy
if [[ -f taxonomy.qza ]]; then
  qiime tools export --input-path taxonomy.qza --output-path "$EXP/taxonomy"
fi

# Rep-seqs
for tag in rep-seqs rep-seqs-clean; do
  if [[ -f "${tag}.qza" ]]; then
    qiime tools export --input-path "${tag}.qza" --output-path "$EXP/${tag}"
  fi
done

# Phylogeny
if [[ -f rooted-tree.qza ]]; then
  qiime tools export --input-path rooted-tree.qza --output-path "$EXP/rooted-tree"
fi

# Denoising stats
if [[ -f denoising-stats.qza ]]; then
  qiime tools export --input-path denoising-stats.qza --output-path "$EXP/denoising-stats"
fi

# Diversity outputs
if [[ -d core-metrics ]]; then
  for m in shannon_vector faith_pd_vector evenness_vector observed_features_vector; do
    [[ -f "core-metrics/${m}.qza" ]] && \
      qiime tools export --input-path "core-metrics/${m}.qza" --output-path "$EXP/${m}"
  done
  for dm in weighted_unifrac_distance_matrix unweighted_unifrac_distance_matrix \
            bray_curtis_distance_matrix jaccard_distance_matrix; do
    [[ -f "core-metrics/${dm}.qza" ]] && \
      qiime tools export --input-path "core-metrics/${dm}.qza" --output-path "$EXP/${dm}"
  done
  for pc in weighted_unifrac_pcoa_results unweighted_unifrac_pcoa_results \
            bray_curtis_pcoa_results jaccard_pcoa_results; do
    [[ -f "core-metrics/${pc}.qza" ]] && \
      qiime tools export --input-path "core-metrics/${pc}.qza" --output-path "$EXP/${pc}"
  done
fi

# Collapsed tables for taxa barplots
for level in 2 3 4 5 6 7; do
  if [[ ! -f "$EXP/table-L${level}/feature-table.tsv" ]]; then
    qiime taxa collapse \
      --i-table table-bio.qza --i-taxonomy taxonomy.qza \
      --p-level $level --o-collapsed-table "$EXP/table-L${level}.qza" 2>/dev/null || continue
    qiime tools export --input-path "$EXP/table-L${level}.qza" --output-path "$EXP/table-L${level}"
    biom convert -i "$EXP/table-L${level}/feature-table.biom" \
                 -o "$EXP/table-L${level}/feature-table.tsv" --to-tsv
  fi
done

echo "==> Export complete. Files under $EXP"
ls "$EXP"
