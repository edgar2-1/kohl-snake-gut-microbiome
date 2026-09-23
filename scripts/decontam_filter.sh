#!/usr/bin/env bash
# decontam_filter.sh -- prevalence-based contaminant filtering using KC/KC2 negatives.
#
# Usage:
#   bash decontam_filter.sh <table.qza> <rep-seqs.qza> <metadata.tsv> <out-dir>
#
# Outputs into <out-dir>:
#   table-decontam.qza
#   rep-seqs-decontam.qza
#   decontam-stats.tsv     (per-ASV: in_neg, in_sample, prevalence-ratio, kept)
#
# Method: an ASV is flagged as a contaminant if its prevalence (proportion of
# samples it appears in) is HIGHER in negative controls than in true samples.
# This is the standard prevalence-based decontam approach (Davis et al. 2018,
# Microbiome) implemented natively against the biom table that QIIME 2 exports.

set -euo pipefail
TABLE="$1"
REPSEQS="$2"
META="$3"
OUTDIR="$4"

cd "$OUTDIR"
rm -rf _decontam_work
mkdir -p _decontam_work
qiime tools export --input-path "$TABLE"   --output-path _decontam_work/table-export
qiime tools export --input-path "$REPSEQS" --output-path _decontam_work/seqs-export
biom convert -i _decontam_work/table-export/feature-table.biom \
             -o _decontam_work/table.tsv --to-tsv

export META_PATH="$META"
python3 - <<'PYEOF'
import pandas as pd, sys, os
work = "_decontam_work"
tbl = pd.read_csv(f"{work}/table.tsv", sep="\t", skiprows=1, index_col=0)
# Columns are sample IDs, rows are feature IDs.
meta = pd.read_csv(os.environ.get("META_PATH"), sep="\t")
meta = meta[meta["sample-id"] != "#q2:types"]
neg_ids = meta.loc[meta["is-negative-control"].astype(str).str.lower() == "true", "sample-id"].tolist()
pos_ids = meta.loc[meta["is-negative-control"].astype(str).str.lower() == "false", "sample-id"].tolist()
neg_in = [s for s in neg_ids if s in tbl.columns]
pos_in = [s for s in pos_ids if s in tbl.columns]

# Prevalence (proportion of samples with non-zero count)
prev_neg = (tbl[neg_in] > 0).sum(axis=1) / max(len(neg_in), 1)
prev_pos = (tbl[pos_in] > 0).sum(axis=1) / max(len(pos_in), 1)

# Davis et al. prevalence score: higher means more likely real
# Use a simple decision: contaminant if prev_neg > prev_pos AND it appears in >=1 neg
is_contam = (prev_neg > prev_pos) & (prev_neg > 0)
stats = pd.DataFrame({
    "in_neg_prevalence": prev_neg.round(3),
    "in_sample_prevalence": prev_pos.round(3),
    "flag_contaminant": is_contam,
})
stats.index.name = "feature-id"
stats.to_csv("decontam-stats.tsv", sep="\t")

kept = stats.index[~stats["flag_contaminant"]].tolist()
with open(f"{work}/keep-ids.tsv", "w") as fh:
    fh.write("feature-id\n")
    for i in kept:
        fh.write(i + "\n")

n_total = len(stats); n_drop = int(stats["flag_contaminant"].sum())
print(f"decontam: {n_drop}/{n_total} features flagged as contaminants ({100*n_drop/max(n_total,1):.1f}%)")
PYEOF

# Filter table + seqs to keep only non-contaminant features, AND drop the
# negative-control samples themselves (we don't carry them into diversity).
qiime feature-table filter-features \
  --i-table "$TABLE" \
  --m-metadata-file _decontam_work/keep-ids.tsv \
  --o-filtered-table _decontam_work/table-filt.qza

qiime feature-table filter-samples \
  --i-table _decontam_work/table-filt.qza \
  --m-metadata-file "$META" \
  --p-where "[is-negative-control]='false'" \
  --o-filtered-table table-decontam.qza

qiime feature-table filter-seqs \
  --i-data "$REPSEQS" \
  --m-metadata-file _decontam_work/keep-ids.tsv \
  --o-filtered-data rep-seqs-decontam.qza

echo "decontam_filter: wrote table-decontam.qza + rep-seqs-decontam.qza + decontam-stats.tsv"
