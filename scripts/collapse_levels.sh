#!/usr/bin/env bash
set -euo pipefail
cd /root/kohl_v1/qiime/16s
EXP=/root/kohl_v1/qiime/exports
for L in 2 3 4 5 6 7; do
  QZA="$EXP/table-L${L}.qza"
  DIR="$EXP/table-L${L}"
  if [[ ! -f "$DIR/feature-table.tsv" ]]; then
    echo "==> collapsing L${L}"
    qiime taxa collapse --i-table table-bio.qza --i-taxonomy taxonomy.qza \
      --p-level "$L" --o-collapsed-table "$QZA"
    qiime tools export --input-path "$QZA" --output-path "$DIR"
    biom convert -i "$DIR/feature-table.biom" -o "$DIR/feature-table.tsv" --to-tsv
  fi
done
echo "Available L tables:"
ls "$EXP" | grep -E "^table-L[0-9]" | sort
