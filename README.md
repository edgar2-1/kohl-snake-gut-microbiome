# Neonate *Nerodia rhombifer* gut microbiome: rearing environment and maternal identity

Data and code for:

> Rosado-Ramos, E. L., B. C. Jayne, K. D. Kohl, S. M. Secor, and T. J. Colston. Rearing environment and maternal identity shape the gut microbiome of neonate Diamond-backed Watersnakes (*Nerodia rhombifer*). *Journal of Herpetology* (submitted).

Captive-born neonates were reared individually for 7–19 d in clean cages (Clean, n = 16) or on pond substrate and water (Muck, n = 10), and the V4 region of the 16S rRNA gene was sequenced from each snake's large intestine (Illumina MiSeq, 2 × 251 bp). Three full-term fetuses and three newborns describe the community at and shortly before birth.

Raw reads: NCBI SRA BioProject **[accession pending]** (37 libraries: 32 biological samples, 4 no-template PCR controls, 1 extraction control).

## Layout

| Path | Contents |
|---|---|
| `meta/metadata.tsv` | Per-sample metadata used by the scripts (36 libraries). The clean-reared group carries its historical label `Sterile`; it is "Clean" in the manuscript, and the cages were washed, not sterilized. |
| `meta/metadata_for_deposit.tsv` | The same metadata as deposited with SRA (37 libraries incl. extraction control SS120; `Sterile` relabelled `Clean`). Built by `scripts/make_deposit_metadata.js`. |
| `meta/manifest_16s.tsv`, `meta/fastq_files_to_deposit.txt` | QIIME 2 import manifest (paths on the analysis host) and the 74 FASTQ files belonging to the 37 libraries. |
| `meta/asv_all.tsv`, `meta/asv_clean.tsv`, `meta/decontam-stats.tsv` | ASV count tables before and after contaminant removal, and per-ASV prevalence-screen results. |
| `meta/pipeline.log` | Full log of the QIIME 2 run that produced `qiime_artifacts/`. |
| `qiime_artifacts/` | QIIME 2 artifacts: feature tables (`table*.qza`), representative sequences, SILVA 138.2 taxonomy, rooted and unrooted trees, denoising statistics, and the eight core-metrics visualizations plus `taxa-bar-plots.qzv` from which `scripts/qzv/` is extracted. |
| `scripts/qzv/` | Distance matrices (Bray–Curtis, Jaccard, unweighted and weighted UniFrac; long format), alpha-diversity vectors, and the genus-level (level-6) table, extracted from the `.qzv` files. |
| `scripts/pcoa_published/` | PCoA ordinations used for Fig. 3. |
| `lean_stats.json` | Every exact within-litter test statistic quoted in the text and in Tables 2–3. |
| `vegan_validation.txt` | Independent check of the PERMANOVA and PERMDISP statistics with vegan 2.7-3 (R 4.5.2). |
| `extraction_control_results.txt` | Sensitivity re-run of the pipeline with the extraction control included (Supplementary Methods). |
| `figures/submission/` | Figs. 1–6 and Tables S1–S5 as submitted; `TableS2_facility_partition.tsv` is Table S2 block (d). |
| `figures/Fig5_DA_genera_table.tsv` | Genus-level differential-abundance results (Table 3). |
| `docs/Supplementary_Methods.md` | Supplementary Methods text. |

## Scripts

**Sequence processing** (QIIME 2 2026.4 on a Linux host; outputs are included, so these need not be re-run to reproduce the statistics)

- `scripts/run_16s.sh` — import, cutadapt primer removal, DADA2 (truncation 220/180), SILVA 138.2 classification (consensus VSEARCH), removal of mitochondria/chloroplast/Eukaryota and phylum-unassigned ASVs, contaminant screen, MAFFT/FastTree phylogeny, core diversity metrics at 5,000 reads.
- `scripts/decontam_filter.sh` — prevalence-based contaminant screen against the four no-template controls (called by `run_16s.sh`).
- `scripts/export_results.sh`, `scripts/collapse_levels.sh` — export artifacts to TSV/Newick/FASTA and collapse to genus level.
- `scripts/rebuild_qzv_inputs.js`, `scripts/export_asv_table.mjs` — regenerate `scripts/qzv/` and the ASV tables from the artifacts in `qiime_artifacts/` without QIIME 2.

**Statistics** (Node.js ≥ 18; Python 3 with numpy and pandas; R with vegan)

- `scripts/exact_stratified_tests.js` — the Clean-versus-Muck analysis: stratified PERMANOVA on the four distance matrices, stratified rank-sum tests on alpha diversity and on 328 genus-level features (Benjamini–Hochberg), the within-litter contrast for female 11.21 (252 assignments), and the maternal-identity PERMANOVA. Reads Table S4 for per-sample metadata and writes `lean_stats.json`, `figures/Fig5_DA_genera_table.tsv`, and Tables S1 and S2 (blocks a–c).
- `scripts/facility_partition.py` — partitions the maternal-identity term into facility and female within facility; writes `TableS2_facility_partition.tsv`, whose rows appear as block (d) at the end of `TableS2.tsv`.
- `scripts/validate_vegan.R` — reproduces the composition statistics with vegan's `adonis2` and `betadisper`; needs the QIIME 2 distance-matrix exports (`PROJECT_DIR=... Rscript scripts/validate_vegan.R`).
- `scripts/build_storms.js` — the STORMS reporting checklist (Table S5).
- `scripts/make_deposit_metadata.js` — builds `meta/metadata_for_deposit.tsv`.

**Figures**

- `scripts/figures/fig1.js` … `fig6.js` (shared primitives in `scripts/figlib.js`) — Figs. 1–6, written to `figures/submission/`.

## Reproducing

```bash
npm install                                  # sharp (figures), h5wasm (.qza export)
node scripts/exact_stratified_tests.js       # lean_stats.json, Table 3, Tables S1 and S2 (a-c)
python scripts/facility_partition.py         # Table S2 block (d)
node scripts/build_storms.js                 # Table S5
for i in 1 2 3 4 5 6; do node scripts/figures/fig$i.js; done
```

`facility_partition.py` uses a fixed seed and reproduces `figures/submission/TableS2_facility_partition.tsv` byte for byte. Table S3 (early-sample ASV overlap at 1,000 reads) and Table S4 (per-sample metadata) are provided as computed; the script that produced Table S3 is not included.

## AI assistance

The statistical analyses and the preparation of this repository were assisted by Claude (Anthropic). The authors reviewed the code and results and take responsibility for them.

## License

Code: MIT (see `LICENSE`). Data, tables and figures: CC BY 4.0.

## Contact

Timothy J. Colston, Department of Biology, University of Puerto Rico–Mayagüez (timothy.colston@upr.edu).
