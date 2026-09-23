# Supplementary Methods — contaminant-screen validation and reproducibility

As a validation, the full pipeline was repeated with the extraction control included, denoising all 37 libraries together; the re-run recovered 1,050 of the 1,052 analyzed ASVs, the remaining two being the features its extended prevalence screen removed. FastTree and rarefaction are stochastic, so phylogenetic-diversity and UniFrac values are not bit-reproducible across runs, and re-run values differed slightly. The two animals of the unsplit litter (dam NR1912) enter the PERMANOVA statistic but are never relabeled under the within-litter enumeration, and contribute nothing to the rank-sum statistic. The exact 211,680 within-litter assignments factor as 252 × 56 × 15 × 1 over the four litters; the smallest attainable exact P is 0.008 within dam 11.21 (252 assignments) and 0.0001 for the dam PERMANOVA (9,999 permutations).

Of the five ASVs in the extraction control, two occurred in biological samples, at 0.001% and 0.004% of biological reads. The 37-library validation screen flagged seven features; removing all seven changed mean per-sample depth by two reads. The prevalence and frequency screens flagged no genus reported as differentially abundant.

## Partition of the dam term

Because each facility held two dams, the dam sum of squares from the PERMANOVA (dam as a four-level factor) was partitioned into facility (1 df) and dam within facility (2 df) on the Gower-centred distance matrices (McArdle and Anderson 2001). Dam within facility was tested by permuting sample labels within facility (9,999 permutations). Results are in Table S2(d).
