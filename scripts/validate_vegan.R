#!/usr/bin/env Rscript
# Independent check of the manuscript's composition statistics against vegan, on the
# same QIIME 2 distance matrices used by results/scripts/exact_stratified_tests.js.
#
#   PROJECT_DIR=/root/kohl_v1 Rscript validate_vegan.R
#
# What is compared (26 reared neonates, Clean = primary-comparison "Sterile", n = 16; Muck, n = 10):
#   environment : adonis2 pseudo-F and R2  (the manuscript's P is exact over the 211,680
#                 within-litter assignments, so adonis2's free-permutation P is printed only for reference)
#   dam         : adonis2 pseudo-F and R2 for mother as a four-level factor, 9,999 permutations
#   dispersion  : betadisper (spatial median) F, permutest 9,999 permutations
suppressPackageStartupMessages({ library(vegan); library(permute) })

PROJECT <- Sys.getenv("PROJECT_DIR", unset = "/root/kohl_v1")
EXP <- file.path(PROJECT, "qiime", "exports")

meta <- read.table(file.path(PROJECT, "meta", "metadata.tsv"), sep = "\t", header = TRUE,
                   comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta[meta[["sample-id"]] != "#q2:types", ]
rownames(meta) <- meta[["sample-id"]]
reared <- meta[meta[["primary-comparison"]] %in% c("Sterile", "Muck"), ]

cat(R.version.string, "\n")
cat("vegan", as.character(packageVersion("vegan")), "| permute", as.character(packageVersion("permute")), "\n")
cat("reared animals:", nrow(reared),
    "| Clean =", sum(reared[["primary-comparison"]] == "Sterile"),
    "| Muck =", sum(reared[["primary-comparison"]] == "Muck"), "\n\n")

dms <- c("Bray-Curtis" = "bray_curtis", "Jaccard" = "jaccard",
         "Unweighted UniFrac" = "unweighted_unifrac", "Weighted UniFrac" = "weighted_unifrac")

# Values in the manuscript (Table 2, Table S2), from results/lean_stats.json.
ref <- list(
  "Bray-Curtis"        = c(envF = 3.657, envR2 = 0.132, damF = 4.53, damR2 = 0.382, dispF = 1.303,  dispP = 0.17),
  "Jaccard"            = c(envF = 2.405, envR2 = 0.091, damF = 3.14, damR2 = 0.300, dispF = 1.295,  dispP = 0.26),
  "Unweighted UniFrac" = c(envF = 5.379, envR2 = 0.183, damF = 2.40, damR2 = 0.246, dispF = 10.885, dispP = 0.005),
  "Weighted UniFrac"   = c(envF = 5.389, envR2 = 0.183, damF = 3.81, damR2 = 0.342, dispF = 1.509,  dispP = 0.17))

set.seed(20260916)
ok <- TRUE
for (nm in names(dms)) {
  m <- as.matrix(read.table(file.path(EXP, paste0(dms[[nm]], "_distance_matrix"), "distance-matrix.tsv"),
                            sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, comment.char = ""))
  ids <- intersect(rownames(reared), rownames(m))
  d <- as.dist(m[ids, ids]); md <- reared[ids, ]
  env <- factor(md[["primary-comparison"]]); dam <- factor(md[["mother"]])

  a_env <- adonis2(d ~ env, permutations = 9999)
  a_dam <- adonis2(d ~ dam, permutations = 9999)
  bd <- betadisper(d, env, type = "median")
  pt <- permutest(bd, permutations = 9999)
  r <- ref[[nm]]

  chk <- function(x, y, tol) abs(x - y) <= tol
  agree <- chk(a_env$F[1], r[["envF"]], 0.0015) && chk(a_env$R2[1], r[["envR2"]], 0.0015) &&
           chk(a_dam$F[1], r[["damF"]], 0.006)  && chk(a_dam$R2[1], r[["damR2"]], 0.0015) &&
           chk(pt$tab$F[1], r[["dispF"]], 0.0015)
  ok <- ok && agree

  cat(sprintf("%-18s n = %d   %s\n", nm, length(ids), if (agree) "AGREES" else "** DIFFERS **"))
  cat(sprintf("  environment  adonis2 F = %.3f  R2 = %.3f   (manuscript F = %.3f, R2 = %.3f; free-permutation P = %.4f, for reference only)\n",
              a_env$F[1], a_env$R2[1], r[["envF"]], r[["envR2"]], a_env$`Pr(>F)`[1]))
  cat(sprintf("  dam          adonis2 F = %.2f   R2 = %.3f   P = %.4f   (manuscript F = %.2f, R2 = %.3f)\n",
              a_dam$F[1], a_dam$R2[1], a_dam$`Pr(>F)`[1], r[["damF"]], r[["damR2"]]))
  cat(sprintf("  dispersion   betadisper F = %.3f  P = %.4f   (manuscript F = %.3f, P = %.3g)\n\n",
              pt$tab$F[1], pt$tab$`Pr(>F)`[1], r[["dispF"]], r[["dispP"]]))
}
cat(if (ok) "ALL FOUR METRICS AGREE with vegan to the reported precision.\n" else "AT LEAST ONE VALUE DIFFERS - do not cite until resolved.\n")
