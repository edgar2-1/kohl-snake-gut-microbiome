#!/usr/bin/env python3
"""Partition the maternal-identity (female) PERMANOVA term into facility and female within facility.

Table S2 block (d). The four litters came from two facilities (females 11.21 and 11.22 at one,
NR1908 and NR1912 at the other), so the 3-df female term is split into facility (1 df) and
female within facility (2 df) on the same 26 reared animals and distance matrices used by
exact_stratified_tests.js.

    python scripts/facility_partition.py [out.tsv]

P for female within facility: 9,999 permutations of female labels within facility.
P for facility: 9,999 free permutations (descriptive only; two females per facility).
Seed 1. Requires numpy and pandas.
"""
import os
import sys

import numpy as np
import pandas as pd

RES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RES, "figures", "submission", "TableS2_facility_partition.tsv")

meta = pd.read_csv(os.path.join(RES, "meta", "metadata.tsv"), sep="\t", comment="#")
# metadata.tsv labels the clean-reared group "Sterile" (historical label; "Clean" in the manuscript)
reared = meta[meta["primary-comparison"].isin(["Sterile", "Clean", "Muck"])].copy()
assert len(reared) == 26, len(reared)
reared["facility"] = reared["mother"].map(lambda d: "Jayne" if str(d).startswith("11.") else "Secor")
ids = list(reared["sample-id"])


def load_dist(name):
    df = pd.read_csv(os.path.join(RES, "scripts", "qzv", f"dist_{name}.tsv"), sep="\t", index_col=0)
    D = pd.DataFrame(0.0, index=ids, columns=ids)
    for a, b, dist in zip(df.SubjectID1, df.SubjectID2, df.Distance):
        if a in D.index and b in D.index:
            D.loc[a, b] = D.loc[b, a] = dist
    return D.values


def gower(D):
    A = -0.5 * D ** 2
    n = len(D)
    J = np.eye(n) - np.ones((n, n)) / n
    return J @ A @ J


def hat(X):
    return X @ np.linalg.pinv(X.T @ X) @ X.T


def design(labels):
    return pd.get_dummies(pd.Series(labels)).values.astype(float)


rng = np.random.default_rng(1)
rows = []
n = 26
B = 9999
fac = reared["facility"].values
fem = reared["mother"].values
for metric in ["bray_curtis", "jaccard", "unweighted_unifrac", "weighted_unifrac"]:
    G = gower(load_dist(metric))
    SST = np.trace(G)
    SSf = np.trace(hat(design(fac)) @ G)
    SSd = np.trace(hat(design(fem)) @ G)
    SSw = SSd - SSf
    Ff = (SSf / 1) / ((SST - SSd) / 22)
    Fw = (SSw / 2) / ((SST - SSd) / 22)
    cf = cw = 0
    for _ in range(B):
        perm = rng.permutation(n)
        SSf_p = np.trace(hat(design(fac[perm])) @ G)
        SSd_p = np.trace(hat(design(fem[perm])) @ G)
        if (SSf_p / 1) / ((SST - SSd_p) / 22) >= Ff - 1e-12:
            cf += 1
        perm2 = np.arange(n)
        for f in np.unique(fac):
            idx = np.where(fac == f)[0]
            perm2[idx] = rng.permutation(idx)
        SSd_p2 = np.trace(hat(design(fem[perm2])) @ G)
        SSf_p2 = np.trace(hat(design(fac[perm2])) @ G)
        if ((SSd_p2 - SSf_p2) / 2) / ((SST - SSd_p2) / 22) >= Fw - 1e-12:
            cw += 1
    rows.append({
        "Metric": {"bray_curtis": "Bray-Curtis", "jaccard": "Jaccard",
                   "unweighted_unifrac": "Unweighted UniFrac", "weighted_unifrac": "Weighted UniFrac"}[metric],
        "R2_dam": SSd / SST,
        "R2_facility": SSf / SST,
        "R2_dam_within_facility": SSw / SST,
        "frac_of_dam_R2_between_facility": SSf / SSd,
        "P_dam_within_facility": (cw + 1) / (B + 1),
        "P_facility_free_perm": (cf + 1) / (B + 1),
    })

pd.DataFrame(rows).round(4).to_csv(OUT, sep="\t", index=False)
print("wrote", os.path.relpath(OUT, RES))
