// exact_stratified_tests.js -- the one Clean-versus-Muck analysis of the manuscript.
//
// Every treatment contrast is evaluated by comparing siblings only with siblings: the test
// statistic is computed under all 211,680 ways the treatments could have been assigned within
// litters (252 x 56 x 15 x 1 for dams 11.21, 11.22, NR1908, NR1912), so P values are exact.
//   * composition:      PERMANOVA pseudo-F on each distance matrix (R2 as effect size)
//   * alpha diversity:  sum over litters of the within-litter Mann-Whitney U (a stratified rank-sum test after van Elteren, 1960, with equal litter weights)
//   * genus screen:     the same statistic on relative abundance, BH across the 328 features
//   * within one litter: the same enumeration restricted to a single litter (Table S2b)
// Maternal identity: PERMANOVA with dam as a four-level factor, 9,999 unbiased seeded permutations.
//
// Inputs (relative to results/):  figures/submission/TableS4.tsv  (dam, treatment, feeding)
//   scripts/qzv/alpha_{shannon,observed_features,faith_pd,evenness}.tsv  (rarefied, 5,000 reads)
//   scripts/qzv/dist_{jaccard,bray_curtis,unweighted_unifrac,weighted_unifrac}.tsv  (pairwise)
//   scripts/qzv/level-6.csv  (unrarefied genus-level counts)
// Outputs: figures/submission/TableS1.tsv, TableS2.tsv, figures/Fig5_DA_genera_table.tsv,
//   and lean_stats.json (every number quoted in the text) in results/.
//   Run:  node scripts/exact_stratified_tests.js
"use strict";
const fs = require("fs"), path = require("path");
const RES = path.resolve(__dirname, "..");
const P = (...p) => path.join(RES, ...p);
const norm = (s) => { s = String(s).trim(); const m = s.match(/^SS0*(\d+)$/); return m ? "SS" + m[1].padStart(3, "0") : s; };

// ---------------------------------------------------------------- design (Table S4)
const s4 = fs.readFileSync(P("figures/submission/TableS4.tsv"), "utf8").split("\n").filter((l) => l && !l.startsWith("#")).map((l) => l.split("\t"));
const H = s4[0], col = (n) => H.indexOf(n);
const reared = s4.slice(1).filter((r) => ["Clean", "Muck"].includes(r[col("treatment")]))
  .map((r) => ({ id: r[col("library")], dam: r[col("dam")], muck: r[col("treatment")] === "Muck", fed: r[col("fed")] === "Yes", days: +r[col("exposure_days")] }));
const dams = [...new Set(reared.map((r) => r.dam))];
const combos = (n, k) => { const out = []; const rec = (s, acc) => { if (acc.length === k) { out.push(acc.slice()); return; } for (let i = s; i < n; i++) { acc.push(i); rec(i + 1, acc); acc.pop(); } }; rec(0, []); return out; };
const litters = dams.map((d) => { const m = reared.filter((r) => r.dam === d); const k = m.filter((r) => r.muck).length;
  return { dam: d, members: m, k, splits: combos(m.length, k), obs: m.map((r, i) => (r.muck ? i : -1)).filter((i) => i >= 0) }; });
const TOTAL = litters.reduce((p, l) => p * l.splits.length, 1);
console.log("design:", litters.map((l) => `${l.dam} ${l.members.length} (${l.k} Muck, ${l.splits.length} splits)`).join("; "), "| assignments:", TOTAL);

// ---------------------------------------------------------------- exact stratified rank-sum test
const uOf = (vals, muckIdx) => { const set = new Set(muckIdx); let u = 0; for (const i of muckIdx) for (let j = 0; j < vals.length; j++) { if (set.has(j)) continue; u += vals[i] > vals[j] ? 1 : vals[i] === vals[j] ? 0.5 : 0; } return u; };
function stratExact(valueOf) {
  let hist = new Map([[0, 1]]), T = 0, E = 0;
  for (const l of litters) {
    const vals = l.members.map((r) => valueOf(r.id)); T += uOf(vals, l.obs); E += l.k * (l.members.length - l.k) / 2;
    const h = new Map(); for (const sp of l.splits) { const u = uOf(vals, sp); h.set(u, (h.get(u) || 0) + 1); }
    const nh = new Map(); for (const [a, ca] of hist) for (const [b, cb] of h) nh.set(a + b, (nh.get(a + b) || 0) + ca * cb); hist = nh;
  }
  const dev = Math.abs(T - E); let c = 0; for (const [t, n] of hist) if (Math.abs(t - E) >= dev - 1e-9) c += n;
  return { T, E, count: c, P: c / TOTAL };
}
const BH = (ps) => { const n = ps.length, idx = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]); const q = new Array(n); let prev = 1; for (let r = n - 1; r >= 0; r--) { const v = Math.min(prev, idx[r][0] * n / (r + 1)); prev = v; q[idx[r][1]] = v; } return q; };
const median = (a) => { const v = a.slice().sort((x, y) => x - y); return (v[(v.length - 1) >> 1] + v[v.length >> 1]) / 2; };

// ---------------------------------------------------------------- alpha diversity
const readAlpha = (f) => { const m = {}; fs.readFileSync(P("scripts/qzv", f), "utf8").split("\n").forEach((l) => { const c = l.split("\t").map((x) => x.trim()); const id = c[0], v = c[c.length - 1]; if (id.startsWith("SS") && !isNaN(parseFloat(v))) m[norm(id)] = parseFloat(v); }); return m; };
const ALPHA = [["Faith's PD", "alpha_faith_pd.tsv"], ["Pielou's evenness", "alpha_evenness.tsv"], ["Observed ASVs", "alpha_observed_features.tsv"], ["Shannon", "alpha_shannon.tsv"]];
const alpha = ALPHA.map(([name, f]) => { const m = readAlpha(f); const r = stratExact((id) => m[id]);
  const med = (g, d) => median(reared.filter((x) => x.muck === g && (!d || x.dam === d)).map((x) => m[x.id]));
  return { name, medianClean: med(false), medianMuck: med(true), T: r.T, E: r.E, count: r.count, P: r.P,
    perLitter: litters.filter((l) => l.k > 0).map((l) => ({ dam: l.dam, clean: med(false, l.dam), muck: med(true, l.dam) })) }; });
BH(alpha.map((a) => a.P)).forEach((q, i) => (alpha[i].q = q));

// ---------------------------------------------------------------- genus screen (unrarefied level-6)
const L6 = fs.readFileSync(P("scripts/qzv/level-6.csv"), "utf8").trim().split("\n").map((l) => l.split(","));
const taxCols = L6[0].map((h, i) => [h, i]).filter(([h, i]) => i > 0 && /^d__|^Unassigned|;/.test(h));
const rows = Object.fromEntries(L6.slice(1).map((r) => [norm(r[0]), r]));
const rel = {}; let minNonzero = Infinity;
for (const r of reared) { const row = rows[r.id]; if (!row) throw new Error("no level-6 row for " + r.id); const tot = taxCols.reduce((s, [, i]) => s + (+row[i] || 0), 0);
  rel[r.id] = Object.fromEntries(taxCols.map(([h, i]) => { const v = (+row[i] || 0) / tot; if (v > 0 && v < minNonzero) minNonzero = v; return [h, v]; })); }
const c = minNonzero / 2; // pseudocount: half the smallest nonzero per-sample relative abundance
const feats = taxCols.map(([h]) => h).filter((h) => reared.some((r) => rel[r.id][h] > 0));
const RANK = ["domain", "phylum", "class", "order", "family", "genus"];
const shortName = (h) => { // Table 3 convention: genus name alone; otherwise the lowest named rank plus a qualifier
  const p = h.split(";").map((x) => x.replace(/^[a-z]__/, ""));
  const named = (x) => x && !/^_+$/.test(x) && !/^Incertae_Sedis$/i.test(x) && !/^uncultured/i.test(x);
  if (named(p[5])) return p[5].replace(/_/g, " ");
  for (let i = 4; i >= 0; i--) if (named(p[i])) return p[i].replace(/_/g, " ") + (/^Incertae_Sedis$/i.test(p[5] || "") ? " (Incertae Sedis)" : " (" + RANK[i] + ")");
  return h; };
const genus = feats.map((h) => { const r = stratExact((id) => rel[id][h]);
  const mc = reared.filter((x) => !x.muck).reduce((s, x) => s + rel[x.id][h], 0) / 16, mm = reared.filter((x) => x.muck).reduce((s, x) => s + rel[x.id][h], 0) / 10;
  return { taxon: h, short: shortName(h), meanClean: mc, meanMuck: mm, detClean: reared.filter((x) => !x.muck && rel[x.id][h] > 0).length, detMuck: reared.filter((x) => x.muck && rel[x.id][h] > 0).length,
    log2FC: Math.log2((mm + c) / (mc + c)), T: r.T, E: r.E, count: r.count, P: r.P }; });
BH(genus.map((g) => g.P)).forEach((q, i) => (genus[i].q = q));
genus.sort((a, b) => a.P - b.P || a.short.localeCompare(b.short));
const sig05 = genus.filter((g) => g.q <= 0.05), sig10 = genus.filter((g) => g.q <= 0.10);
const byShort = (n) => genus.find((g) => g.short === n) || genus.find((g) => g.taxon.endsWith(";g__" + n));

// ---------------------------------------------------------------- distances, PERMANOVA
const ids = reared.map((r) => r.id), N = ids.length;
const readDist = (f) => { const D = {}; fs.readFileSync(P("scripts/qzv", f), "utf8").split("\n").slice(1).forEach((l) => { const x = l.split("\t").map((v) => v.trim()); if (x.length < 6) return; const a = norm(x[1]), b = norm(x[2]), v = parseFloat(x[5]); if (isNaN(v)) return; (D[a] = D[a] || {})[b] = v; (D[b] = D[b] || {})[a] = v; }); for (const k in D) D[k][k] = 0;
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) if (D[ids[i]] == null || D[ids[i]][ids[j]] == null) throw new Error(f + ": missing pair " + ids[i] + "-" + ids[j]); return D; };
function pseudoF(D, groupOf, a, sub) { // groupOf(id) -> group label; sub = subset of ids
  const S = sub || ids, n = S.length; let ssT = 0; const within = {}, size = {};
  for (const id of S) size[groupOf(id)] = (size[groupOf(id)] || 0) + 1;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const d2 = D[S[i]][S[j]] ** 2; ssT += d2; const g = groupOf(S[i]); if (g === groupOf(S[j])) within[g] = (within[g] || 0) + d2; }
  ssT /= n; let ssW = 0; for (const g in size) ssW += (within[g] || 0) / size[g]; const ssA = ssT - ssW;
  return { F: (ssA / (a - 1)) / (ssW / (n - a)), R2: ssA / ssT };
}
const METRICS = [["Jaccard", "dist_jaccard.tsv"], ["Bray-Curtis", "dist_bray_curtis.tsv"], ["Unweighted UniFrac", "dist_unweighted_unifrac.tsv"], ["Weighted UniFrac", "dist_weighted_unifrac.tsv"]];
// unbiased seeded RNG for the dam test (mulberry32; floats from the full 32-bit word)
let seed = 20260914 >>> 0; const rnd = () => { seed = (seed + 0x6D2B79F5) >>> 0; let t = seed; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
const shuffle = (a) => { const b = a.slice(); for (let i = b.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [b[i], b[j]] = [b[j], b[i]]; } return b; };
const NPERM = 9999;
const composition = METRICS.map(([name, f]) => {
  const D = readDist(f); const muckSet = new Set(reared.filter((r) => r.muck).map((r) => r.id));
  const obs = pseudoF(D, (id) => (muckSet.has(id) ? "M" : "C"), 2);
  let ge = 0; const rec = (li, cur) => { if (li === litters.length) { if (pseudoF(D, (id) => (cur.has(id) ? "M" : "C"), 2).F >= obs.F - 1e-12) ge++; return; } const l = litters[li]; for (const sp of l.splits) { const s = new Set(cur); sp.forEach((i) => s.add(l.members[i].id)); rec(li + 1, s); } }; rec(0, new Set());
  // the same test inside each litter that received both treatments
  const within = litters.filter((l) => l.k > 0).map((l) => { const S = l.members.map((r) => r.id); const oSet = new Set(l.obs.map((i) => l.members[i].id));
    const o = pseudoF(D, (id) => (oSet.has(id) ? "M" : "C"), 2, S); let g = 0; for (const sp of l.splits) { const s = new Set(sp.map((i) => l.members[i].id)); if (pseudoF(D, (id) => (s.has(id) ? "M" : "C"), 2, S).F >= o.F - 1e-12) g++; } return { dam: l.dam, n: `${l.members.length - l.k} Clean / ${l.k} Muck`, R2: o.R2, F: o.F, count: g, total: l.splits.length, P: g / l.splits.length }; });
  // dam as a four-level factor
  const damOf = Object.fromEntries(reared.map((r) => [r.id, r.dam])); const dObs = pseudoF(D, (id) => damOf[id], dams.length); let dge = 0;
  for (let p = 0; p < NPERM; p++) { const perm = shuffle(ids); const m = Object.fromEntries(ids.map((id, i) => [id, damOf[perm[i]]])); if (pseudoF(D, (id) => m[id], dams.length).F >= dObs.F - 1e-12) dge++; }
  return { metric: name, environment: { F: obs.F, R2: obs.R2, count: ge, total: TOTAL, P: ge / TOTAL }, within, dam: { F: dObs.F, R2: dObs.R2, P: (dge + 1) / (NPERM + 1), permutations: NPERM } };
});
// sibling versus non-sibling mean distances (descriptive)
const sibDist = METRICS.map(([name, f]) => { const D = readDist(f); let w = [], b = []; for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) (reared[i].dam === reared[j].dam ? w : b).push(D[ids[i]][ids[j]]); const mn = (a) => a.reduce((s, v) => s + v, 0) / a.length; return { metric: name, within: mn(w), between: mn(b) }; });

// ---------------------------------------------------------------- treatment balance: body size and age, same test
const trait = (name) => { const c0 = col(name); const m = {}; s4.slice(1).forEach((r) => { if (["Clean", "Muck"].includes(r[col("treatment")])) m[r[col("library")]] = +r[c0]; }); return m; };
const svl = trait("SVL_cm"), tl = trait("total_length_cm"), mass = trait("mass_g"), age = trait("exposure_days");
const lx = ids.map((id) => Math.log(svl[id])), ly = ids.map((id) => Math.log(mass[id])); const mx = lx.reduce((a, b) => a + b) / N, my = ly.reduce((a, b) => a + b) / N;
const slope = lx.reduce((s, x, k) => s + (x - mx) * (ly[k] - my), 0) / lx.reduce((s, x) => s + (x - mx) ** 2, 0); const cond = Object.fromEntries(ids.map((id, k) => [id, ly[k] - (my + slope * (lx[k] - mx))]));
const balance = [["SVL", svl], ["total length", tl], ["body mass", mass], ["body condition", cond], ["age at sampling", age]].map(([name, m]) => { const r = stratExact((id) => m[id]); return { name, T: r.T, P: r.P, medianClean: median(reared.filter((x) => !x.muck).map((x) => m[x.id])), medianMuck: median(reared.filter((x) => x.muck).map((x) => m[x.id])) }; });

// ---------------------------------------------------------------- descriptive numbers used in the text
const A = byShort("Aeromonas"), M = byShort("Morganella");
const perLitterAero = litters.filter((l) => l.k > 0).map((l) => ({ dam: l.dam, muck: `${l.members.filter((r) => r.muck && rel[r.id][A.taxon] > 0).length}/${l.k}`, clean: `${l.members.filter((r) => !r.muck && rel[r.id][A.taxon] > 0).length}/${l.members.length - l.k}` }));
const aeroByDays = [[7, 9], [14, 14], [19, 19]].map(([lo, hi]) => { const s = reared.filter((r) => r.muck && r.days >= lo && r.days <= hi); return { days: lo === hi ? String(lo) : `${lo}-${hi}`, n: s.length, meanPct: s.reduce((a, r) => a + rel[r.id][A.taxon], 0) / s.length * 100 }; });
const morgPerLitter = litters.filter((l) => l.k > 0).map((l) => ({ dam: l.dam, cleanPct: l.members.filter((r) => !r.muck).reduce((a, r) => a + rel[r.id][M.taxon], 0) / (l.members.length - l.k) * 100, muckPct: l.members.filter((r) => r.muck).reduce((a, r) => a + rel[r.id][M.taxon], 0) / l.k * 100 }));
const envGain = sig05.filter((g) => g.meanMuck > g.meanClean && !/Bacteroidales|Parabacteroides/.test(g.short));
const summedEnvGainPct = envGain.reduce((s, g) => s + g.meanMuck, 0) * 100;

// ---------------------------------------------------------------- write tables
const pct = (v) => (v === 0 ? "0" : v * 100 < 0.01 ? (v * 100).toFixed(4) : (v * 100).toFixed(3));
const fmtP = (p) => (p < 0.0001 ? p.toExponential(1).replace("e-", " × 10^-") : p < 0.001 ? p.toFixed(4) : p.toFixed(3));
// Table S1: all features
{ const out = ["# Table S1. Genus-level differential abundance between clean- and muck-reared siblings: all " + genus.length + " genus-level features with reads in the 26 reared animals.",
  "# Statistic T = sum over litters of the within-litter Mann-Whitney U (Muck > Clean pairs; ties 0.5), evaluated under all " + TOTAL + " within-litter treatment assignments; P is exact, two-sided; q = Benjamini-Hochberg across the " + genus.length + " features.",
  "# Means are percentages of reads (unrarefied filtered table); carriers = animals with at least one read (of 16 Clean, 10 Muck); log2FC = log2[(mean Muck + c)/(mean Clean + c)] with c = " + c.toExponential(3) + ".",
  ["feature", "taxonomy", "mean_clean_pct", "mean_muck_pct", "carriers_clean", "carriers_muck", "log2FC", "T", "T_null_mean", "assignments_as_extreme", "P_exact", "q"].join("\t")];
  genus.forEach((g) => out.push([g.short, g.taxon, (g.meanClean * 100).toFixed(5), (g.meanMuck * 100).toFixed(5), g.detClean, g.detMuck, g.log2FC.toFixed(2), g.T, g.E, g.count, g.P.toExponential(3), g.q.toFixed(4)].join("\t")));
  fs.writeFileSync(P("figures/submission/TableS1.tsv"), out.join("\n") + "\n"); }
// Fig. 5/6 source table (same columns the figure scripts read)
{ const out = ["taxon\tmean_sterile\tmean_muck\tU\tp_value\tq_value\tlog2FC_muck_over_sterile\tshort"];
  genus.forEach((g) => out.push([g.taxon, g.meanClean, g.meanMuck, g.T, g.P, g.q, g.log2FC, g.short].join("\t")));
  fs.writeFileSync(P("figures/Fig5_DA_genera_table.tsv"), out.join("\n") + "\n"); }
// Table S2
{ const out = ["# Table S2. The composition, dispersion and maternal-identity analyses on all four distance metrics, and the treatment contrast within each litter.",
  "# (a) Environment: PERMANOVA pseudo-F(1,24) and R2 for Clean vs Muck (n = 26), P exact from all " + TOTAL + " within-litter treatment assignments. The last two columns are mean pairwise distances between siblings and between non-siblings. Dispersion: PERMDISP on distances to the group spatial median, 9,999 permutations (scripts/permdisp_validated.js). Dam: PERMANOVA with dam as a four-level factor, F(3,22), 9,999 permutations (smallest attainable P = 0.0001).",
  "metric\tenv_F\tenv_R2\tenv_assignments_as_extreme\tenv_P_exact\tpermdisp_F\tpermdisp_P\tdam_F\tdam_R2\tdam_P\tmean_dist_within_litter\tmean_dist_between_litters"];
  const permdisp = { "Jaccard": [1.30, 0.26], "Bray-Curtis": [1.30, 0.17], "Unweighted UniFrac": [10.89, 0.005], "Weighted UniFrac": [1.51, 0.17] };
  composition.forEach((m) => out.push([m.metric, m.environment.F.toFixed(3), m.environment.R2.toFixed(3), m.environment.count, m.environment.P.toExponential(2), permdisp[m.metric][0], permdisp[m.metric][1], m.dam.F.toFixed(3), m.dam.R2.toFixed(3), m.dam.P.toFixed(4), sibDist.find((s) => s.metric === m.metric).within.toFixed(3), sibDist.find((s) => s.metric === m.metric).between.toFixed(3)].join("\t")));
  out.push("#", "# (b) The same PERMANOVA inside each litter that received both treatments, P exact from every assignment within that litter (252, 56 and 15 assignments; smallest attainable P = 0.0079, 0.018 and 0.067). Values are R2 (P).", "dam\tallocation\t" + METRICS.map((m) => m[0]).join("\t"));
  composition[0].within.forEach((w, i) => out.push([w.dam, w.n, ...composition.map((m) => `${m.within[i].R2.toFixed(2)} (${m.within[i].P.toFixed(4)})`)].join("\t")));
  out.push("#", "# (c) Alpha diversity: medians by treatment within each litter that received both treatments (rarefied to 5,000 reads).", "dam\t" + alpha.map((a) => a.name + "_Clean\t" + a.name + "_Muck").join("\t"));
  alpha[0].perLitter.forEach((_, i) => out.push([alpha[0].perLitter[i].dam, ...alpha.flatMap((a) => [a.perLitter[i].clean.toFixed(3), a.perLitter[i].muck.toFixed(3)])].join("\t")));
  fs.writeFileSync(P("figures/submission/TableS2.tsv"), out.join("\n") + "\n"); }

// ---------------------------------------------------------------- every number quoted in the text
const table3 = sig05.map((g) => ({ feature: g.short, clean_pct: pct(g.meanClean), muck_pct: pct(g.meanMuck), carriers: `${g.detClean}, ${g.detMuck}`, log2FC: (g.log2FC > 0 ? "+" : "−") + Math.abs(g.log2FC).toFixed(2), q: g.q.toFixed(3), P: g.P }));
const stats = { assignments: TOTAL, pseudocount: c, litters: litters.map((l) => ({ dam: l.dam, n: l.members.length, muck: l.k, splits: l.splits.length })),
  composition, sibDist, balance, alpha, genus_counts: { tested: genus.length, q05: sig05.length, q10: sig10.length, named_genera_q05: sig05.filter((g) => /;g__[A-Za-z]/.test(g.taxon) && !/Incertae/.test(g.short)).length },
  table3, aeromonas: { P: A.P, q: A.q, meanClean: A.meanClean, meanMuck: A.meanMuck, log2FC: A.log2FC, carriers: [A.detClean, A.detMuck], perLitter: perLitterAero, muckRangePct: reared.filter((r) => r.muck && rel[r.id][A.taxon] > 0).map((r) => rel[r.id][A.taxon] * 100).sort((a, b) => a - b), fedClean: reared.filter((r) => !r.muck && r.fed).map((r) => rel[r.id][A.taxon] * 100), byDays: aeroByDays },
  morganella: { P: M.P, q: M.q, meanClean: M.meanClean, meanMuck: M.meanMuck, log2FC: M.log2FC, perLitter: morgPerLitter }, summedEnvGainPct, envGain: envGain.map((g) => g.short),
  named: Object.fromEntries(["Salmonella", "Edwardsiella", "Cetobacterium", "Bacteroides", "Parabacteroides", "Chryseobacterium", "Pedobacter", "Fuscovulum"].map((n) => { const g = byShort(n); return [n, g && { P: g.P, q: g.q, meanClean: g.meanClean, meanMuck: g.meanMuck, det: [g.detClean, g.detMuck] }]; })) };
fs.writeFileSync(P("lean_stats.json"), JSON.stringify(stats, null, 1));

// ---------------------------------------------------------------- console summary
console.log("\nComposition (exact stratified P; dam PERMANOVA):"); composition.forEach((m) => console.log(`  ${m.metric.padEnd(19)} F=${m.environment.F.toFixed(3)} R2=${m.environment.R2.toFixed(3)} P=${m.environment.P.toExponential(2)} (${m.environment.count}/${TOTAL}) | dam F=${m.dam.F.toFixed(2)} R2=${m.dam.R2.toFixed(3)} P=${m.dam.P.toFixed(4)} | within: ${m.within.map((w) => `${w.dam} ${w.R2.toFixed(2)} (${w.P.toFixed(4)})`).join(", ")}`));
console.log("Balance (stratified exact):", balance.map((b) => `${b.name} ${b.medianClean} vs ${b.medianMuck} P=${b.P.toFixed(3)}`).join("; "));
console.log("Sibling vs non-sibling mean distance:", sibDist.map((s) => `${s.metric} ${s.within.toFixed(3)} vs ${s.between.toFixed(3)}`).join("; "));
console.log("\nAlpha (exact stratified):"); alpha.forEach((a) => console.log(`  ${a.name.padEnd(18)} median C ${a.medianClean.toFixed(3)} M ${a.medianMuck.toFixed(3)} | T=${a.T}/48 P=${a.P.toFixed(4)} q=${a.q.toFixed(4)} | per litter (C vs M): ${a.perLitter.map((p) => `${p.dam} ${p.clean.toFixed(2)} vs ${p.muck.toFixed(2)}`).join("; ")}`));
console.log(`\nGenus screen: ${genus.length} tested, ${sig05.length} at q<=0.05 (${stats.genus_counts.named_genera_q05} named genera), ${sig10.length} at q<=0.10; pseudocount c=${c.toExponential(3)}`);
table3.forEach((t) => console.log(`  ${t.feature.padEnd(30)} C ${t.clean_pct.padStart(6)}% M ${t.muck_pct.padStart(7)}% carriers ${t.carriers.padEnd(6)} log2FC ${t.log2FC.padStart(6)} q ${t.q}  P=${t.P.toExponential(2)}`));
console.log(`Aeromonas per litter (muck | clean): ${perLitterAero.map((p) => `${p.dam} ${p.muck} | ${p.clean}`).join("; ")}; fed Clean = ${stats.aeromonas.fedClean.join(",")}; by days: ${aeroByDays.map((d) => `${d.days} d n=${d.n} mean ${d.meanPct.toFixed(2)}%`).join("; ")}`);
console.log(`Morganella: P=${M.P.toExponential(2)} q=${M.q.toFixed(3)} log2FC=${M.log2FC.toFixed(2)}; per litter: ${morgPerLitter.map((m) => `${m.dam} ${m.cleanPct.toFixed(1)}% -> ${m.muckPct.toFixed(2)}%`).join("; ")}`);
console.log(`Summed Muck mean of significantly muck-enriched environmental features (${envGain.map((g) => g.short).join(", ")}): ${summedEnvGainPct.toFixed(2)}%`);
console.log("Named:", Object.entries(stats.named).map(([n, g]) => `${n} P=${g.P.toFixed(3)} q=${g.q.toFixed(3)}`).join("; "));
console.log("\nwrote figures/submission/TableS1.tsv, TableS2.tsv, figures/Fig5_DA_genera_table.tsv, lean_stats.json");
