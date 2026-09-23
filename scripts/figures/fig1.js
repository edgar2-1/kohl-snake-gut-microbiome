// Fig. 1 — Read retention through DADA2 denoising (32 biological samples + 4 blanks).
//
// Rebuilt on figlib. The superseded version printed the five group labels on a single
// shared baseline, so the three narrow groups (Fetus n = 3, Newborn n = 3, Blank n = 4)
// collided with each other ("Fetus (n =Ne)wborn (n = B)lank (n  4)") and ran straight
// through the rotated sample-ID tick labels, obliterating the four Blank IDs. Here the
// "(n = N)" suffix is moved to its own second line under each group name, the whole
// group band is pushed clear of the tick labels by a coloured span rule, and the
// floating "merged / denoised (%)" string inside the plot is replaced by a real legend
// with a bar swatch and a point swatch.
//
//   node scripts/figures/fig1.js

const fs = require('fs');
const path = require('path');
const F = require(path.join(__dirname, '..', 'figlib'));

const OUT = path.join(F.RES, 'figures', 'submission', 'Fig1.png');
const SRC = path.join(F.RES, 'figures', 'Fig1_dada2_retention.tsv');

// ---------------------------------------------------------------- data
// "Sterile" is the retracted label for what the manuscript calls Clean.
const LABEL = { Sterile: 'Clean', Muck: 'Muck', Fetus: 'Fetus', Newborn: 'Newborn', Blank: 'Blank' };
const ORDER = ['Sterile', 'Muck', 'Fetus', 'Newborn', 'Blank'];
const COLOR = {
  Sterile: F.C.CLEAN, Muck: F.C.MUCK, Fetus: F.C.FETUS,
  Newborn: F.C.NEWBORN, Blank: F.C.BLANK,
};

const lines = fs.readFileSync(SRC, 'utf8').replace(/\s+$/, '').split(/\r?\n/);
const H = lines[0].split('\t');
const ix = (name) => H.indexOf(name);
const rows = lines.slice(1).map((l) => {
  const p = l.split('\t');
  return {
    id: p[0],
    input: +p[ix('input')],
    denoised: +p[ix('denoised')],
    merged: +p[ix('merged')],
    nonchim: +p[ix('non-chimeric')],
    pctNonchim: +p[ix('percentage of input non-chimeric')],
    group: p[ix('group')],
  };
});
rows.forEach((r) => { r.mergeEff = r.denoised > 0 ? (100 * r.merged) / r.denoised : 0; });

// Groups in fixed order; within group, descending nonchimeric retention.
const sorted = ORDER.flatMap((g) =>
  rows.filter((r) => r.group === g).sort((a, b) => b.pctNonchim - a.pctNonchim));
const counts = {};
ORDER.forEach((g) => { counts[g] = rows.filter((r) => r.group === g).length; });

// ---------------------------------------------------------------- checks printed to console
const med = (v) => {
  const s = [...v].sort((a, b) => a - b);
  const m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
};
const reared = rows.filter((r) => r.group === 'Sterile' || r.group === 'Muck');
const early = rows.filter((r) => r.group === 'Fetus' || r.group === 'Newborn');
const rng = (v) => `${Math.min(...v).toFixed(1)}–${Math.max(...v).toFixed(1)}`;
console.log(`  n libraries            : ${rows.length} (${reared.length} reared, ${early.length} early, ${counts.Blank} blanks)`);
console.log(`  median retention       : ${med(rows.map((r) => r.pctNonchim)).toFixed(1)} %   [caption: 78.4]`);
console.log(`  reared merge efficiency: ${rng(reared.map((r) => r.mergeEff))} %   [caption: 94.1-99.1]`);
console.log(`  early  merge efficiency: ${rng(early.map((r) => r.mergeEff))} %   [caption: 8.8-78.8]`);
console.log(`  early libraries < 5,000 nonchimeric reads: ${early.filter((r) => r.nonchim < 5000).map((r) => `${r.id} (${r.nonchim})`).join(', ')}`);

// ---------------------------------------------------------------- geometry
const WMM = 175, HMM = 76;
const f = F.fig({ wmm: WMM, hmm: HMM });

const ML = 128, MR = 24, MT = 88, PH = 452;
const PW = WMM * 10 - ML - MR;          // 1598
const box = { x: ML, y: MT, w: PW, h: PH };
const AXY = MT + PH;                     // baseline of the x axis

const GAP = 46;                          // gap between groups
const SW = (PW - GAP * (ORDER.length - 1)) / sorted.length;  // slot per library
const BW = SW - 5;                       // drawn bar width

const slots = [];
let cx = box.x, prev = null;
for (const r of sorted) {
  if (prev !== null && r.group !== prev) cx += GAP;
  slots.push(cx);
  cx += SW;
  prev = r.group;
}

// ---------------------------------------------------------------- axes
const ax = F.axes(f, box, {
  xlim: [0, 1],
  ylim: [0, 105],
  yticks: [0, 20, 40, 60, 80, 100],
  ylabel: 'Percentage of reads (%)',
  categorical: true,
  xticks: [],
  grid: true,
});

// ---------------------------------------------------------------- bars, points, tick labels
const TICKY = AXY + 13;                  // top of the rotated sample-ID labels
sorted.forEach((r, i) => {
  const c = COLOR[r.group];
  const x = slots[i] + 2.5;
  const yTop = ax.sy(r.pctNonchim);
  f.rect(x, yTop, BW, AXY - yTop, { fill: c, opacity: 0.85, stroke: c, sw: F.S.bar });
  f.circle(slots[i] + SW / 2, ax.sy(r.mergeEff), 6.6,
    { fill: '#ffffff', stroke: F.C.INK, sw: 1.8 });
  f.text(slots[i] + SW / 2 + F.T.tick * 0.36, TICKY, F.esc(r.id),
    { size: F.T.tick, fill: F.C.MUTED, anchor: 'start', rotate: 90 });
});

// ---------------------------------------------------------------- group band
// Two baselines: name, then "(n = N)" below it. Nothing shares a line, and the whole
// band sits below the deepest tick label ("Blank4", 6 characters at 1.6 mm).
const RULEY = AXY + 96;
const NAMEY = RULEY + 30;
const NY = NAMEY + 26;

let start = 0;
sorted.forEach((r, i) => {
  const last = i === sorted.length - 1;
  if (!last && sorted[i + 1].group === r.group) return;
  const g = r.group;
  const x1 = slots[start] + 2.5;
  const x2 = slots[i] + 2.5 + BW;
  const mid = (x1 + x2) / 2;
  f.line(x1, RULEY, x2, RULEY, { stroke: COLOR[g], sw: 4.2 });
  f.text(mid, NAMEY, F.esc(LABEL[g]), { size: 18, anchor: 'middle', weight: 700, fill: COLOR[g] });
  f.text(mid, NY, `(<tspan font-style="italic">n</tspan> = ${counts[g]})`,
    { size: F.T.tick, anchor: 'middle', fill: F.C.MUTED });
  start = i + 1;
});

f.text(box.x + box.w / 2, NY + 34, 'Sequencing library',
  { size: F.T.axis, anchor: 'middle' });

// ---------------------------------------------------------------- legend
// Replaces the floating "merged / denoised (%)" string that used to sit in the plot.
const LY = 34;                            // text baseline
const swx = box.x, sww = 34, swh = 22;
ORDER.forEach((g, k) => {
  f.rect(swx + (k * sww) / ORDER.length, LY - swh + 5, sww / ORDER.length + 0.4, swh,
    { fill: COLOR[g], opacity: 0.85 });
});
f.rect(swx, LY - swh + 5, sww, swh, { fill: 'none', stroke: F.C.INK, sw: 1.2 });
f.text(swx + sww + 14, LY, 'Bars, nonchimeric reads (% of input), colored by group',
  { size: F.T.legend, fill: F.C.INK });

const cx2 = swx + 600;
f.circle(cx2, LY - 5.5, 6.6, { fill: '#ffffff', stroke: F.C.INK, sw: 1.8 });
f.text(cx2 + 16, LY, 'Points, merged reads (% of denoised)',
  { size: F.T.legend, fill: F.C.INK });

// ---------------------------------------------------------------- render
F.render(f, OUT, { dpi: 600 })
  .then(() => console.log('Fig1 written - group labels on their own two-line band, real legend'))
  .catch((e) => { console.error(e.message); process.exit(1); });
