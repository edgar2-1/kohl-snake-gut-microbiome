// Fig. 3 — Principal-coordinates analysis of the 26 reared animals on four
// beta-diversity metrics (Bray-Curtis, Jaccard, unweighted UniFrac, weighted UniFrac).
//
// Regenerated from the QIIME 2 ordination.txt files. The superseded version carried
// no capital panel letters, which the journal requires on every part of a multi-part
// figure; A-D are now on the four metric panels, the two rearing treatments share a
// single legend, and the whole thing is authored at 175 mm and rasterised at 600 dpi.
//
// The ordinations contain 30 samples (16 Clean, 10 Muck, 2 Fetus, 2 Newborn); only the
// 26 reared animals are plotted, selected with F.rearedGroups() (i.e. the metadata's
// `primary-comparison` column, NOT `environment`, which lumps Fetus in with Clean).
// Axis percentages come from the 30-sample "Proportion explained" row, as published.
//
//   node scripts/figures/fig3.js

const fs = require('fs');
const path = require('path');
const F = require(path.join(__dirname, '..', 'figlib'));

const OUT = path.join(F.RES, 'figures', 'submission', 'Fig3.png');
const PCOA = path.join(F.RES, 'scripts', 'pcoa_published');

// ---------------------------------------------------------------- ordination parser

// QIIME 2 ordination text format: an "Eigvals<TAB>k" line, a "Proportion explained<TAB>k"
// line followed by one row of k proportions, and a "Site<TAB>n<TAB>k" line followed by
// n rows of "sampleID<TAB>coord1<TAB>coord2<TAB>...".
function ordination(metric) {
  const lines = fs.readFileSync(path.join(PCOA, metric, 'ordination.txt'), 'utf8')
    .replace(/\r/g, '').split('\n');
  let prop = null;
  const sites = [];
  for (let i = 0; i < lines.length; i++) {
    if (lines[i].startsWith('Proportion explained\t')) {
      prop = lines[i + 1].trim().split('\t').map(Number);
    } else if (lines[i].startsWith('Site\t')) {
      const nrow = +lines[i].split('\t')[1];
      for (let j = i + 1; j <= i + nrow; j++) {
        const p = lines[j].split('\t');
        sites.push({ id: p[0], pc1: +p[1], pc2: +p[2] });
      }
    }
  }
  if (!prop || sites.length === 0) throw new Error('could not parse ' + metric);
  return { sites, pct1: prop[0] * 100, pct2: prop[1] * 100 };
}

// ---------------------------------------------------------------- data

const groups = F.rearedGroups();               // 26 reared animals: Clean 16, Muck 10
const nClean = Object.values(groups).filter((v) => v === 'Clean').length;
const nMuck = Object.values(groups).filter((v) => v === 'Muck').length;

const PANELS = [
  { metric: 'bray_curtis', letter: 'A', title: 'Bray&#8211;Curtis' },
  { metric: 'jaccard', letter: 'B', title: 'Jaccard' },
  { metric: 'unweighted_unifrac', letter: 'C', title: 'Unweighted UniFrac' },
  { metric: 'weighted_unifrac', letter: 'D', title: 'Weighted UniFrac' },
];

// ---------------------------------------------------------------- layout

const WMM = 175, HMM = 172;                     // 17.5 x 17.2 cm, inside 21.5 x 28
const f = F.fig({ wmm: WMM, hmm: HMM });

const ML = 136, MR = 26, MT = 64, CGAP = 152, RGAP = 176;
const PW = (WMM * 10 - ML - MR - CGAP) / 2;
const PH = 660;

// ---------------------------------------------------------------- panels

PANELS.forEach((panel, k) => {
  const col = k % 2, row = (k / 2) | 0;
  const box = { x: ML + col * (PW + CGAP), y: MT + row * (PH + RGAP), w: PW, h: PH };

  const ord = ordination(panel.metric);
  const pts = ord.sites
    .filter((s) => groups[s.id])                // reared animals only
    .map((s) => ({ ...s, grp: groups[s.id] }));

  const xs = pts.map((p) => p.pc1), ys = pts.map((p) => p.pc2);
  const xlo = Math.min(...xs), xhi = Math.max(...xs);
  const ylo = Math.min(...ys), yhi = Math.max(...ys);
  const xpad = (xhi - xlo) * 0.08, ypad = (yhi - ylo) * 0.08;

  const ax = F.axes(f, box, {
    xlim: [xlo - xpad, xhi + xpad],
    ylim: [ylo - ypad, yhi + ypad],
    xlabel: `PCo1 (${ord.pct1.toFixed(1)}%)`,
    ylabel: `PCo2 (${ord.pct2.toFixed(1)}%)`,
    nxt: 5, nyt: 5,
  });

  // Muck under Clean so the smaller group is never buried; white rim separates overlaps.
  const order = ['Clean', 'Muck'];
  for (const g of order) {
    for (const p of pts.filter((q) => q.grp === g)) {
      f.circle(ax.sx(p.pc1), ax.sy(p.pc2), 8,
        { fill: g === 'Clean' ? F.C.CLEAN : F.C.MUCK, stroke: '#ffffff', sw: 1.6 });
    }
  }

  f.text(box.x + box.w / 2, box.y - 16, panel.title,
    { size: F.T.title, anchor: 'middle', weight: 700 });
  f.panelLetter(box.x - 108, box.y - 8, panel.letter);
});

// ---------------------------------------------------------------- one shared legend

const LEG = [
  { color: F.C.CLEAN, label: `Clean (<tspan font-style="italic">n</tspan> = ${nClean})`, chars: 14 },
  { color: F.C.MUCK, label: `Muck (<tspan font-style="italic">n</tspan> = ${nMuck})`, chars: 13 },
];
const legY = MT + 2 * PH + RGAP + 124;          // below the row-2 axis title
const est = (c) => c * F.T.legend * 0.53;       // Helvetica advance, good enough to centre
const entryW = LEG.map((e) => 16 + 22 + est(e.chars));
const GAP = 120;
let lx = (WMM * 10 - (entryW[0] + entryW[1] + GAP)) / 2;
LEG.forEach((e, i) => {
  f.circle(lx + 8, legY - F.T.legend * 0.34, 8, { fill: e.color, stroke: '#ffffff', sw: 1.6 });
  f.text(lx + 38, legY, e.label, { size: F.T.legend, anchor: 'start' });
  lx += entryW[i] + GAP;
});

// ---------------------------------------------------------------- render

F.render(f, OUT, { dpi: 600 })
  .then(() => console.log('Fig3 written — 2x2 PCoA, panel letters A-D, one shared legend, n = 26'))
  .catch((e) => { console.error(e.message); process.exit(1); });
