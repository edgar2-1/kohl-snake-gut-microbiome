// Fig. 4 — Alpha diversity by rearing treatment (n = 26).
//
// Re-laid from the superseded 4-across strip (28.0 x 8.9 cm, over the journal's
// 21.5 cm limit and unreadable after column reduction) into a 2 x 2 grid at final
// print width, with the capital panel letters the journal requires and the
// Mann-Whitney P values annotated so the figure stands without the caption.
//
//   node scripts/figures/fig4.js

const path = require('path');
const F = require(path.join(__dirname, '..', 'figlib'));

const OUT = path.join(F.RES, 'figures', 'submission', 'Fig4.png');

const groups = F.rearedGroups();
const ids = Object.keys(groups);
const clean = ids.filter((i) => groups[i] === 'Clean');
const muck = ids.filter((i) => groups[i] === 'Muck');

// Exact P values from lean_stats.json (stratified rank-sum test over all 211,680 within-litter assignments).
const PANELS = [
  { key: 'observed', letter: 'A', ylab: 'Observed ASVs', p: '0.069', sig: false },
  { key: 'faith', letter: 'B', ylab: "Faith's PD", p: '0.0008', sig: true },
  { key: 'shannon', letter: 'C', ylab: 'Shannon <tspan font-style="italic">H</tspan>&#8242; (bits)', p: '0.017', sig: true },
  { key: 'evenness', letter: 'D', ylab: "Pielou's <tspan font-style=\"italic\">J</tspan>", p: '0.006', sig: true },
];

const WMM = 175, HMM = 128;
const f = F.fig({ wmm: WMM, hmm: HMM });

const ML = 118, MR = 30, MT = 50, MB = 92, CGAP = 158, RGAP = 150;
const PW = (WMM * 10 - ML - MR - CGAP) / 2;
const PH = (HMM * 10 - MT - MB - RGAP) / 2;

PANELS.forEach((panel, k) => {
  const col = k % 2, row = (k / 2) | 0;
  const box = { x: ML + col * (PW + CGAP), y: MT + row * (PH + RGAP), w: PW, h: PH };

  const vals = F.alpha(panel.key);
  const cv = clean.map((i) => vals[i]);
  const mv = muck.map((i) => vals[i]);
  const all = [...cv, ...mv];
  const lo = Math.min(...all), hi = Math.max(...all);
  const pad = (hi - lo) * 0.12;

  const ax = F.axes(f, box, {
    xlim: [0, 1],
    ylim: [lo - pad, hi + pad],
    categorical: true,
    xticks: [],
    ylabel: panel.ylab,
    grid: true,
    nyt: 5,
  });

  const centres = [0.29, 0.71].map((t) => box.x + t * box.w);

  [[cv, F.C.CLEAN, 'Clean', clean.length], [mv, F.C.MUCK, 'Muck', muck.length]]
    .forEach(([v, colr, name, nn], gi) => {
      F.boxplot(f, ax, centres[gi], v, { width: PW * 0.30, color: colr });
      F.strip(f, ax, centres[gi], v, { color: colr, r: 6.5, amp: PW * 0.075 });
      f.text(centres[gi], box.y + box.h + 40, name, { size: F.T.tick, anchor: 'middle' });
      f.text(centres[gi], box.y + box.h + 40 + F.T.tick * 1.25, `(<tspan font-style="italic">n</tspan> = ${nn})`,
        { size: F.T.tick, anchor: 'middle' });
    });

  // Exact stratified P, bold where it survives at 0.05.
  f.text(box.x + box.w - 6, box.y + F.T.annot * 1.1,
    `<tspan font-style="italic">P</tspan> = ${panel.p}`,
    { size: F.T.annot, anchor: 'end', weight: panel.sig ? 700 : 400, fill: panel.sig ? F.C.INK : F.C.MUTED });

  f.panelLetter(box.x - 96, box.y - 12, panel.letter);
});

F.render(f, OUT, { dpi: 600 })
  .then(() => console.log('Fig4 written — 2x2, panel letters A-D, P values annotated'))
  .catch((e) => { console.error(e.message); process.exit(1); });
