// Fig. 6 — The 14 genus-level features differing between
// clean- and muck-reared snakes at q <= 0.05, ranked by log2 fold-change.
//
// Fixes over the superseded version (submission_ORIGINAL_20260828/Fig6.png):
//   * that file was 17.6 cm wide only because it was 0.1 cm over the 21.5 cm limit
//     in height/width bookkeeping; this one is authored at exactly 175 mm;
//   * one label read "Hyphomicrobiales_Incertae_Sedis g. Incertae Sedis" — raw SILVA
//     underscores plus a duplicated rank;
//   * the figure and Table 3 named the same six above-genus features differently.
//     Labels are now generated to Table 3's convention exactly:
//        Steroidobacteraceae (Incertae Sedis)   Pirellulaceae (Incertae Sedis)
//        SC-I-84 (Incertae Sedis)               Anaerolineaceae (Incertae Sedis)
//        Hyphomicrobiales (Incertae Sedis)      Parachlamydiaceae (family)
//     Those six are above genus and are set ROMAN; the twelve genera stay ITALIC.
//
//   node scripts/figures/fig6.js

const path = require('path');
const F = require(path.join(__dirname, '..', 'figlib'));

const SRC = path.join(F.RES, 'figures', 'Fig5_DA_genera_table.tsv');
const OUT = path.join(F.RES, 'figures', 'submission', 'Fig6.png');

// ---------------------------------------------------------------- labels

// Build a display label from the full SILVA lineage, following Table 3:
//   named genus                      -> italic genus
//   genus == Incertae_Sedis          -> "<finest named parent> (Incertae Sedis)", roman
//   genus rank absent altogether     -> "<family> (family)", roman
// The parent of an "<Order>_Incertae_Sedis" family is the order itself, so the rank
// is never duplicated in the label.
function label(taxon) {
  const r = {};
  for (const part of taxon.split(';')) {
    const m = part.match(/^([a-z])__(.*)$/);
    if (m) r[m[1]] = m[2];
  }
  const g = r.g;

  if (g && g !== 'Incertae_Sedis') {
    return { text: g.replace(/_/g, ' '), italic: true };
  }
  if (g === 'Incertae_Sedis') {
    let parent = r.f;
    // "Hyphomicrobiales_Incertae_Sedis" is a placeholder family: name the order.
    if (!parent || parent === 'Incertae_Sedis' || /_Incertae_Sedis$/.test(parent)) {
      parent = r.o || r.c || r.p;
    }
    return { text: `${parent.replace(/_/g, ' ')} (Incertae Sedis)`, italic: false };
  }
  // No g__ field at all (e.g. "...;f__Parachlamydiaceae;__").
  const fam = r.f || r.o || r.c || r.p;
  return { text: `${fam.replace(/_/g, ' ')} (family)`, italic: false };
}

// ---------------------------------------------------------------- data

const rows = F.tsv(SRC).filter((d) => +d.q_value <= 0.05);
if (rows.length !== 14) throw new Error(`expected 14 features at q <= 0.05, got ${rows.length}`);

const feats = rows
  .map((d) => {
    const lab = { text: d.short, italic: /;g__[A-Za-z]/.test(d.taxon) && !d.short.includes("(") };
    return { ...lab, lfc: +d.log2FC_muck_over_sterile, q: +d.q_value };
  })
  .sort((a, b) => b.lfc - a.lfc);

const nItal = feats.filter((d) => d.italic).length;
if (nItal !== 12) throw new Error(`expected 12 italic genus labels, got ${nItal}`);

// ---------------------------------------------------------------- geometry

const WMM = 175, HMM = 99;                 // 17.5 x 9.9 cm, inside the 21.5 x 28 cm limit
const f = F.fig({ wmm: WMM, hmm: HMM });

const ROW = 46;                            // 4.6 mm per feature row
const ML = 340, MR = 40, MT = 66;          // ML sized to the longest label + tick gap
const box = { x: ML, y: MT, w: WMM * 10 - ML - MR, h: ROW * feats.length };

const LIM = 7.6;                           // symmetric so the zero line sits centred
const ax = F.axes(f, box, {
  xlim: [-LIM, LIM],
  ylim: [0, 1],
  xticks: [-6, -4, -2, 0, 2, 4, 6],
  drawY: false,
  // real subscript 2, kept at 1.5 mm so it stays inside the journal's type floor
  xlabel: 'log<tspan dy="6" font-size="15">2</tspan><tspan dy="-6" dx="7">fold-change (Muck / Clean)</tspan>',
});

const rowY = (i) => box.y + (i + 0.5) * ROW;

// dashed zero reference
f.line(ax.sx(0), box.y, ax.sx(0), box.y + box.h, { stroke: F.C.MUTED, sw: 1.6, dash: '11 9' });

// y axis: a spine with a tick per feature; the taxon names are the tick labels
f.line(box.x, box.y, box.x, box.y + box.h, { stroke: F.C.RULE, sw: F.S.axis });

const BH = ROW * 0.66;
feats.forEach((d, i) => {
  const yc = rowY(i);
  const colour = d.lfc > 0 ? F.C.MUCK : F.C.CLEAN;
  const x0 = ax.sx(0), x1 = ax.sx(d.lfc);
  f.rect(Math.min(x0, x1), yc - BH / 2, Math.abs(x1 - x0), BH, { fill: colour });
  f.line(box.x - 8, yc, box.x, yc, { stroke: F.C.RULE, sw: F.S.axis });
  f.text(box.x - 15, yc + F.T.tick * 0.35, F.esc(d.text),
    { size: F.T.tick, anchor: 'end', italic: d.italic });
});

// ---------------------------------------------------------------- legend

// Placed low and right of zero, where every bar in that band is negative.
const LX = box.x + box.w - 235, SW = 34, SH = 20;
[['enriched in Muck', F.C.MUCK], ['enriched in Clean', F.C.CLEAN]].forEach(([txt, col], i) => {
  const ly = box.y + box.h - 96 + i * 44;
  f.rect(LX, ly - SH / 2, SW, SH, { fill: col });
  f.text(LX + SW + 14, ly + F.T.legend * 0.36, F.esc(txt), { size: F.T.legend });
});

// No in-figure title: Journal of Herpetology figures carry none, because the legend states
// the content. The title that used to sit here duplicated the caption's first clause.
// Fig. 1 was built the same way, so the set is consistent.

// ---------------------------------------------------------------- render

F.render(f, OUT, { dpi: 600 })
  .then(() => {
    console.log(`  14 features, ${nItal} italic genera + ${14 - nItal} above-genus (roman)`);
    console.log(`  top ${feats[0].text} ${feats[0].lfc.toFixed(2)} | bottom ${feats[13].text} ${feats[13].lfc.toFixed(2)}`);
  })
  .catch((e) => { console.error(e.message); process.exit(1); });
