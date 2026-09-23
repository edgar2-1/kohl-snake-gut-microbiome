// figlib.js — shared SVG figure primitives for the Kohl J Herpetol submission.
//
// Replaces the matplotlib pipeline (make_figures_final.py), which cannot run on this
// machine — there is no Python with matplotlib. Everything here is Node + sharp.
//
// Geometry contract: the SVG viewBox is in TENTHS OF A MILLIMETRE, so a font-size of
// 16 is 1.6 mm on the printed page. Journal of Herpetology requires lettering of
// 1.5-2.0 mm after reduction, so figures are authored at final print width and are
// not reduced. Output width must never exceed 215 mm (21.5 cm).
//
//   const F = require('./figlib');
//   const f = F.fig({ wmm: 175, hmm: 90 });
//   f.rect(...); f.text(...);
//   await F.render(f, 'out.png', { dpi: 600 });

const fs = require('fs');
const path = require('path');
const sharp = require(path.join(__dirname, '..', 'node_modules', 'sharp'));

// ---------------------------------------------------------------- constants

const MAX_WMM = 215;          // journal hard limit, 21.5 cm
const MAX_HMM = 280;          // journal hard limit, 28 cm

const FONT = 'Helvetica, Arial, Liberation Sans, sans-serif';

// Okabe-Ito, colourblind-safe. CLEAN/MUCK are the paper's two rearing treatments.
const C = {
  CLEAN: '#0072B2',
  MUCK: '#D55E00',
  FETUS: '#009E73',
  NEWBORN: '#CC79A7',
  BLANK: '#7F7F7F',
  INK: '#1A1A1A',
  MUTED: '#5A5A5A',
  GRID: '#D5D5D5',
  RULE: '#1A1A1A',
};

// Type scale, in 0.1 mm. Nothing readable may go below 15 (1.5 mm).
const T = {
  tick: 16,        // 1.6 mm
  axis: 19,        // 1.9 mm
  panel: 26,       // 2.6 mm, bold — panel letters A, B, C
  legend: 17,      // 1.7 mm
  annot: 16,       // 1.6 mm
  title: 20,       // 2.0 mm
};

// Stroke widths, in 0.1 mm.
const S = { axis: 2.6, grid: 1.2, marker: 1.6, bar: 1.4, whisker: 2.0, median: 3.2 };

const esc = (s) => String(s)
  .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
  .replace(/"/g, '&quot;');

const n = (v) => (Math.round(v * 100) / 100);

// ---------------------------------------------------------------- figure

function fig({ wmm, hmm }) {
  if (wmm > MAX_WMM) throw new Error(`width ${wmm} mm exceeds the ${MAX_WMM} mm journal limit`);
  if (hmm > MAX_HMM) throw new Error(`height ${hmm} mm exceeds the ${MAX_HMM} mm journal limit`);
  const W = wmm * 10, H = hmm * 10;
  const parts = [];

  const o = {
    wmm, hmm, W, H, parts,
    raw: (s) => { parts.push(s); return o; },

    rect(x, y, w, h, { fill = 'none', stroke = 'none', sw = S.bar, opacity = 1 } = {}) {
      parts.push(`<rect x="${n(x)}" y="${n(y)}" width="${n(w)}" height="${n(h)}" fill="${fill}" fill-opacity="${opacity}" stroke="${stroke}" stroke-width="${sw}"/>`);
      return o;
    },

    line(x1, y1, x2, y2, { stroke = C.INK, sw = S.axis, dash = null } = {}) {
      parts.push(`<line x1="${n(x1)}" y1="${n(y1)}" x2="${n(x2)}" y2="${n(y2)}" stroke="${stroke}" stroke-width="${sw}"${dash ? ` stroke-dasharray="${dash}"` : ''} stroke-linecap="round"/>`);
      return o;
    },

    circle(cx, cy, r, { fill = C.INK, stroke = 'none', sw = S.marker, opacity = 1 } = {}) {
      parts.push(`<circle cx="${n(cx)}" cy="${n(cy)}" r="${n(r)}" fill="${fill}" fill-opacity="${opacity}" stroke="${stroke}" stroke-width="${sw}"/>`);
      return o;
    },

    // anchor: start | middle | end ; baseline handled by caller via y
    text(x, y, s, { size = T.tick, fill = C.INK, anchor = 'start', weight = 400, italic = false, rotate = 0, family = FONT } = {}) {
      const tf = rotate ? ` transform="rotate(${rotate} ${n(x)} ${n(y)})"` : '';
      parts.push(`<text x="${n(x)}" y="${n(y)}" font-family="${family}" font-size="${size}" fill="${fill}" text-anchor="${anchor}" font-weight="${weight}"${italic ? ' font-style="italic"' : ''}${tf}>${s}</text>`);
      return o;
    },

    // Panel letter in the journal's required capital-letter form.
    panelLetter(x, y, letter) {
      return o.text(x, y, esc(letter), { size: T.panel, weight: 700, fill: C.INK });
    },

    svg() {
      return `<svg xmlns="http://www.w3.org/2000/svg" width="__PXW__" height="__PXH__" viewBox="0 0 ${W} ${H}" font-family="${FONT}">`
        + `<rect width="${W}" height="${H}" fill="#ffffff"/>`
        + parts.join('')
        + '</svg>';
    },
  };
  return o;
}

// ---------------------------------------------------------------- scales & axes

// Nice tick locator: returns ticks covering [lo, hi] with roughly `target` steps.
function ticks(lo, hi, target = 5) {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) {
    out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  }
  return out;
}

function fmt(v, step) {
  const d = step === undefined ? 0 : Math.max(0, -Math.floor(Math.log10(step)));
  const s = v.toFixed(Math.min(6, d));
  return s.replace(/^-/, '−'); // real minus sign
}

// A rectangular plotting area with linear scales and drawn axes.
//   box = { x, y, w, h }  in viewBox units (0.1 mm)
function axes(f, box, opts = {}) {
  const {
    xlim, ylim, xlabel = null, ylabel = null,
    xticks = null, yticks = null,
    xtickLabels = null, ytickLabels = null,
    grid = false, nxt = 5, nyt = 5,
    xtickRotate = 0, drawX = true, drawY = true,
    categorical = false,
  } = opts;

  const [x0, x1] = xlim, [y0, y1] = ylim;
  const sx = (v) => box.x + ((v - x0) / (x1 - x0)) * box.w;
  const sy = (v) => box.y + box.h - ((v - y0) / (y1 - y0)) * box.h;

  const xt = xticks || (categorical ? [] : ticks(x0, x1, nxt));
  const yt = yticks || ticks(y0, y1, nyt);
  const xstep = xt.length > 1 ? Math.abs(xt[1] - xt[0]) : undefined;
  const ystep = yt.length > 1 ? Math.abs(yt[1] - yt[0]) : undefined;

  if (grid) {
    for (const v of yt) if (v >= y0 && v <= y1) f.line(box.x, sy(v), box.x + box.w, sy(v), { stroke: C.GRID, sw: S.grid });
  }

  // y axis
  if (drawY) {
    f.line(box.x, box.y, box.x, box.y + box.h, { stroke: C.RULE, sw: S.axis });
    yt.forEach((v, i) => {
      if (v < y0 - 1e-9 || v > y1 + 1e-9) return;
      const yy = sy(v);
      f.line(box.x - 8, yy, box.x, yy, { stroke: C.RULE, sw: S.axis });
      const lab = ytickLabels ? ytickLabels[i] : fmt(v, ystep);
      f.text(box.x - 13, yy + T.tick * 0.35, esc(lab), { size: T.tick, anchor: 'end' });
    });
    if (ylabel) {
      f.text(box.x - 74, box.y + box.h / 2, ylabel, { size: T.axis, anchor: 'middle', rotate: -90 });
    }
  }

  // x axis
  if (drawX) {
    f.line(box.x, box.y + box.h, box.x + box.w, box.y + box.h, { stroke: C.RULE, sw: S.axis });
    xt.forEach((v, i) => {
      if (v < x0 - 1e-9 || v > x1 + 1e-9) return;
      const xx = sx(v);
      f.line(xx, box.y + box.h, xx, box.y + box.h + 8, { stroke: C.RULE, sw: S.axis });
      const lab = xtickLabels ? xtickLabels[i] : fmt(v, xstep);
      if (xtickRotate) {
        f.text(xx, box.y + box.h + 16, esc(lab), { size: T.tick, anchor: 'end', rotate: xtickRotate });
      } else {
        f.text(xx, box.y + box.h + 16 + T.tick * 0.8, esc(lab), { size: T.tick, anchor: 'middle' });
      }
    });
    if (xlabel) {
      const dy = xtickRotate ? 88 : 46;
      f.text(box.x + box.w / 2, box.y + box.h + dy + T.axis, xlabel, { size: T.axis, anchor: 'middle' });
    }
  }

  return { sx, sy, box, xt, yt };
}

// ---------------------------------------------------------------- marks

// Deterministic jitter so repeat runs are byte-identical.
function jitter(i, amp) {
  const s = Math.sin(i * 12.9898) * 43758.5453;
  return (s - Math.floor(s) - 0.5) * 2 * amp;
}

// Box-and-whisker from raw values. Returns the five-number summary.
function boxplot(f, ax, xCentre, values, { width = 90, color = C.INK, fill = null } = {}) {
  const v = [...values].sort((a, b) => a - b);
  const q = (p) => {
    const idx = (v.length - 1) * p, lo = Math.floor(idx), hi = Math.ceil(idx);
    return v[lo] + (v[hi] - v[lo]) * (idx - lo);
  };
  const q1 = q(0.25), med = q(0.5), q3 = q(0.75);
  const iqr = q3 - q1;
  const lo = v.find((x) => x >= q1 - 1.5 * iqr);
  const hi = [...v].reverse().find((x) => x <= q3 + 1.5 * iqr);
  const half = width / 2;

  f.line(xCentre, ax.sy(lo), xCentre, ax.sy(hi), { stroke: color, sw: S.whisker });
  f.line(xCentre - half * 0.45, ax.sy(lo), xCentre + half * 0.45, ax.sy(lo), { stroke: color, sw: S.whisker });
  f.line(xCentre - half * 0.45, ax.sy(hi), xCentre + half * 0.45, ax.sy(hi), { stroke: color, sw: S.whisker });
  f.rect(xCentre - half, ax.sy(q3), width, ax.sy(q1) - ax.sy(q3),
    { fill: fill || color, opacity: fill ? 1 : 0.18, stroke: color, sw: S.bar });
  f.line(xCentre - half, ax.sy(med), xCentre + half, ax.sy(med), { stroke: color, sw: S.median });
  return { q1, med, q3, lo, hi };
}

// Jittered strip of points over a category centre.
function strip(f, ax, xCentre, values, { color = C.INK, r = 7, amp = 26 } = {}) {
  values.forEach((v, i) => {
    f.circle(xCentre + jitter(i, amp), ax.sy(v), r, { fill: color, stroke: '#ffffff', sw: 1.2 });
  });
}

// ---------------------------------------------------------------- data loading

const RES = path.resolve(__dirname, '..');

function tsv(file, { comment = '#q2:types' } = {}) {
  const raw = fs.readFileSync(file, 'utf8').replace(/\s+$/, '');
  const lines = raw.split(/\r?\n/).filter((l) => l !== '' && !l.startsWith(comment));
  const head = lines[0].split('\t');
  return lines.slice(1).map((l) => {
    const p = l.split('\t');
    const o = {};
    head.forEach((h, i) => { o[h] = p[i]; });
    return o;
  });
}

// The 26 reared animals with treatment, keyed by sample id.
// NOTE the metadata's `environment` column lumps the 3 Fetus samples in with the
// clean-reared ("Sterile" = 19). `primary-comparison` is the correct column:
// Sterile = the 16 clean-reared animals only.
function rearedGroups() {
  const rows = tsv(path.join(RES, 'scripts', 'qzv', 'alpha_shannon.tsv'));
  const g = {};
  for (const r of rows) {
    const pc = r['primary-comparison'];
    if (pc === 'Sterile') g[r.id] = 'Clean';
    else if (pc === 'Muck') g[r.id] = 'Muck';
  }
  return g;
}

function alpha(metric) {
  const map = {
    shannon: ['alpha_shannon.tsv', 'shannon_entropy'],
    observed: ['alpha_observed_features.tsv', 'observed_features'],
    faith: ['alpha_faith_pd.tsv', 'faith_pd'],
    evenness: ['alpha_evenness.tsv', 'pielou_evenness'],
  }[metric];
  if (!map) throw new Error('unknown alpha metric ' + metric);
  const rows = tsv(path.join(RES, 'scripts', 'qzv', map[0]));
  const out = {};
  for (const r of rows) out[r.id] = +r[map[1]];
  return out;
}

// ---------------------------------------------------------------- render

async function render(f, outPath, { dpi = 600 } = {}) {
  const pxw = Math.round(f.wmm / 25.4 * dpi);
  const pxh = Math.round(f.hmm / 25.4 * dpi);
  const svg = f.svg().replace('__PXW__', pxw).replace('__PXH__', pxh);
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  const info = await sharp(Buffer.from(svg))
    .png({ compressionLevel: 9 })
    .withMetadata({ density: dpi })
    .toFile(outPath);
  const wcm = (info.width / dpi) * 2.54, hcm = (info.height / dpi) * 2.54;
  console.log(`  ${path.basename(outPath)}  ${info.width}x${info.height} px  ${dpi} dpi  ${wcm.toFixed(1)} x ${hcm.toFixed(1)} cm`);
  if (wcm > 21.5001 || hcm > 28.0001) throw new Error(`${outPath} exceeds the journal size limit`);
  return info;
}

module.exports = {
  fig, axes, ticks, fmt, boxplot, strip, jitter, render, tsv, esc, n,
  rearedGroups, alpha,
  C, T, S, FONT, RES, MAX_WMM, MAX_HMM,
};
