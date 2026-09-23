// Rebuild scripts/qzv/ (distance matrices, alpha vectors, genus table) directly from
// the QIIME 2 artifacts in results/qiime_artifacts/. The original export directory was
// lost; this makes reanalysis_node_run.js and host_trait_maternal.js runnable again.
//   node rebuild_qzv_inputs.js
const fs = require('fs'), path = require('path'), zlib = require('zlib');

const QA = path.resolve(__dirname, '..', 'qiime_artifacts');
const OUT = path.join(__dirname, 'qzv');

function zipEntries(file) {
  const b = fs.readFileSync(file);
  let eocd = -1;
  for (let i = b.length - 22; i >= 0 && i > b.length - 66000; i--) {
    if (b.readUInt32LE(i) === 0x06054b50) { eocd = i; break; }
  }
  if (eocd < 0) throw new Error('no EOCD: ' + file);
  const n = b.readUInt16LE(eocd + 10);
  let off = b.readUInt32LE(eocd + 16);
  const list = [];
  for (let k = 0; k < n; k++) {
    if (b.readUInt32LE(off) !== 0x02014b50) break;
    const e = {
      method: b.readUInt16LE(off + 10),
      csize: b.readUInt32LE(off + 20),
      nameLen: b.readUInt16LE(off + 28),
      extraLen: b.readUInt16LE(off + 30),
      cmtLen: b.readUInt16LE(off + 32),
      lho: b.readUInt32LE(off + 42),
    };
    e.name = b.slice(off + 46, off + 46 + e.nameLen).toString('utf8');
    list.push(e);
    off += 46 + e.nameLen + e.extraLen + e.cmtLen;
  }
  return {
    read(suffix) {
      const e = list.find((x) => x.name.endsWith(suffix));
      if (!e) throw new Error('entry not found: ' + suffix + ' in ' + path.basename(file));
      const lnl = b.readUInt16LE(e.lho + 26), lel = b.readUInt16LE(e.lho + 28);
      const raw = b.slice(e.lho + 30 + lnl + lel, e.lho + 30 + lnl + lel + e.csize);
      return (e.method === 0 ? raw : zlib.inflateRawSync(raw)).toString('utf8');
    },
  };
}

fs.mkdirSync(OUT, { recursive: true });
const report = [];

// ---- distance matrices (long form: SubjectID1 / SubjectID2 / Distance) ----
const METRICS = ['bray_curtis', 'jaccard', 'unweighted_unifrac', 'weighted_unifrac'];
for (const m of METRICS) {
  const src = path.join(QA, 'core-metrics', `${m}_distance_matrix-environment-permanova.qzv`);
  const txt = zipEntries(src).read('data/raw_data.tsv');
  const dst = path.join(OUT, `dist_${m}.tsv`);
  fs.writeFileSync(dst, txt, 'utf8');
  const lines = txt.replace(/\s+$/, '').split(/\r?\n/);
  const ids = new Set();
  const h = lines[0].split('\t');
  const i1 = h.indexOf('SubjectID1'), i2 = h.indexOf('SubjectID2');
  for (let i = 1; i < lines.length; i++) { const p = lines[i].split('\t'); ids.add(p[i1]); ids.add(p[i2]); }
  report.push(`dist_${m}.tsv        ${String(lines.length - 1).padStart(5)} pairs, ${ids.size} samples`);
}

// ---- alpha vectors ----
const ALPHA = { shannon: 'shannon', observed_features: 'observed_features', faith_pd: 'faith_pd', evenness: 'evenness' };
for (const [out, vec] of Object.entries(ALPHA)) {
  const src = path.join(QA, 'core-metrics', `${vec}_vector-group-significance.qzv`);
  const txt = zipEntries(src).read('data/metadata.tsv');
  const dst = path.join(OUT, `alpha_${out}.tsv`);
  fs.writeFileSync(dst, txt, 'utf8');
  const lines = txt.trim().split(/\r?\n/).filter((l) => !l.startsWith('#q2:types'));
  const hdr = lines[0].split('\t');
  report.push(`alpha_${out}.tsv`.padEnd(28) + `${lines.length - 1} samples, value column "${hdr[hdr.length - 1]}"`);
}

// ---- genus (level-6) table ----
const l6 = zipEntries(path.join(QA, 'taxa-bar-plots.qzv')).read('data/level-6.csv');
fs.writeFileSync(path.join(OUT, 'level-6.csv'), l6, 'utf8');
report.push(`level-6.csv                 ${l6.trim().split(/\r?\n/).length - 1} samples`);

console.log('Rebuilt ' + OUT + '\n');
report.forEach((r) => console.log('  ' + r));
