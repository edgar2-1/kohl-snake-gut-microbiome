// Produce a deposit-ready sample metadata file for SRA/ENA submission (panel item B12).
//   node make_deposit_metadata.js
//
// Two defects in meta/metadata.tsv would otherwise be published verbatim:
//   1. The control group is labelled "Sterile" although the Methods explicitly disclaim
//      sterility ("Clean cages were washed, not sterilized, and neither treatment was axenic").
//   2. The `environment` column codes the three Fetus samples as "Sterile", so a downloader
//      counting that column gets n = 19 controls instead of 16 — the origin of a historical
//      error in this project.
// The working file is left untouched so the analysis scripts keep running.
const fs = require('fs'), path = require('path');
const RES = path.resolve(__dirname, '..');
const SRC = path.join(RES, 'meta', 'metadata.tsv');
const DST = path.join(RES, 'meta', 'metadata_for_deposit.tsv');

const lines = fs.readFileSync(SRC, 'utf8').replace(/\s+$/, '').split(/\r?\n/);
const H = lines[0].split('\t');
const iEnv = H.indexOf('environment'), iPri = H.indexOf('primary-comparison'), iCls = H.indexOf('sample-class');
if (iEnv < 0 || iPri < 0) throw new Error('expected environment and primary-comparison columns');

const changes = { env: 0, pri: 0, fetus: 0 };
const out = [lines[0]];
for (let i = 1; i < lines.length; i++) {
  const p = lines[i].split('\t');
  if (p[0] === '#q2:types') { out.push(lines[i]); continue; }

  const pri = p[iPri];
  // Fetus samples are mis-coded as "Sterile" in `environment`; recode from primary-comparison.
  if (p[iEnv] === 'Sterile' && pri === 'Fetus') { p[iEnv] = 'Fetus'; changes.fetus++; }
  else if (p[iEnv] === 'Sterile') { p[iEnv] = 'Clean'; changes.env++; }
  if (p[iPri] === 'Sterile') { p[iPri] = 'Clean'; changes.pri++; }
  out.push(p.join('\t'));
}
fs.writeFileSync(DST, out.join('\n') + '\n', 'utf8');

// verify
const chk = fs.readFileSync(DST, 'utf8').replace(/\s+$/, '').split(/\r?\n/);
const count = (idx) => {
  const c = {};
  for (let i = 1; i < chk.length; i++) { const p = chk[i].split('\t'); if (p[0] === '#q2:types') continue; c[p[idx] || '(empty)'] = (c[p[idx] || '(empty)'] || 0) + 1; }
  return c;
};
console.log('wrote meta/metadata_for_deposit.tsv');
console.log(`  environment: Sterile -> Clean on ${changes.env} rows; Sterile -> Fetus on ${changes.fetus} rows`);
console.log(`  primary-comparison: Sterile -> Clean on ${changes.pri} rows`);
console.log('  environment counts now:        ' + JSON.stringify(count(iEnv)));
console.log('  primary-comparison counts now: ' + JSON.stringify(count(iPri)));
const remaining = fs.readFileSync(DST, 'utf8').match(/Sterile/g);
console.log('  remaining occurrences of "Sterile": ' + (remaining ? remaining.length : 0));
