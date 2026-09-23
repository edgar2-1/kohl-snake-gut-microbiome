// Export the ASV feature table from a QIIME 2 .qza (HDF5 BIOM v2.1) to TSV,
// without QIIME/Python. Run:  node export_asv_table.mjs <in.qza> <out.tsv>
import * as h5 from 'h5wasm/node';
import fs from 'fs';
import path from 'path';
import os from 'os';
import { execFileSync } from 'child_process';

const QZA = process.argv[2];
const OUT = process.argv[3];

// unzip the .qza and pull out the .biom payload
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'qza-'));
execFileSync('powershell', ['-NoProfile', '-Command',
  `Add-Type -AssemblyName System.IO.Compression.FileSystem; ` +
  `$z=[System.IO.Compression.ZipFile]::OpenRead('${QZA}'); ` +
  `$e=$z.Entries | Where-Object { $_.FullName -like '*/data/feature-table.biom' } | Select-Object -First 1; ` +
  `[System.IO.Compression.ZipFileExtensions]::ExtractToFile($e, '${path.join(tmp,'t.biom').replace(/\\/g,'\\\\')}', $true); $z.Dispose()`
]);

await h5.ready;
const f = new h5.File(path.join(tmp, 't.biom'), 'r');

const obsIds = f.get('observation/ids').value;          // ASV ids
const smpIds = f.get('sample/ids').value;               // sample ids
// BIOM stores CSR by sample: sample/matrix/{data,indices,indptr}
const data = f.get('sample/matrix/data').value;
const indices = f.get('sample/matrix/indices').value;   // observation index
const indptr = f.get('sample/matrix/indptr').value;
f.close();

const nObs = obsIds.length, nSmp = smpIds.length;
// dense matrix: rows = ASVs, cols = samples
const M = Array.from({ length: nObs }, () => new Float64Array(nSmp));
for (let s = 0; s < nSmp; s++) {
  for (let k = Number(indptr[s]); k < Number(indptr[s + 1]); k++) {
    M[Number(indices[k])][s] = Number(data[k]);
  }
}

const lines = ['#OTU ID\t' + smpIds.join('\t')];
for (let i = 0; i < nObs; i++) lines.push(obsIds[i] + '\t' + Array.from(M[i]).join('\t'));
fs.writeFileSync(OUT, lines.join('\n') + '\n');

let nz = 0; for (let i = 0; i < nObs; i++) for (let j = 0; j < nSmp; j++) if (M[i][j] > 0) nz++;
console.log(`wrote ${OUT}: ${nObs} ASVs x ${nSmp} samples, ${nz} non-zero cells`);
console.log('samples:', smpIds.join(', '));
fs.rmSync(tmp, { recursive: true, force: true });
