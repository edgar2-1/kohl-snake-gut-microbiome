// build_storms.js — produce Table S5, the STORMS reporting checklist.
//
//   node scripts/build_storms.js
//
// The item numbers, names and recommendations come from the official STORMS v1.03
// template (Mirzayi et al. 2021, Nat Med 27:1885-1892; editable template at
// https://zenodo.org/records/5703117), vendored as scripts/storms_template_1.03.json
// so the checklist regenerates without network access.
//
// STORMS was written for HUMAN microbiome research. This is a captive-animal study, so
// several participant-facing items are genuinely not applicable; those are answered NA
// with the reason given rather than left blank, which is what the guideline asks for.
//
// Answers marked [[AUTHOR]] cannot be derived from the manuscript and must be filled by
// the first author before submission. They correspond to the same gaps flagged as
// [[MISSING]] in the manuscript itself.

const fs = require('fs');
const path = require('path');
const RES = path.resolve(__dirname, '..');
const template = require(path.join(__dirname, 'storms_template_1.03.json'));

const A = {
  // -------- Abstract
  '1.0': ['Yes', 'Abstract reports background, methods, results and conclusions (244 words).'],
  '1.1': ['Yes', 'Abstract: postnatal rearing experiment, clean vs muck cages, n = 26 neonates from four dams.'],
  '1.2': ['Yes', 'Abstract and Methods (Sequence Processing): 16S rRNA V4 amplicon sequencing, Illumina MiSeq.'],
  '1.3': ['Yes', 'Large-intestinal lumen swab, collected at necropsy. Methods (DNA Extraction, Amplification, and Sequencing).'],

  // -------- Introduction
  '2.0': ['Yes', 'Introduction: vertical vs environmental acquisition of the neonate gut microbiota in squamates.'],
  '2.1': ['Yes', 'Introduction states the objective: test whether postnatal environmental exposure shapes the neonate gut community, and separate that from maternal background. Exploratory with respect to which taxa differ.'],

  // -------- Methods: study design
  '3.0': ['Yes', 'Methods (Study Animals and Sampling): randomized postnatal rearing experiment, two treatments, individually housed 7-19 d, plus two baseline groups (Fetus, Newborn) characterizing the starting community.'],
  '3.1': ['Yes', 'Adapted: subjects are captive-born Nerodia rhombifer neonates, not human participants. 26 reared neonates from four dams; Table 1 gives the allocation by dam and treatment.'],
  '3.2': ['[[AUTHOR]]', 'Holding facility and the source pond for the muck inoculum are recorded as [[MISSING]] in Methods. Supply the institution and the general locality.'],
  '3.3': ['Partly', 'Births in 2019; dams 11.21 and 11.22 gave birth 21 and 30 August 2019. Parturition dates for NR1908, NR1910 and NR1912 were not recorded (Table 1 legend). Per-animal birth and sampling dates are in Table S4.'],
  '3.4': ['Yes', 'Methods and Table 1: all liveborn neonates of the four dams that survived to sampling were included. Dam NR1910 contributed only the three Fetus samples.'],
  '3.5': ['Yes', 'No antibiotics were administered. Animals were captive-born and held 7-19 d before sampling.'],
  '3.6': ['Yes', 'Methods and Table 1: 37 libraries sequenced (32 biological, 4 no-template controls, 1 extraction control); 36-library feature table after decontamination; 26 reared animals in the primary contrast; 30 samples retained after rarefaction to 5,000 reads.'],
  '3.7': ['NA', 'Not longitudinal. Each animal was sampled once, terminally.'],
  '3.8': ['NA', 'Not a matched design. Allocation to treatment was randomized within dam where litter size allowed; dam 11.21 was the only brood with balanced allocation (5 Clean / 5 Muck).'],
  '3.9': ['[[AUTHOR]]', 'Methods cites the HACC (2004) guidelines and AVMA (2020) euthanasia guidelines, but the IACUC protocol number, approving institution, and collection/holding permit numbers are [[MISSING]].'],

  // -------- Methods: laboratory
  '4.0': ['[[AUTHOR]]', 'Extraction laboratory and sequencing centre are [[MISSING]] in Methods.'],
  '4.1': ['Yes', 'Methods: swab of the large-intestinal lumen taken at necropsy. Sampling was terminal because the lumen cannot be sampled regionally in vivo at this body size.'],
  '4.2': ['NA', 'Samples were not shipped between collection and extraction; they were frozen in place.'],
  '4.3': ['Yes', 'Methods: swabs frozen at -80 degrees C until DNA extraction.'],
  '4.4': ['[[AUTHOR]]', 'Extraction kit and version are [[MISSING]]. Post-extraction concentrations are reported (7.5-53.6 ng/uL for the 32 biological samples; 1.4 ng/uL for the extraction control).'],
  '4.5': ['NA', 'No host-DNA depletion or microbial enrichment was performed; amplicon sequencing was used.'],
  '4.6': ['Yes', 'Methods: 16S rRNA V4 region, primers 515F (GTGYCAGCMGCCGCGGTAA; Parada et al. 2016) and 806R (GGACTACNVGGGTWTCTAAT; Apprill et al. 2015).'],
  '4.7': ['No', 'No mock community was included. Noted as a limitation.'],
  '4.8': ['Yes', 'Four no-template PCR negative controls carried through library preparation and sequencing, plus a reagents-only extraction control (SS120). Read yields for all five are in Table S4.'],
  '4.9': ['Yes', 'Methods (Sequence Processing): prevalence and frequency screens at score threshold 0.1 following decontam (Davis et al. 2018); differentially abundant genera additionally screened against the Salter et al. (2014) reagent-contaminant catalog. Results report the extraction control (141 reads in five ASVs, three absent from every biological sample) and flag which Clean-enriched genera are recurrent reagent taxa.'],
  '4.10': ['No', 'No technical or biological replicates were sequenced. The pipeline was, however, re-run with the extraction control included (all 37 libraries denoised together) as a validation; it recovered 1,050 of the 1,052 analyzed ASVs.'],
  '4.11': ['Yes', 'Amplicon (16S rRNA V4) sequencing. Methods.'],
  '4.12': ['Yes', 'Relative abundance only; no quantitative (QMP, cell-count or spike-in) profiling was performed.'],
  '4.13': ['Yes', 'All 37 libraries were prepared and sequenced on a single Illumina MiSeq run (2 x 251 bp), so batch is not confounded with treatment. Methods.'],
  '4.14': ['NA', 'No metatranscriptomics.'],
  '4.15': ['NA', 'No metaproteomics.'],
  '4.16': ['NA', 'No metabolomics.'],

  // -------- Methods: variables and design
  '5.0': ['Yes', 'Methods (Study Animals and Sampling): birth date, treatment, exposure duration, feeding, sex, SVL, total length, mass and age at sampling recorded per animal; body condition is the residual of log mass on log SVL. Per-animal values in Table S4.'],
  '6.0': ["Yes", "Methods: dam is the principal confounder. Every treatment contrast compares siblings only with siblings, evaluated under all 211,680 within-litter treatment assignments, so litter is built into the test; the same test within the one evenly split litter, whose animals were all 9 d old and unfed, gives the contrast free of age and diet. Feeding is partly confounded with dam and treatment (Table 1)."],
  '6.1': ['Yes', 'Methods and Results: two of six early-timepoint libraries fell below the rarefaction threshold, so the Fetus and Newborn groups are analyzed qualitatively and separately. No animals were lost to follow-up.'],

  // -------- Methods: bioinformatics and statistics
  '7.0': ["Yes", "Methods (Statistical Analyses): diversity analyses on a table rarefied to 5,000 reads per sample; genus-level differential abundance on relative abundances from the unrarefied filtered table."],
  '7.1': ['Yes', 'Methods (Sequence Processing): cutadapt primer removal; DADA2 denoising with forward and reverse reads truncated to 220 and 180 bp; chimera filtering; removal of ASVs unassigned at phylum level and those classified as mitochondria, chloroplast or Eukaryota. Per-library retention is Fig. 1 and Table S4.'],
  '7.2': ['Yes', 'Methods: consensus VSEARCH classification in q2-feature-classifier against SILVA 138.2 curated with RESCRIPt. No functional profiling is reported (the PICRUSt2 analysis was withdrawn before submission).'],
  '7.3': ["Yes", "Methods (Statistical Analyses): one test for every Clean-versus-Muck contrast, evaluated under all 211,680 within-litter treatment assignments so P is exact: PERMANOVA pseudo-F for composition, and the stratified rank-sum statistic (sum of within-litter Mann-Whitney U, after van Elteren, 1960, with equal litter weights) for alpha diversity and genus-level abundance. PERMDISP on the spatial median (9,999 permutations) for dispersion; PERMANOVA with dam as a four-level factor (9,999 permutations) for maternal identity; Benjamini-Hochberg FDR across the four alpha metrics and across the 328 features. All tests two-sided."],
  '7.4': ['NA', 'Not longitudinal.'],
  '7.5': ["Yes", "Methods and Table S2: the same treatment test applied within each litter that received both treatments (exact P), with alpha-diversity medians by treatment within those litters."],
  '7.6': ['Yes', 'Parturition dates for three dams were not recorded (Table 1 legend); cage dimensions, water regime, temperature and photoperiod were not recorded (Methods, Captive Conditions). No imputation was performed; affected animals are excluded only from analyses requiring the missing variable.'],
  '7.7': ["NA", "None beyond the within-litter contrasts in Table S2; the primary test conditions on litter by design."],
  '7.8': ['Yes', 'Methods: features reported at q <= 0.05 in Table 3 and Fig. 6; the full q <= 0.10 set plus five below-threshold taxa discussed in the text is in Table S1. Aeromonas prevalence is a post-selection confirmatory contrast and is identified as such.'],
  '7.9': ['Yes', 'QIIME 2 v.2026.4, cutadapt, DADA2, VSEARCH, q2-feature-classifier, SILVA 138.2 with RESCRIPt, MAFFT, FastTree, R 4.5.2 with vegan 2.7-3 (adonis2, betadisper); custom Node.js scripts for the exact within-litter tests.'],

  // -------- Methods: access
  '8.0': ["Yes", "Methods (Data Availability): analysis code (scripts/exact_stratified_tests.js, scripts/permdisp_validated.js), exported distance matrices, diversity vectors and intermediate result files are in the project repository. The treatment tests enumerate every within-litter assignment and use no random numbers; the dam test and PERMDISP use a seeded generator; pseudo-F and PERMDISP statistics were validated against QIIME 2 and vegan."],
  '8.1': ['[[AUTHOR]]', 'All 37 libraries to be deposited in the NCBI SRA; the BioProject accession is [[MISSING]] and must exist before submission.'],
  '8.2': ['Yes', 'Feature tables, distance matrices and diversity vectors are in the repository at https://github.com/edgar2-1/kohl-snake-gut-microbiome'],
  '8.3': ['Yes', 'Adapted: per-animal metadata (dam, treatment, exposure, feeding, sex, morphometrics, sequencing depth) for all 37 libraries is Table S4. No human participants, so no de-identification is required.'],
  '8.4': ['[[AUTHOR]]', 'Code is in the repository above; the archived Zenodo DOI is [[MISSING]].'],
  '8.5': ["Yes", "Table S1 gives the genus-level results for all 328 features tested; Table S2 the composition, dispersion and dam analyses on all four metrics and the within-litter contrasts; Table S3 the early-timepoint ASV counts. All are computer-readable TSV."],

  // -------- Results
  '9.0': ['Yes', 'Table 1: study animals by dam and treatment, with SVL and mass medians, exposure durations, feeding and sequencing depth. Table S4 gives per-animal values.'],
  '10.0': ["Yes", "Results: alpha diversity (Fig. 4), beta diversity and PERMANOVA (Fig. 3, Table 2), taxonomic composition (Fig. 2), and maternal identity (Table 2)."],
  '10.1': ['Yes', 'SILVA 138.2 taxonomy; full taxonomy strings for every reported feature are in Table S1. Features unresolved at genus are named by their finest resolved rank.'],
  '10.2': ["Yes", "Table 3 and Fig. 6 give the 14 genus-level features at q <= 0.05 with group means, carriers, log2 fold-change and q; Table S1 reports all 328 features tested."],
  '10.3': ['NA', 'No functional potential, MAG assembly or RNA-seq is reported. The PICRUSt2 predicted-function analysis was withdrawn before submission and is not referenced.'],
  '10.4': ["Yes", "PERMDISP; the early-timepoint ASV membership comparison (Table S3)."],

  // -------- Discussion
  '11.0': ['Yes', 'Discussion opens with the environment effect and its size relative to maternal identity.'],
  '12.0': ['Yes', 'Discussion: the effect is interpreted as environmental acquisition without claiming establishment; detection is distinguished from colonization throughout.'],
  '13.0': ['Yes', 'Discussion (limits paragraph): one cohort and one source pond; feeding and exposure duration partly confounded with litter; no substrate, water or prey sample sequenced; 16S V4 does not resolve species or strain. The low-biomass early timepoints are qualified in Results.'],
  '13.1': ["Yes", "Results report the extraction control and identify which Clean-enriched genera are recurrent reagent taxa; the dam confound is addressed by the within-litter design of the primary test and by the dam model in Table 2."],
  '13.2': ['Yes', 'Discussion (limits paragraph): the animals came from one cohort and one source pond; no claim is made about wild neonates or other species.'],
  '14.0': ['Yes', 'Discussion identifies replication of the exposure across litters and substrate sources as the next step, and notes that non-terminal cloacal sampling would enable longitudinal designs. Establishing the pond as the source will require a new study that samples the rearing material alongside the animals: no environmental or dietary samples were retained from this cohort (Limitations).'],

  // -------- Other information
  '15.0': ['[[AUTHOR]]', 'Funding sources and the role of funders are [[MISSING]] in the Acknowledgments.'],
  '15.1': ['Yes', 'Acknowledgments, which include the dedication to the late S. M. Secor. People to thank are still to be named ([[MISSING]] in the manuscript).'],
  '15.2': ['[[AUTHOR]]', 'A conflicts-of-interest statement must be added before submission.'],
  '16.0': ["Yes", "Supplementary Material comprises Table S1 (full differential abundance), Table S2 (all-metric composition, dispersion and dam analyses, and within-litter contrasts), Table S3 (early-timepoint ASV overlap), Table S4 (per-sample metadata) and Table S5 (this checklist). Hosted by BioOne with the article."],
  '17.0': ["Yes", "Tables S1-S4 are computer-readable TSV files. Table S1 reports every one of the 328 features tested with full SILVA taxonomy strings, exact P and BH q; Table S2 reports the environment, dispersion and dam analyses on all four distance metrics and the within-litter contrasts."],
};

const HEAD = ['Number', 'Item', 'Recommendation', 'Section', 'Yes/No/NA', 'Comments or location in manuscript'];

const out = [
  '# Table S5. STORMS reporting checklist (Mirzayi et al., 2021; template v1.03).',
  '# Adapted for a captive-animal study: STORMS was written for human microbiome research,',
  '# so participant-facing items are answered NA with the reason. Entries marked [[AUTHOR]]',
  '# require information not present in the manuscript and must be completed before submission.',
  HEAD.join('\t'),
];

let section = '';
let filled = 0, todo = 0, missingItems = [];
for (const row of template) {
  const [num, item, rec] = row;
  if (!num || num === 'Number') continue;
  if (!/^\d+\.\d+$/.test(num)) { section = num; continue; }   // section header row
  const ans = A[num];
  if (!ans) { missingItems.push(num); continue; }
  if (ans[0] === '[[AUTHOR]]') todo++; else filled++;
  out.push([num, item, rec, section, ans[0], ans[1]].map((c) => String(c).replace(/\t/g, ' ')).join('\t'));
}

if (missingItems.length) {
  console.error('ERROR: no answer authored for items: ' + missingItems.join(', '));
  process.exit(1);
}

const dst = path.join(RES, 'figures', 'submission', 'TableS5.tsv');
fs.writeFileSync(dst, out.join('\n') + '\n', 'utf8');
console.log(`wrote figures/submission/TableS5.tsv`);
console.log(`  ${filled + todo} items: ${filled} answered, ${todo} awaiting author input`);
console.log('  awaiting author: ' + Object.entries(A).filter(([, v]) => v[0] === '[[AUTHOR]]').map(([k]) => k).join(', '));
