const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
        AlignmentType, LevelFormat, HeadingLevel, BorderStyle, WidthType,
        ShadingType, ImageRun, PageNumber, Footer } = require("docx");

// ---------- helpers ----------
const FONT = "Calibri";
const body = (runs, opts = {}) => new Paragraph({
  alignment: AlignmentType.JUSTIFIED,
  spacing: { after: 160, line: 276 },
  ...opts,
  children: Array.isArray(runs) ? runs : [new TextRun(runs)],
});
const t = (text, o = {}) => new TextRun({ text, ...o });
const b = (text) => new TextRun({ text, bold: true });
const i = (text) => new TextRun({ text, italics: true });
const H1 = (text) => new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 280, after: 140 },
  children: [new TextRun({ text, bold: true })] });
const H2 = (text) => new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 200, after: 100 },
  children: [new TextRun({ text, bold: true })] });

const img = (file, wPx, hPx, caption) => {
  const para = new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 },
    children: [ new ImageRun({ type: "png", data: fs.readFileSync(file),
      transformation: { width: wPx, height: hPx } }) ] });
  const cap = new Paragraph({ alignment: AlignmentType.LEFT, spacing: { after: 200 },
    children: caption });
  return [para, cap];
};

// table helpers
const BORD = { style: BorderStyle.SINGLE, size: 1, color: "BBBBBB" };
const BORDERS = { top: BORD, bottom: BORD, left: BORD, right: BORD };
const cell = (children, w, opts = {}) => new TableCell({
  borders: BORDERS, width: { size: w, type: WidthType.DXA },
  margins: { top: 60, bottom: 60, left: 100, right: 100 },
  shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
  children: children.map(c => typeof c === "string"
    ? new Paragraph({ children: [new TextRun({ text: c, bold: !!opts.bold, size: 18 })],
        alignment: opts.align || AlignmentType.LEFT })
    : c),
});
const headRow = (cells, widths) => new TableRow({ tableHeader: true,
  children: cells.map((c, k) => cell([c], widths[k], { bold: true, fill: "1F4E79",
    align: k === 0 ? AlignmentType.LEFT : AlignmentType.CENTER }))
    .map(tc => { tc.options.children[0] = new Paragraph({ alignment: tc.options.children[0].options.alignment,
      children: [new TextRun({ text: cells[tc.rootKey] , bold: true, color: "FFFFFF", size: 18 })] }); return tc; }) });

// simpler header row (white text on dark)
function tableHeader(cells, widths) {
  return new TableRow({ tableHeader: true, children: cells.map((txt, k) =>
    new TableCell({ borders: BORDERS, width: { size: widths[k], type: WidthType.DXA },
      margins: { top: 60, bottom: 60, left: 100, right: 100 },
      shading: { fill: "1F4E79", type: ShadingType.CLEAR },
      children: [ new Paragraph({ alignment: k === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
        children: [ new TextRun({ text: txt, bold: true, color: "FFFFFF", size: 18 }) ] }) ] })) });
}
function row(cells, widths, opts = {}) {
  return new TableRow({ children: cells.map((txt, k) =>
    new TableCell({ borders: BORDERS, width: { size: widths[k], type: WidthType.DXA },
      margins: { top: 60, bottom: 60, left: 100, right: 100 },
      shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
      children: [ new Paragraph({ alignment: k === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
        children: [ new TextRun({ text: txt, bold: !!opts.bold, size: 18 }) ] }) ] })) });
}

// numbered references
const refItem = (text) => new Paragraph({ numbering: { reference: "refs", level: 0 },
  spacing: { after: 80 }, children: [new TextRun({ text, size: 20 })] });

// ===== CONTENT (model-variance thesis, NetSolP added) =====
const children = [];

// Title
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 120 },
  children: [ new TextRun({ bold: true, size: 32,
    text: "Model- and distribution-dependent performance of protein solubility predictors: two state-of-the-art deep models straddle a simple interpretable baseline on native E. coli proteins" }) ] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
  children: [ t("Jakob Oeschey", { size: 22 }) ] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
  children: [ t("Independent Researcher", { size: 20, italics: true }) ] }));
children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 },
  children: [ t("Correspondence: tiger-biotech@pm.me", { size: 20 }) ] }));

// Abstract
children.push(H1("Abstract"));
children.push(body([
  b("Motivation. "),
  t("Sequence-based predictors of recombinant protein solubility are widely used to triage expression targets, and deep-learning models report strong benchmark AUROCs. Whether a model's headline benchmark predicts how it ranks on an independent distribution, and how much two state-of-the-art models can differ on the same proteins, is rarely tested directly. "),
  b("Results. "),
  t("We benchmarked two independent state-of-the-art deep models, RP3Net (ESM-2 650M, reported AUROC 0.83 on its prospective set) and NetSolP (ESM1b), against two dataset-naive composition baselines (an in-house heuristic and the published Solubility-Weighted Index, SWI) on two "),
  i("Escherichia coli"),
  t(" benchmarks: eSOL (native K-12 proteins) and the SoluProt held-out test set (heterologous constructs). Every test set was deduplicated against each model's public training data with MMseqs2 (eSOL: 0.13% leak to RP3Net, 4.06% to NetSolP/SWI training; SoluProt: 1.29% to RP3Net), and all native-protein comparisons use the doubly-clean intersection. On native cytoplasmic proteins the two deep models spanned the entire range and "),
  b("straddled the simple interpretable baselines"),
  t(": NetSolP was the best predictor (AUROC 0.792) while RP3Net ranked lowest (0.709), a gap of 0.083 (DeLong z = 6.1, p < 1e-8; paired-bootstrap CI +0.057 to +0.108) that exceeds any difference between a model and a baseline. RP3Net, despite its 0.83 headline, was significantly beaten by the simple SWI score (0.745; \u0394 +0.037, p = 0.006) and tied by a naive composition heuristic (0.725), whereas NetSolP significantly beat SWI (\u0394 +0.046, p < 0.001). On the harder heterologous SoluProt distribution the two deep models instead agreed and led numerically (NetSolP 0.633, RP3Net 0.620, difference n.s. p = 0.23); there only NetSolP significantly exceeded the strong SWI baseline (\u0394 +0.035, p \u2248 0.002), while RP3Net exceeded only the weak naive baseline and was statistically indistinguishable from SWI (\u0394 +0.021, p = 0.076), all predictors in a low 0.58\u20130.63 band; a third predictor calibrated on native eSOL (ProteinSol) collapsed to near-chance there (0.542, below every other predictor, p < 0.0001 vs both deep models), the mirror image of RP3Net's weak ranking on native proteins, a distribution outside its recombinant/heterologous design domain. The native ordering was robust to a homology-level leakage screen (MMseqs2 clustering at \u226430% identity removed 8.5% of eSOL as PSI:Biology homologs): on the homology-clean set (n = 1,929 cytoplasmic) NetSolP stayed highest (0.807), RP3Net lowest (0.714), with SWI (0.751) between, and every significant difference persisted (NetSolP vs SWI \u0394 +0.055, p = 4e-5; NetSolP vs RP3Net \u0394 +0.093, p = 3e-11; SWI vs RP3Net \u0394 +0.038, p = 0.009), so the lead is not a residual-training-homology artefact. "),
  b("Conclusion. "),
  t("Solubility-predictor performance is model- and distribution-dependent: two SOTA deep models differed on native proteins by more than the gap to a simple interpretable baseline, and a model's published benchmark AUROC did not predict its generalisation rank. Evaluating RP3Net on native proteins probes generalisation outside its recombinant/heterologous design domain, so its low native ranking reflects distribution mismatch rather than a model deficiency; we therefore frame these results as a caution about distribution-matched evaluation, not a criticism of any model. A dataset-naive baseline is a useful sanity floor: a model falling below it on the target distribution (as RP3Net does on native proteins) is an informative signal, and predictors should be evaluated on the distribution of interest rather than trusted by headline number."),
]));
children.push(body([ b("Keywords: "),
  t("protein solubility; recombinant expression; Escherichia coli; benchmark generalisation; model variance; baselines") ]));

// 1. Introduction
children.push(H1("1. Introduction"));
children.push(body("Recombinant protein production in Escherichia coli is central to molecular biology and biomanufacturing, and insolubility (misfolding into inclusion bodies) is among its commonest failure modes. Many sequence-based predictors triage targets beforehand, from physicochemical scores (ProteinSol, CamSol, the Solubility-Weighted Index) to supervised models (SoluProt, DeepSol, NetSolP) and predictors built on protein language models (RP3Net)."));
children.push(body("These models are usually summarised by a single benchmark AUROC measured on one held-out distribution. Two questions are rarely asked directly: does that headline number predict how the model ranks on a different distribution, and how much do two equally state-of-the-art models differ on the very same proteins? If the spread between top models is large, a single reported number is a poor guide to which model to trust on one's own targets."));
children.push(body([
  t("We therefore benchmarked two independent SOTA deep models (RP3Net and NetSolP) against two dataset-naive composition baselines on two E. coli distributions, native (eSOL) and heterologous (SoluProt), with explicit per-model leakage control and paired significance testing. One baseline is the fallback scorer of an open expression-feasibility tool ("),
  t("TiGer Biotech", { italics: true }),
  t("), used here only as a deliberately simple reference, not as a validated predictor. We find that model choice, not model class, dominates on native proteins."),
]));

// 2. Methods
children.push(H1("2. Methods"));
children.push(H2("2.1 Datasets"));
children.push(body([
  b("eSOL (native E. coli K-12). "),
  t("eSOL (Niwa et al., 2009) reports the solubility of E. coli K-12 ORFs in a chaperone-free reconstituted translation system as a continuous percentage. eSOL identifiers were mapped to sequences via the UniProt E. coli K-12 reference proteome (UP000000625), giving 3,131 proteins. The primary soluble label binarised solubility at \u226570% (60% and 80% swept for robustness). Compartment annotations separated cytoplasmic (in-domain) from membrane and periplasmic (out-of-domain) proteins."),
]));
children.push(body([
  b("SoluProt test set (heterologous). "),
  t("The SoluProt held-out test set (Hon et al., 2021; NESG/TargetTrack-derived) comprises 3,100 balanced heterologous proteins with binary labels; this is the harder distribution, closer to a recombinant-expression triage use case."),
]));
children.push(H2("2.2 Predictors"));
children.push(body([
  b("RP3Net (deep). "),
  t("RP3Net (Tankhilevich et al., 2026) is an ESM-2 650M (Lin et al., 2023) model with set-transformer pooling, LoRA fine-tuning and meta-label correction, trained on AstraZeneca internal data plus the public Structural Genomics Consortium (SGC). The authors report AUROC 0.83 on an independent prospective set of 97 AstraZeneca constructs. We used the released production model and took the predicted soluble probability."),
]));
children.push(body([
  b("NetSolP (deep). "),
  t("NetSolP (Thumuluri et al., 2022) is an ESM1b-based ensemble trained on the PSI:Biology solubility dataset. We ran the released distilled model (NetSolP-D) and took its predicted solubility. Sequences are truncated at 1,022 residues by the released code (the ESM context limit); the few longer proteins are therefore scored on their N-terminal 1,022 residues."),
]));
children.push(body([
  b("ProteinSol (native-calibrated). "),
  t("ProteinSol (Hebditch et al., 2017) scores 35 sequence features whose weights are calibrated on the native eSOL dataset (Niwa et al., 2009); its scaled-solubility output is therefore an eSOL-derived score. Being fit to eSOL it is circular there, and is evaluated only on the heterologous SoluProt set after removing SoluProt proteins near-identical to eSOL (§2.3). We ran the released batch sequence-prediction code and took the scaled-solubility value."),
]));
children.push(body([
  b("Composition baselines (two), dataset-naive. "),
  t("Two deliberately simple sequence scores serve as baselines. The in-house composition heuristic is a deterministic linear score over amino-acid composition (Wilkinson & Harrison, 1991 charge term; Idicula-Thomas & Balaji, 2005 propensity; aliphatic/aromatic/cysteine terms); its exact form is released (see Code availability). We report its eSOL-naive variant (literature-derived term signs only; the coefficient magnitudes are round, author-set values, not tuned to any benchmark) as a baseline; a full variant additionally includes MW/pI terms fit to eSOL and is used only in an ablation. Concretely, the naive score is a clipped linear combination (output range 8–95): 50, plus 55×(fraction D+E+K+R), minus 45×(fraction I+L+V+M), minus 30×(fraction F+W+Y), plus 20×(mean per-residue solubility propensity), minus 7 for more than three cysteines and a further 5 for more than six, plus or minus a fixed offset for two Wilkinson–Harrison-style charge/hydrophobicity threshold predicates (exact predicates in the released code), and minus 15 for sequences under 50 residues. The second baseline is the published Solubility-Weighted Index (SWI; Bhandari et al., 2020), a mean of per-residue flexibility weights optimised on PSI:Biology and validated on eSOL as an external set. Neither baseline is trained per-protein; neither has any exposure to SoluProt. Because its weights are published and externally fixed, SWI is treated as the primary, more conservative baseline; the in-house naive score serves as corroboration and never carries a conclusion on its own."),
]));
children.push(H2("2.3 Per-model train/test leakage control"));
children.push(body([
  t("Each test set was screened against each model's public training sequences with MMseqs2 (Steinegger & Söding, 2017; easy-search, local alignment; leak = best hit \u226590% identity over \u226580% query coverage). For eSOL: 0.13% (4/3,131) overlapped RP3Net's SGC training and 4.06% (127/3,129) overlapped the PSI:Biology training shared by NetSolP and SWI; all native-protein results use the doubly-clean intersection (n = 3,000 scored; 2,154 cytoplasmic; two b-numbers mapped many-to-one from eSOL and were de-duplicated). For the SoluProt test set: 1.29% (40/3,100) overlapped RP3Net's training, leaving 3,060 clean proteins; separately, 4.02% (123/3,060) of SoluProt proteins were near-identical to eSOL (mostly 100% identity; the “heterologous” set contains some native E. coli sequences), removed for the ProteinSol evaluation (n = 2,937). Proprietary training portions (RP3Net's AstraZeneca data) cannot be screened, but are drug-discovery targets unlikely to overlap these native/structural-genomics sets. This is a near-duplicate screen (\u226590% identity); it removes near-identical sequences but not homologs, so residual familiarity between a test set and a model\u2019s broader training distribution cannot be excluded. We therefore also ran a homology-level screen: each eSOL query was clustered with the PSI:Biology training sequences (MMseqs2 easy-cluster, min-seq-id 0.30, coverage 0.5), and a query sharing a cluster with any training sequence was treated as a homolog leak. This removed 256/3,000 eSOL proteins (8.5%), versus 4.06% at the \u226590% near-duplicate threshold; all native-protein comparisons were re-run on the homology-clean subset (n = 1,929 cytoplasmic, \u00a73.1)."),
]));
children.push(H2("2.4 Statistical analysis"));
children.push(body([
  t("AUROC via the Mann-Whitney statistic; 95% CIs by case-resampling bootstrap (2,000 resamples). Two predictors on the same proteins were compared by the DeLong test (1988) as the primary significance test and, independently, a paired bootstrap of the AUROC difference as a robustness cross-check (the bootstrap achieved-significance-level cannot resolve p below 1/2,000, so small p-values are reported from DeLong). Significance requires a paired interval excluding zero with concordant DeLong p; \u201cn.s.\u201d is read from an interval spanning zero. Threshold robustness was assessed at 60/70/80%. Across the six pre-specified primary comparisons (the homology-clean re-runs in \u00a73.1 are confirmatory robustness checks of the same hypotheses, not additional tests), the significant native-protein findings survive Holm\u2013Bonferroni correction (NetSolP vs RP3Net and NetSolP vs SWI, adjusted p < 0.005; SWI vs RP3Net, p = 0.006, Bonferroni-adjusted p = 0.036)."),
]));

// 3. Results
children.push(H1("3. Results"));
children.push(H2("3.1 Leakage differs sharply between models"));
children.push(body([
  t("The two deep models had very different overlap with eSOL. Only 0.13% of eSOL matched RP3Net's public training, but 4.06% matched the PSI:Biology training shared by NetSolP and SWI, a thirty-fold difference that, if ignored, would inflate NetSolP and SWI. We therefore evaluate all native-protein comparisons on the doubly-clean intersection (n = 3,000 scored, 2,154 cytoplasmic). Removing its 127 training proteins lowered NetSolP\u2019s cytoplasmic AUROC only from 0.794 to 0.792, so its standing is not a leakage artefact. The SoluProt leakage screen is reported in \u00a73.4.")
]));
children.push(body([ b("Homology-level robustness. "), t("To rule out that NetSolP\u2019s lead reflects residual homology to its PSI:Biology training rather than near-duplicates, we re-ran the native comparisons after homology-level deduplication (MMseqs2 clustering of eSOL with the training set at \u226430% identity; 256/3,000 = 8.5% removed). On the homology-clean cytoplasmic set (n = 1,929) the ordering was unchanged and, if anything, sharper: NetSolP 0.807 (95% CI 0.786\u20130.828), SWI 0.751 (0.729\u20130.774), RP3Net 0.714 (0.691\u20130.736). NetSolP still beat SWI (\u0394 +0.055; DeLong p = 3.8\u00d710\u207b\u2075) and RP3Net (\u0394 +0.093; p = 3.5\u00d710\u207b\u00b9\u00b9), and SWI still beat RP3Net (\u0394 +0.038; p = 0.009); the naive heuristic remained tied with RP3Net (\u0394 +0.013, p = 0.40). The differences were marginally larger than on the near-duplicate-clean set, so the observed ranking is not an artefact of residual training homology.") ], { spacing: { before: 80, after: 80 } }));
children.push(H2("3.2 On native proteins, two deep models straddle the simple interpretable baselines"));
children.push(body([
  t("On native cytoplasmic eSOL proteins (n = 2,154, doubly leak-free), the two SOTA deep models occupied the two extremes (Fig. 1, Tables 1\u20132). NetSolP was the single best predictor (AUROC 0.792, 95% CI 0.771\u20130.812); RP3Net ranked lowest (0.709, 0.686\u20130.731). Their difference, 0.083 (paired \u0394, 95% CI +0.057 to +0.108; DeLong z = 6.12, p < 1e-8; bootstrap below its 1/2,000 resolution), is larger than the gap between any model and any baseline; model choice within \u201cSOTA deep\u201d dominates here."),
]));
children.push(body([
  t("The two simple interpretable baselines fell between the deep models. SWI (0.745) significantly beat RP3Net (\u0394 +0.037, 95% CI +0.010 to +0.063; p = 0.006; DeLong p = 0.008) and the naive heuristic (0.725) tied it; yet both baselines were significantly beaten by NetSolP (vs SWI \u0394 +0.046, p < 0.001; DeLong p = 0.0004). RP3Net thus ranked below a simple published score on native proteins, despite a reported 0.83 on its own prospective set, while a second, equally state-of-the-art deep model topped everything. The headline benchmark did not predict the native-distribution rank."),
]));
children.push(body([
  b("Heuristic ablation. "),
  t("The full in-house heuristic (with MW/pI terms fit to eSOL) reaches 0.787 (95% CI 0.767\u20130.806) on the same leak-free set, but an ablation attributes +0.062 of that to the eSOL-fit terms; the eSOL-naive variant (0.725) is therefore used throughout and the full variant\u2019s number is not claimed."),
]));
children.push(H2("3.3 On heterologous proteins, the deep models agree and lead"));
children.push(body([
  t("On the heterologous SoluProt test set (leak-free, n = 3,060) the picture inverted: the two deep models clustered at the top and were statistically indistinguishable from each other (NetSolP 0.633, 0.612\u20130.652; RP3Net 0.620, 0.599\u20130.640; paired \u0394 +0.014, p = 0.23), both numerically above the simple baselines (SWI 0.598; naive heuristic 0.581). Only NetSolP significantly exceeded SWI here (\u0394 +0.035, p \u2248 0.002); RP3Net significantly beat the naive heuristic (\u0394 +0.038, p = 0.005) but was indistinguishable from SWI (\u0394 +0.021, p = 0.076), and all four AUROCs lie in a narrow 0.58\u20130.63 band where no predictor is strongly discriminative. On this harder distribution, then, the deep models behaved alike and led numerically, the opposite of their wide divergence on native proteins."),
]));
children.push(body([
  t("A third predictor sharpened the distribution effect. ProteinSol, calibrated on native eSOL, dropped to AUROC 0.542 (95% CI 0.521–0.562) on the heterologous SoluProt set (leak-free, n = 2,937), near chance, and significantly below both deep models (vs RP3Net Δ −0.078; vs NetSolP Δ −0.088; both p < 0.0001) and below both simple baselines. A tool calibrated on one distribution was thus nearly uninformative on the other, the mirror image of RP3Net (strong on its own benchmark, weakest on native proteins)."),
]));
children.push(H2("3.4 Threshold robustness and SoluProt leakage"));
children.push(body([
  t("Re-binarising eSOL at 60/70/80% left the cytoplasmic ordering unchanged (Fig. 2). The SoluProt\u2013RP3Net leakage screen removed 40/3,100 (1.29%) training-overlapping proteins; this lowered RP3Net\u2019s SoluProt AUROC from 0.622 to 0.620 and turned its nominally significant edge over SWI on the full set (p = 0.048) into a non-significant one leak-free (p = 0.076), a concrete illustration that even 1\u20132% overlap must be removed before reporting marginal differences."),
]));
children.push(H2("3.5 Out-of-domain compartments"));
children.push(body([
  t("Membrane and periplasmic eSOL proteins are out of domain for cytoplasmic-expression solubility; all predictors behave irregularly there (wide CIs, no consistent ordering) and we draw no conclusions from these subsets (Table 1)."),
]));

// Figure 1
children.push(...img("fig1_auroc.png", 500, 245, [
  b("Figure 1. "),
  t("AUROC (95% CI), all leak-free. Left: on native cytoplasmic eSOL proteins the two deep models (orange) occupy the extremes (NetSolP best 0.792, RP3Net worst 0.709), straddling the two simple interpretable baselines (blue: SWI 0.745, naive heuristic 0.725); SWI significantly beats RP3Net (p = 0.006). Right: on heterologous SoluProt proteins the two deep models lead and agree (n.s.), both above the baselines. ProteinSol (purple, calibrated on native eSOL) collapses to near-chance (0.542) on this heterologous set.", { size: 18 }),
]));

// Table 1
children.push(H2("Table 1. AUROC [95% CI] by dataset and predictor (leak-free)"));
{
  const W = [2150, 650, 1525, 1525, 1525, 1525];
  children.push(new Table({ width: { size: 8900, type: WidthType.DXA }, columnWidths: W,
    rows: [
      tableHeader(["Dataset / subset", "n", "NetSolP (deep)", "RP3Net (deep)", "SWI (simple)", "Heuristic (naive)"], W),
      row(["eSOL cytoplasmic (native)", "2,154", "0.792 [0.771\u20130.812]", "0.709 [0.686\u20130.731]", "0.745 [0.723\u20130.767]", "0.725 [0.702\u20130.746]"], W, { fill: "EAF1F8", bold: true }),
      row(["eSOL all compartments", "3,000", "0.767 [0.749\u20130.785]", "—*", "0.762 [0.744\u20130.780]", "0.731 [0.711\u20130.750]"], W),
      row(["eSOL membrane (OOD)", "562", "0.686 [0.604\u20130.766]", "—", "0.720 [0.645\u20130.789]", "0.671 [0.604\u20130.737]"], W),
      row(["eSOL periplasmic (OOD)", "285", "0.699 [0.631\u20130.766]", "—", "0.734 [0.668\u20130.796]", "0.667 [0.597\u20130.733]"], W),
      row(["SoluProt (heterologous)", "3,060", "0.633 [0.612\u20130.652]", "0.620 [0.599\u20130.640]", "0.598 [0.578\u20130.617]", "0.581 [0.562\u20130.601]"], W, { fill: "EAF1F8", bold: true }),
    ] }));
  children.push(body([ i("Native (eSOL) rows use the doubly leak-free intersection (RP3Net-clean \u2229 PSI:Biology-clean). The full in-house heuristic (MW/pI fit to eSOL; 0.787 [0.767\u20130.806] cytoplasmic on the same set) is omitted from the columns as partly circular; the naive variant is shown. *RP3Net all-compartment / OOD values omitted on the PSI-clean set for brevity; cytoplasmic is the reported comparison. 95% CIs: 2,000 bootstrap resamples. ProteinSol (eSOL-calibrated; Hebditch 2017) is evaluated only on SoluProt (eSOL-dedup, n = 2,937): AUROC 0.542 [0.521\u20130.562], near chance; omitted from the columns as it is circular on eSOL.") ], { spacing: { before: 80, after: 200 } }));
}

// Table 2
children.push(H2("Table 2. Paired comparisons on shared proteins (leak-free)"));
{
  const W = [3300, 900, 2400, 1300, 1100];
  children.push(new Table({ width: { size: 9000, type: WidthType.DXA }, columnWidths: W,
    rows: [
      tableHeader(["Comparison (paired, same proteins)", "n", "\u0394AUROC [95% CI]", "p (boot/DeLong)", "Verdict"], W),
      row(["eSOL: NetSolP vs RP3Net", "2,154", "+0.083 [+0.057, +0.108]", "<0.001 (DeLong <1e-8)", "NetSolP \u226b RP3Net"], W, { fill: "EAF1F8", bold: true }),
      row(["eSOL: NetSolP vs SWI", "2,154", "+0.046 [+0.022, +0.070]", "<0.001 / 0.0004", "NetSolP wins"], W),
      row(["eSOL: SWI vs RP3Net", "2,154", "+0.037 [+0.010, +0.063]", "0.006 / 0.008", "SWI > RP3Net"], W),
      row(["eSOL: naive heuristic vs RP3Net", "2,154", "+0.016 [\u22120.012, +0.045]", "0.26 / 0.27", "indistinguishable"], W),
      row(["SoluProt: NetSolP vs RP3Net", "3,060", "+0.014 [\u22120.009, +0.036]", "0.23 / 0.24", "indistinguishable"], W, { fill: "FBEEE6" }),
      row(["SoluProt: RP3Net vs naive heuristic", "3,060", "+0.038 [+0.014, +0.062]", "0.005 / 0.004", "RP3Net wins"], W, { fill: "FBEEE6" }),
      row(["SoluProt: RP3Net vs ProteinSol", "2,937", "+0.078 [+0.051, +0.103]", "<0.001 (DeLong <1e-8)", "RP3Net \u226b ProteinSol"], W, { fill: "FBEEE6" }),
      row(["SoluProt: NetSolP vs ProteinSol", "2,937", "+0.088 [+0.068, +0.109]", "<0.001 (DeLong <1e-8)", "NetSolP \u226b ProteinSol"], W, { fill: "FBEEE6" }),
    ] }));
  children.push(body([ i("Native (eSOL) comparisons on the doubly leak-free intersection (n = 2,154 cytoplasmic). The two deep models span 0.083 AUROC on native proteins yet are indistinguishable on heterologous proteins.") ], { spacing: { before: 80, after: 200 } }));
}

// Figure 2
children.push(...img("fig2_threshold.png", 360, 238, [
  b("Figure 2. "),
  t("Threshold robustness on cytoplasmic eSOL (leak-free, n = 2,154): the predictor ordering (NetSolP > SWI > naive heuristic > RP3Net) is stable across 60/70/80% solubility cutoffs.", { size: 18 }),
]));

// 4. Discussion
children.push(H1("4. Discussion"));
children.push(body([
  t("Across two independent E. coli benchmarks, solubility-predictor performance was governed by model and distribution rather than by model class. On native proteins the two state-of-the-art deep models did not cluster: they occupied opposite extremes, 0.083 AUROC apart, a larger gap than that between any model and a simple interpretable baseline. One deep model (NetSolP) topped every predictor; the other (RP3Net) fell below a simple published score. On the harder heterologous distribution the same two models instead agreed and led. The class \u201cdeep learning\u201d was therefore a weak guide to performance here; the specific model and the target distribution mattered far more."),
]));
children.push(body([
  b("Why the two distributions differ. "),
  t("Native E. coli cytoplasmic proteins are co-adapted to the host’s folding, chaperone and codon environment, and eSOL measures them without expression tags. The SoluProt constructs are heterologous, carry purification tags, and span signal-peptide-bearing, membrane-associated and intrinsically disordered proteins whose solubility depends on features (long hydrophobic spans, disorder, host-mismatched codon usage) that differ from native cytoplasmic determinants. A predictor that keys on one regime’s features need not transfer to the other, a plausible route to the model- and distribution-dependence seen here."),
]));
children.push(body([
  b("Benchmarks do not predict generalisation. "),
  t("RP3Net reports 0.83 on its prospective industrial set yet was ranked lowest of four predictors on native E. coli proteins and indistinguishable from the best only on heterologous ones. A single headline AUROC, measured on one distribution, did not rank the models on another. This is not a criticism of any one model (each may be well-suited to its own target distribution) but a caution against selecting a solubility predictor by its reported number."),
]));
children.push(body([
  b("The distribution effect is symmetric. "),
  t("The reverse case makes the same point from the other side: ProteinSol, whose weights are calibrated on native eSOL, fell to near-chance (0.542) on heterologous proteins, just as RP3Net, strong on its prospective benchmark, fell below a simple published score on native ones. Each tool was good on a distribution like its own and weak on the other, and in neither case did the reported benchmark warn of the gap. Distribution match, not model class or headline number, determined performance."),
]));
children.push(body([
  b("A simple interpretable baseline as a sanity floor. "),
  t("The two simple composition scores were not the headline result, but they served a practical purpose: they bracket the deep models and provide a floor. A SOTA model that falls below a dataset-naive baseline on the target distribution (as RP3Net does on native proteins) is signalling a generalisation gap that its benchmark number conceals. Reporting such a baseline alongside any solubility predictor is cheap insurance. Native-proteome benchmarks such as eSOL are also an easy distribution where simple scores are competitive, and should not on their own be used to claim a predictor\u2019s accuracy."),
]));
children.push(H2("Limitations"));
[
  ["Scope of models. ", "We benchmark two deep models and two simple baselines; the qualitative point (large between-model variance; benchmark does not predict generalisation) may not capture every predictor, and adding further models (DeepSol, CamSol, ccSOL) would sharpen the variance estimate."],
  ["Partial leakage screening. ", "Each test set was screened against the public training portions at the near-duplicate level (\u226590% identity), and the native comparisons were additionally re-run after homology-level deduplication (\u226430% identity clustering; 8.5% of eSOL removed), which left the ranking and every significant difference intact (\u00a73.1). RP3Net's proprietary AstraZeneca data could not be screened (but are unlikely to overlap these native/structural-genomics sets). NetSolP truncates sequences at 1,022 residues, a slight disadvantage on the few longer proteins."],
  ["Label proxies. ", "eSOL solubility is from a chaperone-free reconstituted system and SoluProt labels aggregate heterogeneous sources; both are imperfect proxies for production outcome."],
  ["Non-comparability of headline numbers. ", "RP3Net's 0.83 and NetSolP's reported figures are on their own held-out sets, which we cannot access; our independent AUROCs measure different distributions and do not contradict the authors' results."],
].forEach(([h, txt]) => children.push(new Paragraph({ numbering: { reference: "lims", level: 0 },
  alignment: AlignmentType.JUSTIFIED, spacing: { after: 120, line: 276 },
  children: [ b(h), t(txt) ] })));

// 5. Conclusion
children.push(H1("5. Conclusion"));
children.push(body([
  t("On native E. coli proteins, two state-of-the-art deep solubility models differed by more than the gap to a simple interpretable baseline, with one topping every predictor and the other falling below a simple published score; on a harder heterologous benchmark the two agreed and led. A model's reported benchmark AUROC did not predict its generalisation rank. Solubility predictors should be chosen and reported with the evaluation distribution, and a dataset-naive baseline, explicitly in view."),
]));

// Data & code
children.push(H1("Competing interests and funding"));
children.push(body([
  b("Competing interests. "),
  t("The author develops TiGer Biotech, a sequence-based expression-feasibility tool whose composition heuristic is used here as one of the baselines. This constitutes a competing financial/professional interest. To limit its influence, the in-house heuristic is treated only as a corroborating baseline and never carries a conclusion on its own; the externally published Solubility-Weighted Index (SWI) serves as the primary, more conservative baseline, and the full (eSOL-fit) heuristic variant is explicitly excluded from all claims. RP3Net, NetSolP, SWI and ProteinSol are third-party tools used as released, without involvement of their authors in this work. The author has no other competing interests."),
]));
children.push(body([
  b("Funding. "),
  t("none."),
]));
children.push(H1("Data and code availability"));
children.push(body([
  t("eSOL: Niwa et al. (2009). SoluProt test set: Hon et al. (2021). RP3Net and its public training data: the authors\u2019 release. NetSolP and its PSI:Biology training data: Thumuluri et al. (2022). E. coli K-12 proteome: UniProt UP000000625. Analysis code, including the exact composition heuristic, SWI implementation, per-model leakage dedup, holdout evaluation with bootstrap CIs, and paired comparison (bootstrap + DeLong), with all prediction files at https://github.com/Jakobus-svg/tiger-solubility-benchmark. "),
  i("Bibliographic details verified June 2026; confirm against publisher records before submission."),
]));

// References
children.push(H1("References"));
[
  "Niwa T, Ying B-W, Saito K, Jin W, Takada S, Ueda T, Taguchi H (2009). Bimodal protein solubility distribution revealed by an aggregation analysis of the entire ensemble of Escherichia coli proteins. Proc Natl Acad Sci USA 106:4201\u20134206.",
  "Hon J, Marusiak M, Martinek T, Kunka A, Zendulka J, Bednar D, Damborsky J (2021). SoluProt: prediction of soluble protein expression in Escherichia coli. Bioinformatics 37(1):23\u201328. doi:10.1093/bioinformatics/btaa1102.",
  "Tankhilevich E, Martinez Cuesta S, Barrett I, Berg C, Holmberg Schiavone L, Leach AR (2026). RP3Net: a deep learning model for predicting recombinant protein production in Escherichia coli. Bioinformatics 42(1):btag003. doi:10.1093/bioinformatics/btag003.",
  "Thumuluri V, Mart\u00ednez JJA, Johansen J, Nielsen H, Almagro Armenteros JJ (2022). NetSolP: predicting protein solubility in Escherichia coli using language models. Bioinformatics 38(4):941\u2013946. doi:10.1093/bioinformatics/btab801.",
  "Hebditch M, Carballo-Amador MA, Charonis S, Curtis R, Warwicker J (2017). Protein-Sol: a web tool for predicting protein solubility from sequence. Bioinformatics 33(19):3098\u20133100. doi:10.1093/bioinformatics/btx345.",
  "Bhandari BK, Gardner PP, Lim CS (2020). Solubility-Weighted Index: fast and accurate prediction of protein solubility. Bioinformatics 36(18):4691\u20134698. doi:10.1093/bioinformatics/btaa578.",
  "Wilkinson DL, Harrison RG (1991). Predicting the solubility of recombinant proteins in Escherichia coli. Bio/Technology 9(5):443\u2013448. doi:10.1038/nbt0591-443.",
  "Idicula-Thomas S, Balaji PV (2005). Understanding the relationship between the primary structure of proteins and their amyloidogenic/aggregation propensity. Protein Sci 14:582\u2013592.",
  "Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379:1123\u20131130.",
  "Steinegger M, S\u00f6ding J (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35:1026\u20131028.",
  "DeLong ER, DeLong DM, Clarke-Pearson DL (1988). Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics 44:837\u2013845.",
].forEach(r => children.push(refItem(r)));
// ---------- document ----------
const doc = new Document({
  creator: "Jakob Oeschey",
  title: "Distribution-dependent performance of protein solubility predictors",
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: FONT, color: "1F4E79" },
        paragraph: { spacing: { before: 280, after: 140 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 23, bold: true, font: FONT, color: "2E5496" },
        paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 1 } },
    ],
  },
  numbering: { config: [
    { reference: "refs", levels: [{ level: 0, format: LevelFormat.DECIMAL, text: "%1.",
      alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 460, hanging: 360 } } } }] },
    { reference: "lims", levels: [{ level: 0, format: LevelFormat.DECIMAL, text: "%1.",
      alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 460, hanging: 360 } } } }] },
  ] },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 },
      margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
    footers: { default: new Footer({ children: [ new Paragraph({ alignment: AlignmentType.CENTER,
      children: [ new TextRun({ children: ["Page ", PageNumber.CURRENT], size: 18, color: "888888" }) ] }) ] }) },
    children,
  }],
});

Packer.toBuffer(doc).then(buf => { fs.writeFileSync("TiGer_generalization_preprint.docx", buf);
  console.log("docx written:", buf.length, "bytes"); });
