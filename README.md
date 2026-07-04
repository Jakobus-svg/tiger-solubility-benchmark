# Model- and distribution-dependent performance of protein solubility predictors

Reproducibility code for the preprint *"Model- and distribution-dependent performance of protein
solubility predictors: two state-of-the-art deep models straddle a simple interpretable baseline on native
E. coli proteins."*

This repository contains **only** the independent benchmark and analysis code — not the
solubility-screening tool whose composition heuristic is used here as one baseline.

## Summary of findings

Two independent E. coli benchmarks (**eSOL**, native K-12; **SoluProt** held-out test set,
heterologous) were used to compare **two** state-of-the-art deep models (**RP3Net**, ESM-2 650M;
**NetSolP**, ESM1b) against two dataset-naive composition baselines (an in-house heuristic and the
published Solubility-Weighted Index, SWI). Every test set was deduplicated against **each model's**
public training data at the **near-duplicate level** (MMseqs2, ≥90% identity / ≥80% coverage;
eSOL: 0.13% to RP3Net, 4.06% to NetSolP/SWI training; SoluProt: 1.29% to RP3Net); native-protein
comparisons use the doubly-clean intersection.

> **Near-duplicate vs homology screening.** The ≥90%/≥80% screen removes near-identical sequences,
> not homologs. As a stronger control we also ran a homology-level screen (`scripts/dedup_homology.py`,
> MMseqs2 clustering at ≤30% identity): 256/3,000 eSOL proteins (8.5%) are PSI:Biology homologs, vs
> 4.06% at the ≥90% threshold. On the homology-clean set (n = 1,929 cytoplasmic) the ranking is
> **unchanged and slightly sharper** — NetSolP **0.807**, SWI **0.751**, RP3Net **0.714**; NetSolP > SWI
> (Δ +0.055, p = 4e-5), NetSolP > RP3Net (Δ +0.093, p = 3e-11), SWI > RP3Net (Δ +0.038, p = 0.009) all
> remain significant, so NetSolP's lead is not a residual-training-homology artefact. The clean
> allowlist (`data/esol_clean_psi_homology.txt`) and screen output (`results/esol_psi_homology_dedup.json`)
> are committed; reproduce with `scripts/subset_by_allowlist.py` (no model re-run needed).

**eSOL cytoplasmic, leak-free (n=2,154):** the two deep models occupy the extremes and straddle the
simple interpretable baselines.

| predictor | type | AUROC |
|---|---|---|
| NetSolP | deep | **0.792** (best) |
| SWI | simple | 0.745 |
| naive heuristic | simple | 0.725 |
| RP3Net | deep | **0.709** (worst) |

| paired comparison (eSOL cyto, leak-free) | ΔAUROC | p (boot / DeLong) |
|---|---|---|
| NetSolP vs RP3Net | +0.083 | <0.0001 |
| NetSolP vs SWI | +0.046 | <0.001 / 0.0004 |
| SWI vs RP3Net | +0.037 | 0.006 / 0.008 |
| naive heuristic vs RP3Net | +0.016 | 0.26 / 0.27 (n.s.) |

The two SOTA deep models differ by **0.083 AUROC** on native proteins — more than the gap between
any model and a baseline. RP3Net (despite its reported 0.83) is beaten by the trivial SWI score;
NetSolP tops everything. On **SoluProt** (leak-free, n=3,060) the deep models instead agree and lead
numerically (NetSolP 0.633, RP3Net 0.620, n.s.). There, only NetSolP significantly beats the strong
baseline SWI (0.598; Δ +0.035, p≈0.002); RP3Net beats the weak naive baseline (0.581) but is
**statistically indistinguishable from SWI** (Δ +0.021, p = 0.076), so the "deep beats simple on
heterologous proteins" claim rests mainly on NetSolP, and all predictors sit in a low band
(0.58–0.63). ProteinSol — calibrated on native eSOL — collapses to **0.542** (near-chance) on
SoluProt, the mirror image of RP3Net's native-protein weakness. A model's headline benchmark does
**not** predict its generalisation rank. The full in-house heuristic's eSOL number
(0.787) is partly circular (MW/pI fit to eSOL; ablation drops it to 0.725) and is **not** claimed.

## Layout

```
scripts/    analysis code (stdlib only; matplotlib optional for figures)
  heuristic.py                  standalone TiGer composition score (naive + full variants)
  swi.py                        Solubility-Weighted Index (Bhandari et al. 2020, exact weights)
  _val_stats.py                 AUROC, bootstrap CI, paired bootstrap, DeLong test
  validate_esol_holdout.py      eSOL evaluation (engines: heuristic | swi | external | tool)
  retrospective_validation.py   SoluProt evaluation (same engines)
  compare_aurocs.py             paired AUROC comparison of two prediction files
  dedup_esol_vs_rp3net.py       eSOL vs RP3Net public training (MMseqs2)
  dedup_soluprot_vs_rp3net.py   SoluProt vs RP3Net public training (MMseqs2)
  dedup_vs_training.py          generic test-vs-training leakage screen
data/       committed allowlists (deduplicated id lists)
results/    prediction files + comparison summaries (incl. RP3Net predictions, which require the
            proprietary model to regenerate)
paper/      manuscript (docx), figures, and their generators
```

## Reproduce the headline numbers (no deep model required)

The heuristic and SWI engines need only Python's standard library and the input datasets; the
RP3Net predictions are provided in `results/` (they require the production model to regenerate).

> **Two clean sets, by design.** The commands below evaluate on the **RP3Net-clean** eSOL set
> (only the 4 RP3Net near-duplicates removed): naive **0.723**, full 0.787, SWI **0.743**. The paper's
> *native-protein comparisons* instead use the **doubly-clean intersection** (RP3Net-clean ∩
> PSI:Biology-clean, n = 3,000 scored; 2,154 cytoplasmic), so that NetSolP/SWI are on equal footing;
> on that set the same scores are naive **0.725** / SWI **0.745**. The ~0.002 shift is the PSI:Biology
> proteins, not a discrepancy. The doubly-clean allowlist is now committed
> (`data/esol_clean_both.txt`, 3,000 b-numbers, derived from the leak-free prediction files); use it
> to reproduce the exact table values directly, no NetSolP-dedup re-run required. (eSOL→UniProt is not
> strictly 1:1: two b-numbers map to multiple eSOL rows and are de-duplicated on load — one, b2094,
> with conflicting solubility across its source rows; the loader now warns and keeps the first.)

```bash
# eSOL — exact TABLE values (doubly-clean set, committed allowlist): naive 0.725 / SWI 0.745
python scripts/validate_esol_holdout.py --allowlist data/esol_clean_both.txt --engine heuristic --heuristic-naive --save-preds results/esol_naive_clean_preds.json
python scripts/validate_esol_holdout.py --allowlist data/esol_clean_both.txt --engine swi                        --save-preds results/esol_swi_clean_preds.json
# (RP3Net-clean set — only the 4 RP3Net near-duplicates removed: naive 0.723 / SWI 0.743)
python scripts/validate_esol_holdout.py --allowlist data/esol_clean_bnumbers.txt --engine heuristic --heuristic-naive --save-preds results/esol_naive_preds.json
python scripts/validate_esol_holdout.py --allowlist data/esol_clean_bnumbers.txt --engine swi                        --save-preds results/esol_swi_preds.json

# paired comparisons vs the provided RP3Net predictions (native table = doubly-clean)
python scripts/compare_aurocs.py results/esol_swi_clean_preds.json   results/esol_rp3_clean2_preds.json --subset cytoplasmic --label-a SWI            --label-b RP3Net
python scripts/compare_aurocs.py results/esol_naive_clean_preds.json results/esol_rp3_clean2_preds.json --subset cytoplasmic --label-a heuristic-naive --label-b RP3Net

# SoluProt — the leak-free comparisons run directly on the committed clean prediction files:
python scripts/compare_aurocs.py results/solu_swi_clean_preds.json   results/solu_rp3_clean_preds.json --label-a SWI            --label-b RP3Net
python scripts/compare_aurocs.py results/solu_naive_clean_preds.json results/solu_rp3_clean_preds.json --label-a heuristic-naive --label-b RP3Net
# To regenerate SoluProt predictions from raw data you need the row-index allowlist. Rebuild it
# offline from the committed clean predictions + raw CSV (no MMseqs2 needed):
python scripts/regen_soluprot_clean_idx.py --clean-preds results/solu_rp3_clean_preds.json --soluprot data/soluprot_test.csv --out data/soluprot_clean_idx.txt
python scripts/retrospective_validation.py --allowlist data/soluprot_clean_idx.txt --engine heuristic --heuristic-naive --save-preds results/solu_naive_clean_preds.json
python scripts/retrospective_validation.py --allowlist data/soluprot_clean_idx.txt --engine swi                        --save-preds results/solu_swi_clean_preds.json
```

## NetSolP (second deep model) — run + ingest

NetSolP runs as an external tool; this repo ingests its CSV output (`sid,fasta,predicted_solubility`).

```bash
# 1) export the leak-free test sets as FASTA
python scripts/validate_esol_holdout.py --allowlist data/esol_clean_bnumbers.txt --export-fasta esol_query.fasta
python scripts/retrospective_validation.py --allowlist data/soluprot_clean_idx.txt --export-fasta soluprot_query.fasta
# 2) run the released NetSolP-D on each FASTA -> esol_netsolp.csv, soluprot_netsolp.csv
#    (NetSolP standalone: predict.py --MODEL_TYPE Distilled --PREDICTION_TYPE S)
# 3) dedup eSOL against NetSolP's PSI:Biology training (covers SWI too; 4.06% leak):
python scripts/dedup_vs_training.py --query esol_query.fasta --train-fasta <NetSolP>/Datasets/PSI_Biology/pET_full_without_his_tag.fa --label esol_vs_psi --mmseqs /path/to/mmseqs --allowlist-out data/esol_clean_psi.txt
#    intersect with the RP3Net-clean list for the doubly-clean native set:
python -c "a=set(open('data/esol_clean_bnumbers.txt').read().split()); b=set(open('data/esol_clean_psi.txt').read().split()); open('data/esol_clean_both.txt','w').write('\n'.join(sorted(a&b)))"
# 4) ingest + compare on the doubly-clean set
python scripts/validate_esol_holdout.py --allowlist data/esol_clean_both.txt --engine external --pred-csv esol_netsolp.csv --pred-label NetSolP --save-preds results/esol_netsolp_clean_preds.json
python scripts/compare_aurocs.py results/esol_netsolp_clean_preds.json results/esol_rp3_clean2_preds.json --subset cytoplasmic --label-a NetSolP --label-b RP3Net
```

## Leakage screens (require MMseqs2)

```bash
# near-duplicate screen (>=90% id / >=80% cov) — the numbers reported in the paper
python scripts/dedup_esol_vs_rp3net.py     --mmseqs /path/to/mmseqs --out results/dedup_summary.json
python scripts/dedup_soluprot_vs_rp3net.py --mmseqs /path/to/mmseqs --out results/soluprot_dedup.json

# homology-level sensitivity screen (<=30% id clustering) — the stronger control for generalisation
python scripts/dedup_homology.py --query esol_query.fasta \
    --train-fasta <NetSolP>/Datasets/PSI_Biology/pET_full_without_his_tag.fa \
    --id 0.30 --mmseqs /path/to/mmseqs --label esol_vs_psi_homology \
    --allowlist-out data/esol_clean_psi_homology.txt --out results/esol_psi_homology_dedup.json
# then re-evaluate NetSolP/SWI on data/esol_clean_psi_homology.txt and re-report the lead.
```

### Homology-clean numbers (already committed — reproduce without re-running any model)

```bash
# filter the committed clean predictions to the 2,744 homology-clean b-numbers and compare:
python scripts/subset_by_allowlist.py data/esol_clean_psi_homology.txt results/esol_netsolp_clean_preds.json results/esol_swi_clean_preds.json  NetSolP SWI     # 0.807 vs 0.751, p=4e-5
python scripts/subset_by_allowlist.py data/esol_clean_psi_homology.txt results/esol_netsolp_clean_preds.json results/esol_rp3_clean2_preds.json NetSolP RP3Net  # 0.807 vs 0.714, p=3e-11
python scripts/subset_by_allowlist.py data/esol_clean_psi_homology.txt results/esol_swi_clean_preds.json     results/esol_rp3_clean2_preds.json SWI     RP3Net  # 0.751 vs 0.714, p=0.009
```

## Environment / versions (pin these for exact reproduction)

Record and cite the exact versions used, since MMseqs2 results and model outputs are version-dependent:

- Python 3.9+ (stdlib only for the analysis code; `matplotlib>=3.5` for figures).
- **MMseqs2**: pin the release (`mmseqs version`); the screens use `easy-search`/`easy-cluster`.
- **NetSolP**: released NetSolP-D, `predict.py --MODEL_TYPE Distilled --PREDICTION_TYPE S`.
- **RP3Net**: the authors' released production model (record the release tag/commit).
- **ProteinSol**: released batch sequence-prediction code (record version).
- Bootstrap/DeLong are deterministic given `seed=0`; final reported figures use `--bootstrap 10000`.
  DeLong is the primary significance test; the paired bootstrap is a robustness cross-check and cannot
  resolve p below `1/n_boot` (do not report "p < 0.0001" from a 2,000-resample bootstrap).
- Hardware: model inference (NetSolP/RP3Net/ProteinSol) was run on <record GPU/CPU + runtime>; the
  stdlib analysis (AUROC/bootstrap/DeLong) runs in seconds on a laptop CPU.

## Input datasets (not redistributed here)

Place the following under `data/` (see `data/README.md` for the exact columns expected):

- `esol.csv` — eSOL solubility data, Niwa et al. (2009) PNAS 106:4201.
- `ecoli_k12_proteome.tsv` — UniProt E. coli K-12 reference proteome (UP000000625).
- `soluprot_test.csv` — SoluProt balanced test set, Hon et al. (2021) Bioinformatics 37(1):23.

## Predictors evaluated (third-party, cite the originals)

- **RP3Net** (deep) — Tankhilevich et al. (2026) Bioinformatics 42(1):btag003.
- **NetSolP** (deep) — Thumuluri et al. (2022) Bioinformatics 38(4):941. Run the released NetSolP-D, then ingest its CSV with `--engine external --pred-csv`.
- **SWI** (simple) — Bhandari, Gardner, Lim (2020) Bioinformatics 36(18):4691. Exact weights reproduced in `scripts/swi.py`.
- **Composition heuristic** (simple) — `scripts/heuristic.py` (naive + full variants), reproduces 0.725 / 0.787 exactly.
- **ProteinSol** (native-calibrated) — Hebditch et al. (2017) Bioinformatics 33(19):3098. Calibrated on eSOL → evaluated only on SoluProt (after eSOL-dedup, n=2,937); collapses to 0.542 (near-chance) there.

## Requirements

Python 3.9+ standard library. Optional: `matplotlib` (figures), MMseqs2 (leakage screens),
Node.js + `docx` (manuscript regeneration). See `requirements.txt`.

## License

MIT (see `LICENSE`) for the analysis code in this repository. Third-party predictors, datasets and the
SWI weights retain their original licenses; cite the original works.

**Redistributed predictions.** `results/` contains *derived numeric outputs* (per-protein soluble
scores/probabilities) of third-party models (RP3Net, NetSolP, ProteinSol) on public test sequences,
included only to make the AUROC comparisons reproducible without re-running each model. No model
weights, training data or code are redistributed. If any predictor's license restricts publication of its
outputs, open an issue and the affected `results/*.json` will be replaced by a regeneration script; the
analysis is otherwise reproducible by running each released model yourself (see the run blocks above).
