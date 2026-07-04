#!/usr/bin/env python3
"""
compare_aurocs.py — is the AUROC difference between two solubility predictors significant?

Takes two per-protein prediction files written by --save-preds (one heuristic run, one RP3Net run
of the SAME validation set), aligns them by protein id, and reports on the SHARED proteins:

  * AUROC of each predictor (paired, identical samples + labels)
  * AUROC(A) - AUROC(B) with a 95% paired-bootstrap CI and two-sided p
  * DeLong (1988) z and p as an independent cross-check

This is the test behind the paper's central claim ("a composition heuristic beats RP3Net on native
E. coli proteins"): if the paired CI excludes 0 / p < 0.05, the difference is significant; if it
straddles 0, the claim must be softened.

USAGE
-----
  # eSOL: run twice (heuristic + RP3Net), saving preds, then compare the cytoplasmic subset
  SOL_FORCE_HEURISTIC=1 python scripts/validate_esol_holdout.py --allowlist data/esol_clean_bnumbers.txt --save-preds esol_heur_preds.json
  python scripts/validate_esol_holdout.py --allowlist data/esol_clean_bnumbers.txt --save-preds esol_rp3_preds.json
  python scripts/compare_aurocs.py esol_heur_preds.json esol_rp3_preds.json --threshold 70 --subset cytoplasmic --label-a heuristic --label-b RP3Net

  # SoluProt:
  SOL_FORCE_HEURISTIC=1 python scripts/retrospective_validation.py --save-preds solu_heur_preds.json
  python scripts/retrospective_validation.py --save-preds solu_rp3_preds.json
  python scripts/compare_aurocs.py solu_heur_preds.json solu_rp3_preds.json --label-a heuristic --label-b RP3Net
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _val_stats as st


def load_preds(path, threshold, subset):
    """Return {id: (pred, label)}. Handles both eSOL and SoluProt --save-preds schemas."""
    data = json.loads(Path(path).read_text())
    out = {}
    for e in data.get("preds", []):
        pid = e.get("bnumber") or e.get("id")
        if pid is None:
            continue
        if subset and subset.lower() not in (e.get("compartment") or "").lower():
            continue
        if "label" in e and e["label"] is not None:           # SoluProt: explicit label
            label = int(e["label"])
        elif e.get("solubility_pct") is not None:             # eSOL: binarise on threshold
            label = 1 if float(e["solubility_pct"]) >= threshold else 0
        else:
            continue
        if e.get("pred") is None:
            continue
        out[pid] = (float(e["pred"]), label)
    return out, data.get("engine") or data.get("active")


def main():
    ap = argparse.ArgumentParser(description="Paired AUROC comparison (bootstrap + DeLong).")
    ap.add_argument("file_a")
    ap.add_argument("file_b")
    ap.add_argument("--threshold", type=float, default=70.0,
                    help="solubility%% cutoff for eSOL files (ignored for SoluProt; default 70)")
    ap.add_argument("--subset", default="",
                    help="optional compartment substring filter for eSOL (e.g. cytoplasmic)")
    ap.add_argument("--label-a", default="A", help="display name for file_a's predictor")
    ap.add_argument("--label-b", default="B", help="display name for file_b's predictor")
    ap.add_argument("--bootstrap", type=int, default=2000, help="bootstrap resamples (default 2000)")
    ap.add_argument("--out", default="", help="optional JSON summary path")
    args = ap.parse_args()

    a, eng_a = load_preds(args.file_a, args.threshold, args.subset)
    b, eng_b = load_preds(args.file_b, args.threshold, args.subset)
    shared = sorted(set(a) & set(b))
    if not shared:
        raise SystemExit("[compare] FATAL: no shared protein ids between the two files "
                         "(were they run on the same set? for eSOL use the same --allowlist).")

    # Labels must agree per protein (same ground truth); warn + trust file_a if they ever differ.
    mism = sum(1 for k in shared if a[k][1] != b[k][1])
    if mism:
        print(f"[compare] WARNING: {mism}/{len(shared)} labels differ between files "
              f"(different threshold/subset?) — using file_a's labels.")

    scores_a = [a[k][0] for k in shared]
    scores_b = [b[k][0] for k in shared]
    labels = [a[k][1] for k in shared]
    npos = sum(labels)

    print(f"[compare] {args.label_a} ({eng_a}) vs {args.label_b} ({eng_b})")
    print(f"[compare] shared proteins: {len(shared)}  (positives {npos}, {npos/len(shared)*100:.1f}%)"
          + (f"  subset='{args.subset}'" if args.subset else ""))

    auc_a, lo_a, hi_a = st.bootstrap_auroc_ci(scores_a, labels, n_boot=args.bootstrap)
    auc_b, lo_b, hi_b = st.bootstrap_auroc_ci(scores_b, labels, n_boot=args.bootstrap)
    print(f"  {args.label_a:<12} AUROC={auc_a:.3f}  [95% CI {lo_a:.3f}-{hi_a:.3f}]")
    print(f"  {args.label_b:<12} AUROC={auc_b:.3f}  [95% CI {lo_b:.3f}-{hi_b:.3f}]")

    diff, dlo, dhi, dp = st.paired_bootstrap_diff(scores_a, scores_b, labels, n_boot=args.bootstrap)
    dp_str = f"< {1.0/args.bootstrap:.1g} (below resolution)" if dp == 0 else f"{dp:.4f}"
    print(f"  Δ ({args.label_a} - {args.label_b}) = {diff:+.3f}  "
          f"[95% CI {dlo:+.3f}..{dhi:+.3f}]  paired-bootstrap p={dp_str}  (robustness cross-check)")
    da, db, z, p = st.delong_test(scores_a, scores_b, labels)
    print(f"  DeLong (primary): z={z:.3f}  p={p:.2e}")
    sig = (dlo > 0 or dhi < 0)
    verdict = (f"{args.label_a} > {args.label_b}" if diff > 0 else f"{args.label_b} > {args.label_a}")
    print(f"  → difference is {'SIGNIFICANT' if sig else 'NOT significant'} at 95% "
          f"({verdict if sig else 'CI straddles 0'}).")

    if args.out:
        Path(args.out).write_text(json.dumps({
            "label_a": args.label_a, "engine_a": eng_a, "auroc_a": round(auc_a, 4),
            "auroc_a_ci95": [round(lo_a, 4), round(hi_a, 4)],
            "label_b": args.label_b, "engine_b": eng_b, "auroc_b": round(auc_b, 4),
            "auroc_b_ci95": [round(lo_b, 4), round(hi_b, 4)],
            "n_shared": len(shared), "n_pos": npos, "subset": args.subset or None,
            "diff_a_minus_b": round(diff, 4), "diff_ci95": [round(dlo, 4), round(dhi, 4)],
            "paired_bootstrap_p": round(dp, 4),
            "delong_z": round(z, 4) if z == z else None,
            "delong_p": round(p, 4) if p == p else None,
            "significant_95": bool(sig),
        }, indent=2))
        print(f"[compare] summary → {args.out}")


if __name__ == "__main__":
    main()
