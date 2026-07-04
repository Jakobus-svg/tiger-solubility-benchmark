#!/usr/bin/env python3
"""
subset_by_allowlist.py — restrict two committed eSOL prediction files to an allowlist of b-numbers
(e.g. the homology-clean set) and report the paired cytoplasmic AUROC comparison, WITHOUT re-running
any model. The predictions are score-invariant to which subset they are read on, so filtering the
existing results/esol_*_clean*.json to data/esol_clean_psi_homology.txt reproduces the homology-clean
numbers exactly.

    python scripts/subset_by_allowlist.py data/esol_clean_psi_homology.txt \
        results/esol_netsolp_clean_preds.json results/esol_swi_clean_preds.json NetSolP SWI

Prints AUROC + 95% CI for each predictor and the paired difference (bootstrap CI + DeLong p) on the
shared cytoplasmic proteins. Optionally writes the two filtered prediction files with --save-a/--save-b.
"""
import argparse, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _val_stats as st


def load(path, allow, subset="cytoplasmic", threshold=70.0):
    d = json.loads(Path(path).read_text())
    out, kept = {}, []
    for e in d.get("preds", []):
        b = e.get("bnumber") or e.get("id")
        if b not in allow:
            continue
        if subset and subset not in (e.get("compartment") or "").lower():
            continue
        if e.get("solubility_pct") is None or e.get("pred") is None:
            continue
        out[b] = (float(e["pred"]), 1 if float(e["solubility_pct"]) >= threshold else 0)
        kept.append(e)
    return out, kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("allowlist")
    ap.add_argument("file_a")
    ap.add_argument("file_b")
    ap.add_argument("label_a", nargs="?", default="A")
    ap.add_argument("label_b", nargs="?", default="B")
    ap.add_argument("--subset", default="cytoplasmic")
    ap.add_argument("--threshold", type=float, default=70.0)
    ap.add_argument("--save-a", default="", help="write filtered file_a predictions here")
    ap.add_argument("--save-b", default="", help="write filtered file_b predictions here")
    args = ap.parse_args()

    allow = {l.strip() for l in open(args.allowlist) if l.strip()}
    A, keptA = load(args.file_a, allow, args.subset, args.threshold)
    B, keptB = load(args.file_b, allow, args.subset, args.threshold)
    common = sorted(set(A) & set(B))
    if not common:
        raise SystemExit("[subset] FATAL: no shared proteins after filtering.")
    lab = [A[k][1] for k in common]
    sa = [A[k][0] for k in common]
    sb = [B[k][0] for k in common]
    for name, s in ((args.label_a, sa), (args.label_b, sb)):
        a, lo, hi = st.bootstrap_auroc_ci(s, lab)
        print(f"  {name:10} AUROC={a:.4f} [{lo:.4f}-{hi:.4f}]")
    diff, dlo, dhi, dp = st.paired_bootstrap_diff(sa, sb, lab)
    _, _, z, p = st.delong_test(sa, sb, lab)
    print(f"  n({args.subset})={len(common)}  Δ({args.label_a}-{args.label_b})={diff:+.4f} "
          f"[{dlo:+.4f},{dhi:+.4f}]  DeLong z={z:.3f} p={p:.2e}")

    for save, kept, src in ((args.save_a, keptA, args.file_a), (args.save_b, keptB, args.file_b)):
        if save:
            d = json.loads(Path(src).read_text())
            d["preds"] = kept
            d["allowlist"] = args.allowlist
            Path(save).parent.mkdir(parents=True, exist_ok=True)
            Path(save).write_text(json.dumps(d, indent=2))
            print(f"  filtered predictions -> {save}")


if __name__ == "__main__":
    main()
