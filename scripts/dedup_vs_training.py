#!/usr/bin/env python3
"""
dedup_vs_training.py
====================

Generic train/test leakage screen: align a QUERY test set (FASTA) against any TRAINING set
(CSV column or FASTA) with MMseqs2 and emit a clean (non-leaked) allowlist of query headers.

Designed for adding more predictors to the generalisation benchmark. SWI, NetSolP and DeepSol are
all trained on the PSI:Biology 'Solubility' dataset (the SoDoPE set), so ONE screen of a test set
against PSI:Biology covers leakage for all three. Build the query FASTA with the validation scripts:

  python scripts/validate_esol_holdout.py --allowlist data/esol_clean_bnumbers.txt --export-fasta esol_query.fasta
  python scripts/retrospective_validation.py --export-fasta soluprot_query.fasta

PSI:Biology training sequences ship with the NetSolP standalone (Downloads tab at
services.healthtech.dtu.dk/services/NetSolP-1.0/ → datasets/PSI_Biology_solubility_trainset.csv)
or the SoDoPE_paper_2020 repo. Point --train-csv at it (auto-detects the sequence column) or pass a
--train-fasta.

USAGE
-----
  python scripts/dedup_vs_training.py --query esol_query.fasta \
      --train-csv PSI_Biology_solubility_trainset.csv --label esol_vs_psibiology \
      --mmseqs ./mmseqs/bin/mmseqs --allowlist-out data/esol_clean_psibiology.txt --out esol_vs_psi.json

The allowlist-out lists the query headers (b-numbers for eSOL, row indices for SoluProt) that are
NOT near-identical to any training sequence, ready for --allowlist on the validation scripts.
"""
import argparse, csv, shutil, subprocess, sys
from pathlib import Path


def csv_to_fasta(csv_path: Path, out_fa: Path) -> int:
    rows = list(csv.DictReader(open(csv_path, encoding="utf-8-sig")))
    if not rows:
        raise SystemExit(f"[dedup] FATAL: {csv_path} empty.")
    cols = list(rows[0].keys()); low = {c.lower().strip(): c for c in cols}
    seq_col = next((low[k] for k in ("sequence", "seq", "seqs", "fasta", "protein_sequence") if k in low), None)
    if seq_col is None:
        raise SystemExit(f"[dedup] FATAL: no sequence column in {csv_path} (cols: {cols}).")
    n = 0
    with open(out_fa, "w") as fh:
        for i, r in enumerate(rows):
            s = (r.get(seq_col) or "").strip().upper()
            if s:
                fh.write(f">train_{i}\n{s}\n"); n += 1
    print(f"[dedup] training: {n} sequences from {csv_path} (col '{seq_col}')")
    return n


def read_query_ids(fa: Path):
    ids = []
    for line in open(fa):
        if line.startswith(">"):
            ids.append(line[1:].strip().split()[0])
    return ids


def run_mmseqs(mmseqs, q, t, m8, tmp, sens, max_seqs, threads):
    if shutil.which(mmseqs) is None and not Path(mmseqs).exists():
        raise SystemExit(f"[dedup] FATAL: MMseqs2 not found ('{mmseqs}'). See https://mmseqs.com")
    cmd = [mmseqs, "easy-search", str(q), str(t), str(m8), str(tmp),
           "--min-seq-id", "0.0", "-c", "0.5", "--cov-mode", "0", "-s", str(sens),
           "--max-seqs", str(max_seqs), "--threads", str(threads),
           "--format-output", "query,target,pident,qcov,tcov,evalue,bits"]
    print(f"[dedup] MMseqs2 easy-search (s={sens}) ...", flush=True)
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-2000:]); raise SystemExit(f"[dedup] FATAL: MMseqs2 exited {r.returncode}")


def best_hits(m8: Path, id_thr, cov_thr):
    """query -> (leak, pident, qcov, target). leak = True iff ANY hit clears both thresholds
    (>= id_thr% identity AND >= cov_thr query coverage); the stored hit is the highest-identity
    qualifying hit when one exists, else the highest-identity hit (for reporting). Scanning all hits
    (not just the single max-identity one) avoids missing a leak whose qualifying hit has lower
    identity but sufficient coverage."""
    best = {}
    for line in open(m8):
        p = line.rstrip("\n").split("\t")
        if len(p) < 7:
            continue
        q, t, pid, qc = p[0], p[1], float(p[2]), float(p[3])
        qualifies = pid >= id_thr and qc >= cov_thr
        cur = best.get(q)
        if cur is None or (qualifies, pid) > (cur[0], cur[1]):
            best[q] = (qualifies, pid, qc, t)
    return best


def main():
    ap = argparse.ArgumentParser(description="Generic test-vs-training leakage screen (MMseqs2).")
    ap.add_argument("--query", required=True, help="query test-set FASTA (from --export-fasta)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--train-csv", help="training set CSV (auto-detects sequence column)")
    g.add_argument("--train-fasta", help="training set FASTA")
    ap.add_argument("--label", default="test_vs_training")
    ap.add_argument("--id", type=float, default=90.0, help="identity%% leak threshold (default 90)")
    ap.add_argument("--cov", type=float, default=0.80, help="query coverage threshold 0-1 (default 0.80)")
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument("--sensitivity", type=float, default=5.0)
    ap.add_argument("--max-seqs", type=int, default=50)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--workdir", default="")
    ap.add_argument("--allowlist-out", required=True, help="where to write the clean (non-leaked) query headers")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    work = Path(args.workdir or f"dedup_{args.label}"); work.mkdir(exist_ok=True)
    query = Path(args.query)
    q_ids = read_query_ids(query)
    if not q_ids:
        raise SystemExit(f"[dedup] FATAL: no headers in {query}.")
    print(f"[dedup] query: {len(q_ids)} sequences from {query}")

    if args.train_fasta:
        target = Path(args.train_fasta)
    else:
        target = work / "train.fasta"
        csv_to_fasta(Path(args.train_csv), target)

    m8 = work / f"{args.label}.m8"
    run_mmseqs(args.mmseqs, query, target, m8, work / "tmp", args.sensitivity, args.max_seqs, args.threads)

    best = best_hits(m8, args.id, args.cov)
    total = len(q_ids)
    leaks = sorted([(q, pid, qc, t) for q, (ql, pid, qc, t) in best.items() if ql],
                   key=lambda x: -x[1])
    leak_ids = {q for q, *_ in leaks}
    clean = [q for q in q_ids if q not in leak_ids]

    print(f"\n[dedup] query rows with ANY hit: {len(best)} ({len(best)/total*100:.1f}%)")
    print(f"[dedup] LEAKS at id>={args.id:.0f}% & qcov>={args.cov*100:.0f}%: {len(leaks)}/{total} = {len(leaks)/total*100:.2f}%")
    for q, pid, qc, t in leaks[:25]:
        print(f"    {q}  id={pid:5.1f}%  qcov={qc:.2f} -> {t}")
    if len(leaks) > 25:
        print(f"    ... and {len(leaks)-25} more")

    Path(args.allowlist_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.allowlist_out).write_text("\n".join(clean) + "\n")
    print(f"[dedup] clean query headers -> {args.allowlist_out}  ({len(clean)} kept)")

    if args.out:
        import json
        Path(args.out).write_text(json.dumps({
            "label": args.label, "n_query": total, "n_with_any_hit": len(best),
            "leak_threshold": {"identity_pct": args.id, "query_coverage": args.cov},
            "n_leaks": len(leaks), "leak_fraction": round(len(leaks)/total, 4),
            "n_clean": len(clean), "allowlist_path": args.allowlist_out,
        }, indent=2))
        print(f"[dedup] JSON -> {args.out}")


if __name__ == "__main__":
    main()
