#!/usr/bin/env python3
"""
dedup_soluprot_vs_rp3net.py
===========================

Quantify sequence overlap between the SoluProt held-out TEST set (NESG/TargetTrack-derived,
n~3,100) and RP3Net's PUBLIC training data (SGC Stockholm+Toronto), so the SoluProt AUROC can be
reported as a genuine holdout for RP3Net rather than a possibly-leaked figure.

WHY THIS MATTERS
----------------
RP3Net was trained on AstraZeneca-internal + SGC data. The SoluProt test set is derived from the
NESG structural-genomics effort. NESG and SGC are different consortia, but both are structural-
genomics target collections, so some overlap is conceivable. This script measures it. (The same
caveat as the eSOL dedup applies: only RP3Net's PUBLIC SGC training portion can be screened; the
AstraZeneca portion is proprietary/unavailable.)

OUTPUT
------
  * stdout report (leak count + identity distribution + the leaked row indices),
  * data/soluprot_clean_idx.txt  -- 0-based indices of the NON-leaked SoluProt rows (an allowlist
    that retrospective_validation.py --allowlist can consume to re-run leak-free),
  * optional JSON summary (--out).

REQUIREMENTS
------------
  * MMseqs2 on PATH (or --mmseqs /path/to/mmseqs).  https://mmseqs.com
  * data/soluprot_test.csv (columns: seqs,labels) -- the same file retrospective_validation.py uses.

USAGE
-----
  python scripts/dedup_soluprot_vs_rp3net.py --mmseqs ./mmseqs/bin/mmseqs --out soluprot_dedup.json
  python scripts/retrospective_validation.py --allowlist data/soluprot_clean_idx.txt
"""
import argparse, csv, gzip, json, shutil, subprocess, sys, urllib.request
from pathlib import Path

RP3NET_FASTA_URL = "https://ftp.ebi.ac.uk/pub/software/RP3Net/v0.1/data/rp3.fasta.gz"


def load_soluprot(path: Path):
    """Return list of (row_index, sequence) for SoluProt test rows."""
    out = []
    with open(path, encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        # tolerate seqs/labels or sequence/label column names
        seq_col = next((c for c in (rd.fieldnames or []) if c.lower() in ("seqs", "sequence", "seq")), None)
        if seq_col is None:
            raise SystemExit(f"[dedup] FATAL: no sequence column in {path} (cols: {rd.fieldnames})")
        for i, r in enumerate(rd):
            s = (r.get(seq_col) or "").strip().upper()
            if s:
                out.append((i, s))
    return out


def ensure_rp3net_fasta(cache: Path) -> Path:
    if cache.exists() and cache.stat().st_size > 1_000_000:
        return cache
    gz = cache.with_suffix(cache.suffix + ".gz")
    print(f"[dedup] downloading RP3Net training sequences -> {gz.name} ...", flush=True)
    urllib.request.urlretrieve(RP3NET_FASTA_URL, gz)
    with gzip.open(gz, "rb") as fi, open(cache, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    return cache


def run_mmseqs(mmseqs, query, target, out_m8, tmp, sens, max_seqs, threads):
    if shutil.which(mmseqs) is None and not Path(mmseqs).exists():
        raise SystemExit(f"[dedup] FATAL: MMseqs2 not found ('{mmseqs}'). See https://mmseqs.com")
    cmd = [mmseqs, "easy-search", str(query), str(target), str(out_m8), str(tmp),
           "--min-seq-id", "0.0", "-c", "0.5", "--cov-mode", "0", "-s", str(sens),
           "--max-seqs", str(max_seqs), "--threads", str(threads),
           "--format-output", "query,target,pident,qcov,tcov,evalue,bits"]
    print(f"[dedup] running MMseqs2 easy-search (s={sens}) ...", flush=True)
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-2000:]); raise SystemExit(f"[dedup] FATAL: MMseqs2 exited {r.returncode}.")


def best_hits(m8: Path, id_thr, cov_thr) -> dict:
    """query -> (leak, pident, qcov, target). leak = True iff ANY hit clears both thresholds
    (>= id_thr% identity AND >= cov_thr query coverage); the stored hit is the highest-identity
    qualifying hit when one exists, else the highest-identity hit (for reporting). Scanning all hits
    (not just the single max-identity one) avoids missing a leak whose qualifying hit has lower
    identity but sufficient coverage."""
    best = {}
    with open(m8) as fh:
        for line in fh:
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
    ap = argparse.ArgumentParser(description="Dedup SoluProt test set vs RP3Net (public SGC) training.")
    ap.add_argument("--soluprot", default="data/soluprot_test.csv")
    ap.add_argument("--id", type=float, default=90.0, help="identity%% leak threshold (default 90)")
    ap.add_argument("--cov", type=float, default=0.80, help="query-coverage threshold 0-1 (default 0.80)")
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument("--sensitivity", type=float, default=5.0)
    ap.add_argument("--max-seqs", type=int, default=50)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--workdir", default="soluprot_dedup")
    ap.add_argument("--allowlist-out", default="data/soluprot_clean_idx.txt")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    work = Path(args.workdir); work.mkdir(exist_ok=True)
    rows = load_soluprot(Path(args.soluprot))
    if not rows:
        raise SystemExit("[dedup] FATAL: 0 SoluProt sequences loaded.")
    print(f"[dedup] SoluProt sequences: {len(rows)}")

    query_fa = work / "soluprot_query.fasta"
    with open(query_fa, "w") as fh:
        for idx, seq in rows:
            fh.write(f">{idx}\n{seq}\n")
    target_fa = ensure_rp3net_fasta(work / "rp3_target.fasta")

    m8 = work / "soluprot_vs_rp3.m8"
    run_mmseqs(args.mmseqs, query_fa, target_fa, m8, work / "tmp_mmseqs",
               args.sensitivity, args.max_seqs, args.threads)

    best = best_hits(m8, args.id, args.cov)
    total = len(rows)
    leaks = sorted([(q, pid, qc, t) for q, (ql, pid, qc, t) in best.items() if ql],
                   key=lambda x: -x[1])
    leak_idx = {q for q, *_ in leaks}
    clean = [str(idx) for idx, _ in rows if str(idx) not in leak_idx]

    print(f"\n[dedup] SoluProt rows with ANY hit: {len(best)} ({len(best)/total*100:.1f}%)")
    print(f"[dedup] LEAKS at identity>={args.id:.0f}% & qcov>={args.cov*100:.0f}%: "
          f"{len(leaks)} / {total} = {len(leaks)/total*100:.2f}%")
    for q, pid, qc, t in leaks[:25]:
        print(f"    row {q}  id={pid:5.1f}%  qcov={qc:.2f}  -> {t}")
    if len(leaks) > 25:
        print(f"    ... and {len(leaks)-25} more")

    Path(args.allowlist_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.allowlist_out).write_text("\n".join(clean) + "\n")
    print(f"[dedup] clean (non-leaked) row indices -> {args.allowlist_out}  ({len(clean)} rows)")
    print(f"[dedup] re-run leak-free: python scripts/retrospective_validation.py --allowlist {args.allowlist_out}")
    print("\n[dedup] SCOPE: vs PUBLIC RP3Net training (SGC) only; AstraZeneca portion proprietary/unavailable.")

    if args.out:
        Path(args.out).write_text(json.dumps({
            "dataset": "SoluProt test set (NESG-derived)",
            "rp3net_training_source": "public SGC (Stockholm+Toronto); AZ portion proprietary/unavailable",
            "n_total": total, "n_with_any_hit": len(best),
            "leak_threshold": {"identity_pct": args.id, "query_coverage": args.cov},
            "n_leaks": len(leaks), "leak_fraction": round(len(leaks)/total, 4),
            "n_clean": len(clean), "allowlist_path": args.allowlist_out,
        }, indent=2))
        print(f"[dedup] JSON summary -> {args.out}")


if __name__ == "__main__":
    main()
