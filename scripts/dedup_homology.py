#!/usr/bin/env python3
"""
dedup_homology.py — HOMOLOGY-level (not just near-duplicate) train/test leakage sensitivity screen.

WHY THIS EXISTS
---------------
The dedup_*_vs_rp3net.py / dedup_vs_training.py screens flag only NEAR-IDENTICAL hits (>=90% identity
over >=80% coverage). That removes duplicates but NOT homologs: two proteins at, say, 40% identity are
distinct sequences yet a model may still transfer between them, so a >=90% screen can leave real
train/test familiarity behind. For a *generalisation* claim, reviewers will (correctly) ask for a
homology-level screen — clustering at <=30% identity, the usual bar for train/test independence.

This script re-runs the eSOL<->training overlap at a configurable identity floor (default 30%) using
MMseqs2 clustering, and re-reports how many test proteins fall in a cluster that also contains a
training sequence. Feed the resulting allowlist back into validate_esol_holdout.py / compare_aurocs.py
to obtain the AUROCs on the homology-clean set. If NetSolP's lead over the baselines survives the 30%
screen, the headline is robust; if it shrinks, it must be reported at the homology-clean level.

    python scripts/dedup_homology.py \
        --query esol_query.fasta \
        --train-fasta <NetSolP>/Datasets/PSI_Biology/pET_full_without_his_tag.fa \
        --id 0.30 --cov 0.5 --mmseqs /path/to/mmseqs \
        --label esol_vs_psi_homology \
        --allowlist-out data/esol_clean_psi_homology.txt \
        --out results/esol_psi_homology_dedup.json

Then:
    python scripts/validate_esol_holdout.py --allowlist data/esol_clean_psi_homology.txt \
        --engine external --pred-csv esol_netsolp.csv --pred-label NetSolP \
        --save-preds results/esol_netsolp_homoclean_preds.json
    python scripts/compare_aurocs.py results/esol_netsolp_homoclean_preds.json \
        results/esol_swi_homoclean_preds.json --subset cytoplasmic --label-a NetSolP --label-b SWI

NOTE: this screen requires the raw training FASTA (not redistributed) and MMseqs2; it produces NO
numbers on its own and fabricates nothing — it is the tool the revision uses to fill in the
homology-clean AUROCs.
"""
import argparse, json, shutil, subprocess, sys
from pathlib import Path


def read_fasta_ids(fa: Path):
    ids = []
    for line in open(fa):
        if line.startswith(">"):
            ids.append(line[1:].strip().split()[0])
    return ids


def write_combined(query_fa: Path, train_fa: Path, out_fa: Path):
    """Concatenate query + train into one FASTA with disjoint, prefixed ids so clusters can be
    attributed back to test vs train. Returns (query_ids, n_train)."""
    q_ids = []
    with open(out_fa, "w") as out:
        for line in open(query_fa):
            if line.startswith(">"):
                qid = line[1:].strip().split()[0]
                q_ids.append(qid)
                out.write(f">Q|{qid}\n")
            else:
                out.write(line)
        n_train = 0
        for line in open(train_fa):
            if line.startswith(">"):
                n_train += 1
                out.write(f">T|{line[1:].strip().split()[0]}\n")
            else:
                out.write(line)
    return q_ids, n_train


def run_cluster(mmseqs, combined_fa: Path, tmp: Path, out_prefix: Path, min_id: float, cov: float,
                cov_mode: int, sens: float, threads: int):
    if shutil.which(mmseqs) is None and not Path(mmseqs).exists():
        raise SystemExit(f"[homology] FATAL: MMseqs2 not found ('{mmseqs}'). See https://mmseqs.com")
    cmd = [mmseqs, "easy-cluster", str(combined_fa), str(out_prefix), str(tmp),
           "--min-seq-id", str(min_id), "-c", str(cov), "--cov-mode", str(cov_mode),
           "-s", str(sens), "--threads", str(threads)]
    print(f"[homology] MMseqs2 easy-cluster (min-seq-id={min_id}, c={cov}, s={sens}) ...", flush=True)
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-2000:])
        raise SystemExit(f"[homology] FATAL: MMseqs2 exited {r.returncode}")
    return Path(f"{out_prefix}_cluster.tsv")


def leaked_queries(cluster_tsv: Path):
    """Read MMseqs cluster TSV (rep<TAB>member). A query id leaks if it shares a cluster with any
    train (T|) member. Returns (set_of_leaked_query_ids, n_clusters)."""
    members = {}  # rep -> list of members
    for line in open(cluster_tsv):
        rep, mem = line.rstrip("\n").split("\t")[:2]
        members.setdefault(rep, []).append(mem)
    leaked = set()
    for rep, mem in members.items():
        has_train = any(m.startswith("T|") for m in mem)
        if not has_train:
            continue
        for m in mem:
            if m.startswith("Q|"):
                leaked.add(m[2:])
    return leaked, len(members)


def main():
    ap = argparse.ArgumentParser(description="Homology-level (default 30% identity) train/test leakage screen.")
    ap.add_argument("--query", required=True, help="query test-set FASTA (from --export-fasta)")
    ap.add_argument("--train-fasta", required=True, help="training set FASTA")
    ap.add_argument("--id", type=float, default=0.30, help="identity floor for clustering (0-1, default 0.30)")
    ap.add_argument("--cov", type=float, default=0.5, help="coverage for clustering (0-1, default 0.5)")
    ap.add_argument("--cov-mode", type=int, default=0)
    ap.add_argument("--sens", type=float, default=7.5, help="MMseqs sensitivity -s (default 7.5; higher = more sensitive)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument("--label", default="homology_screen")
    ap.add_argument("--tmp", default="", help="MMseqs tmp dir (default: <label>_tmp)")
    ap.add_argument("--allowlist-out", default="", help="write the homology-CLEAN query ids (one per line)")
    ap.add_argument("--out", default="", help="write a JSON summary")
    args = ap.parse_args()

    query_fa, train_fa = Path(args.query), Path(args.train_fasta)
    for f in (query_fa, train_fa):
        if not f.exists():
            raise SystemExit(f"[homology] FATAL: {f} not found.")
    workdir = Path(f"{args.label}_work"); workdir.mkdir(exist_ok=True)
    combined = workdir / "combined.fasta"
    q_ids, n_train = write_combined(query_fa, train_fa, combined)
    tmp = Path(args.tmp) if args.tmp else workdir / "tmp"
    cluster_tsv = run_cluster(args.mmseqs, combined, tmp, workdir / args.label,
                              args.id, args.cov, args.cov_mode, args.sens, args.threads)
    leaked, n_clusters = leaked_queries(cluster_tsv)
    clean = [q for q in q_ids if q not in leaked]
    pct = 100.0 * len(leaked) / len(q_ids) if q_ids else 0.0
    print(f"[homology] query={len(q_ids)} train={n_train} clusters={n_clusters}")
    print(f"[homology] LEAKS at identity>={args.id*100:.0f}% (homolog-level): "
          f"{len(leaked)}/{len(q_ids)} ({pct:.2f}%) — vs the ~4% seen at the >=90% near-duplicate screen")
    print(f"[homology] homology-clean queries: {len(clean)}")
    if args.allowlist_out:
        Path(args.allowlist_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.allowlist_out).write_text("\n".join(sorted(clean)) + "\n")
        print(f"[homology] clean allowlist -> {args.allowlist_out}")
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({
            "label": args.label, "identity_floor": args.id, "coverage": args.cov,
            "n_query": len(q_ids), "n_train": n_train, "n_clusters": n_clusters,
            "n_leaked": len(leaked), "pct_leaked": round(pct, 3),
            "n_clean": len(clean), "leaked_ids": sorted(leaked),
        }, indent=2))
        print(f"[homology] summary -> {args.out}")


if __name__ == "__main__":
    main()
