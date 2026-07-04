#!/usr/bin/env python3
"""
dedup_esol_vs_rp3net.py
=======================

Quantify sequence overlap between the eSOL benchmark (native E. coli K-12 ORFs, used as an
independent solubility holdout in validate_esol_holdout.py) and RP3Net's PUBLIC training data,
so the eSOL AUROC can be reported honestly as a holdout rather than a possibly-leaked figure.

WHAT IT MEASURES
----------------
For every eSOL protein (those with a numeric solubility value) it finds the closest match in the
RP3Net training/validation/test sequences by local alignment (MMseqs2), then flags any eSOL
protein whose best hit exceeds an identity + coverage threshold (default >=90% identity over
>=80% of the eSOL sequence) as a potential leak. It writes:

  * a human-readable report to stdout (leak count, identity distribution, the leaked b-numbers),
  * data/esol_clean_bnumbers.txt   -- the NON-leaked eSOL b-numbers (an allowlist that
    validate_esol_holdout.py --allowlist can consume to re-run the holdout leak-free),
  * an optional JSON summary (--out).

IMPORTANT SCOPE / HONESTY CAVEAT
--------------------------------
RP3Net was trained on AstraZeneca (AZ) INTERNAL data + the Structural Genomics Consortium (SGC,
Stockholm + Toronto). Only the SGC portion is in the public release used here
(https://ftp.ebi.ac.uk/pub/software/RP3Net/v0.1/data/); the AZ portion is proprietary and is NOT
available, so this dedup is against the SGC portion only. That is the realistic leak source for a
native E. coli K-12 benchmark: AZ targets are drug-discovery proteins (diverse, frequently
human/eukaryotic) and are very unlikely to coincide with the E. coli K-12 cytoplasmic proteome.
State this scope whenever citing the deduped result.

REQUIREMENTS
------------
  * MMseqs2 binary on PATH (or pass --mmseqs /path/to/mmseqs). Get a static build from
    https://mmseqs.com (e.g. mmseqs-linux-avx2.tar.gz) -- no compilation needed.
  * The eSOL table (data/esol.csv, bundled) and the E. coli K-12 proteome lookup
    (data/ecoli_k12_proteome.tsv; auto-downloaded from UniProt if absent, same as
    validate_esol_holdout.py -- nothing is trained on it, it only maps b-number -> sequence).

USAGE
-----
  python scripts/dedup_esol_vs_rp3net.py
  python scripts/dedup_esol_vs_rp3net.py --id 90 --cov 0.8 --out dedup_summary.json
  python scripts/dedup_esol_vs_rp3net.py --mmseqs /tmp/mmseqs/bin/mmseqs

Then re-run the holdout on the clean set:
  python scripts/validate_esol_holdout.py --allowlist data/esol_clean_bnumbers.txt
"""
import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

# ── data sources ────────────────────────────────────────────────────────────
RP3NET_FASTA_URL = "https://ftp.ebi.ac.uk/pub/software/RP3Net/v0.1/data/rp3.fasta.gz"
PROTEOME_URL = ("https://rest.uniprot.org/uniprotkb/stream?"
                + urllib.parse.urlencode({
                    "query": "proteome:UP000000625",
                    "fields": "accession,gene_oln,sequence",
                    "format": "tsv"}))


# ── proteome map (b-number -> sequence); identical contract to validate_esol_holdout ─────────
def _proteome_looks_valid(path: Path) -> bool:
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            if not fh.readline().lower().startswith("entry\t"):
                return False
            return sum(1 for _ in fh) > 1000
    except OSError:
        return False


def load_proteome_map(cache: Path) -> dict:
    """b-number (ordered locus name, e.g. b0003) -> protein sequence (UniProt lookup, not trained on)."""
    if cache.exists() and not _proteome_looks_valid(cache):
        print(f"[dedup] cached {cache.name} looks invalid — re-downloading", flush=True)
        cache.unlink()
    if not cache.exists():
        print(f"[dedup] downloading E. coli K-12 proteome → {cache.name} ...", flush=True)
        urllib.request.urlretrieve(PROTEOME_URL, cache)
        if not _proteome_looks_valid(cache):
            raise SystemExit(f"[dedup] FATAL: proteome download to {cache} is invalid. "
                             f"Download manually to that path:\n  {PROTEOME_URL}")
    bmap = {}
    with open(cache, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for row in rd:
            if len(row) < 3:
                continue
            _, oln, seq = row[0], row[1], (row[2] or "").strip().upper()
            if not seq:
                continue
            for tok in oln.split():               # "b0003 JW..." → take the b-number token
                tok = tok.strip()
                if tok.startswith("b") and tok[1:].split(".")[0].isdigit():
                    bmap.setdefault(tok.split(".")[0], seq)
    if len(bmap) < 1000:
        raise SystemExit(f"[dedup] FATAL: only {len(bmap)} b-number→sequence entries parsed "
                         f"from {cache} (expected ~4400). Delete it and re-run to re-download.")
    return bmap


def load_esol_queries(esol_csv: Path, bmap: dict) -> dict:
    """{b-number: sequence} for eSOL rows that carry a numeric solubility (the scored set)."""
    out = {}
    rows = list(csv.DictReader(open(esol_csv, encoding="latin-1")))
    for r in rows:
        b = (r.get("B number") or "").strip()
        sol = (r.get("Solubility (%)") or "").strip()
        seq = bmap.get(b)
        if not seq or not sol or b in out:
            continue
        try:
            float(sol)
        except ValueError:
            continue
        out[b] = seq
    return out


def ensure_rp3net_fasta(cache: Path) -> Path:
    """Download + decompress the public RP3Net training FASTA (SGC Stockholm+Toronto) if absent."""
    if cache.exists() and cache.stat().st_size > 1_000_000:
        return cache
    gz = cache.with_suffix(cache.suffix + ".gz")
    print(f"[dedup] downloading RP3Net training sequences → {gz.name} ...", flush=True)
    urllib.request.urlretrieve(RP3NET_FASTA_URL, gz)
    with gzip.open(gz, "rb") as fi, open(cache, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    n = sum(1 for line in open(cache) if line.startswith(">"))
    if n < 10000:
        raise SystemExit(f"[dedup] FATAL: RP3Net FASTA looks too small ({n} seqs). "
                         f"Check the download:\n  {RP3NET_FASTA_URL}")
    print(f"[dedup] RP3Net training sequences: {n}", flush=True)
    return cache


def write_fasta(path: Path, seqmap: dict) -> None:
    with open(path, "w") as fh:
        for k, seq in seqmap.items():
            fh.write(f">{k}\n{seq}\n")


def run_mmseqs(mmseqs: str, query: Path, target: Path, out_m8: Path, tmp: Path,
               sensitivity: float, max_seqs: int, threads: int) -> None:
    if shutil.which(mmseqs) is None and not Path(mmseqs).exists():
        raise SystemExit(
            f"[dedup] FATAL: MMseqs2 not found ('{mmseqs}'). Install a static build from "
            f"https://mmseqs.com and put it on PATH, or pass --mmseqs /path/to/mmseqs.")
    cmd = [mmseqs, "easy-search", str(query), str(target), str(out_m8), str(tmp),
           "--min-seq-id", "0.0", "-c", "0.5", "--cov-mode", "0",
           "-s", str(sensitivity), "--max-seqs", str(max_seqs), "--threads", str(threads),
           "--format-output", "query,target,pident,qcov,tcov,evalue,bits"]
    print(f"[dedup] running MMseqs2 easy-search (s={sensitivity}, threads={threads}) ...", flush=True)
    res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stderr[-2000:])
        raise SystemExit(f"[dedup] FATAL: MMseqs2 exited {res.returncode}.")


def best_hits(m8: Path) -> dict:
    """query b-number -> (pident%, qcov, tcov, target) for the single best (max-identity) hit."""
    best = {}
    with open(m8) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            q, t, pid, qcov, tcov = p[0], p[1], float(p[2]), float(p[3]), float(p[4])
            if q not in best or pid > best[q][0]:
                best[q] = (pid, qcov, tcov, t)
    return best


def main():
    ap = argparse.ArgumentParser(description="Dedup eSOL vs RP3Net (public SGC) training set.")
    ap.add_argument("--esol", default="data/esol.csv", help="eSOL table (default: data/esol.csv)")
    ap.add_argument("--id", type=float, default=90.0,
                    help="identity%% leak threshold (default 90)")
    ap.add_argument("--cov", type=float, default=0.80,
                    help="query-coverage leak threshold 0-1 (default 0.80)")
    ap.add_argument("--mmseqs", default="mmseqs", help="MMseqs2 binary (default: 'mmseqs' on PATH)")
    ap.add_argument("--sensitivity", type=float, default=5.0, help="MMseqs2 -s (default 5.0)")
    ap.add_argument("--max-seqs", type=int, default=50, help="MMseqs2 --max-seqs (default 50)")
    ap.add_argument("--threads", type=int, default=4, help="MMseqs2 threads (default 4)")
    ap.add_argument("--workdir", default="rp3_dedup", help="scratch dir for FASTA/MMseqs output")
    ap.add_argument("--allowlist-out", default="data/esol_clean_bnumbers.txt",
                    help="where to write the clean (non-leaked) b-number allowlist")
    ap.add_argument("--out", default="", help="optional path for a JSON summary")
    args = ap.parse_args()

    esol_path = Path(args.esol)
    work = Path(args.workdir)
    work.mkdir(exist_ok=True)

    bmap = load_proteome_map(esol_path.parent / "ecoli_k12_proteome.tsv")
    print(f"[dedup] proteome map: {len(bmap)} b-number→sequence entries")
    esol = load_esol_queries(esol_path, bmap)
    if not esol:
        raise SystemExit("[dedup] FATAL: 0 eSOL queries built. Check column names / proteome map.")
    print(f"[dedup] eSOL queries (numeric solubility): {len(esol)}")

    query_fa = work / "esol_query.fasta"
    target_fa = ensure_rp3net_fasta(work / "rp3_target.fasta")
    write_fasta(query_fa, esol)

    m8 = work / "esol_vs_rp3.m8"
    run_mmseqs(args.mmseqs, query_fa, target_fa, m8, work / "tmp_mmseqs",
               args.sensitivity, args.max_seqs, args.threads)

    best = best_hits(m8)
    total = len(esol)
    leaks = sorted(
        [(q, pid, qc, t) for q, (pid, qc, tc, t) in best.items() if pid >= args.id and qc >= args.cov],
        key=lambda x: -x[1])
    leak_bs = {q for q, *_ in leaks}
    clean = sorted(b for b in esol if b not in leak_bs)

    # ── report ──────────────────────────────────────────────────────────────
    print(f"\n[dedup] eSOL queries with ANY hit: {len(best)} ({len(best)/total*100:.1f}%)")
    print(f"[dedup] LEAKS at identity>={args.id:.0f}% & qcov>={args.cov*100:.0f}%: "
          f"{len(leaks)} / {total} = {len(leaks)/total*100:.2f}%")
    for q, pid, qc, t in leaks:
        print(f"    {q}  id={pid:5.1f}%  qcov={qc:.2f}  -> {t}")
    print(f"[dedup] clean (non-leaked) eSOL proteins: {len(clean)} / {total}")

    Path(args.allowlist_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.allowlist_out).write_text("\n".join(clean) + "\n")
    print(f"[dedup] clean b-number allowlist → {args.allowlist_out}")
    print("[dedup] re-run the holdout leak-free with:")
    print(f"        python scripts/validate_esol_holdout.py --allowlist {args.allowlist_out}")

    print("\n[dedup] SCOPE: overlap measured vs the PUBLIC RP3Net training data (SGC Stockholm+Toronto)")
    print("        only; the AstraZeneca portion is proprietary and unavailable. AZ targets are")
    print("        drug-discovery proteins, very unlikely to overlap the native E. coli K-12 proteome.")

    if args.out:
        summary = {
            "rp3net_training_source": "public SGC (Stockholm+Toronto); AZ portion proprietary/unavailable",
            "esol_total": total,
            "esol_with_any_hit": len(best),
            "leak_threshold": {"identity_pct": args.id, "query_coverage": args.cov},
            "n_leaks": len(leaks),
            "leak_fraction": round(len(leaks) / total, 4),
            "leaks": [{"bnumber": q, "identity_pct": round(pid, 1),
                       "query_coverage": round(qc, 3), "rp3net_target": t} for q, pid, qc, t in leaks],
            "n_clean": len(clean),
            "allowlist_path": args.allowlist_out,
        }
        Path(args.out).write_text(json.dumps(summary, indent=2))
        print(f"[dedup] JSON summary → {args.out}")


if __name__ == "__main__":
    main()
