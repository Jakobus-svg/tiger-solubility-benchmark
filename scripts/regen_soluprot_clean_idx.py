#!/usr/bin/env python3
"""
regen_soluprot_clean_idx.py — rebuild data/soluprot_clean_idx.txt (the RP3Net-clean row-index
allowlist) from the committed clean prediction file + the raw SoluProt test CSV.

The committed results/solu_*_clean_preds.json store per-protein ids as md5(sequence)[:12] (not row
indices), so the row-index allowlist cannot be recovered from them alone. This helper maps
md5(sequence) -> row index using the raw soluprot_test.csv, so you can regenerate the exact allowlist
that the committed leak-free predictions correspond to — without re-running MMseqs2. It fabricates
nothing: every emitted index is a row whose sequence hashes to an id present in the clean prediction
file.

    python scripts/regen_soluprot_clean_idx.py \
        --clean-preds results/solu_rp3_clean_preds.json \
        --soluprot data/soluprot_test.csv \
        --out data/soluprot_clean_idx.txt

(If you have MMseqs2 and RP3Net's training FASTA, dedup_soluprot_vs_rp3net.py regenerates the same
list from scratch; this helper is the offline path.)
"""
import argparse, csv, hashlib, json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean-preds", default="results/solu_rp3_clean_preds.json")
    ap.add_argument("--soluprot", default="data/soluprot_test.csv")
    ap.add_argument("--out", default="data/soluprot_clean_idx.txt")
    args = ap.parse_args()

    clean_ids = {e["id"] for e in json.loads(Path(args.clean_preds).read_text())["preds"]}
    if not clean_ids:
        raise SystemExit(f"[regen] FATAL: no ids in {args.clean_preds}.")

    rows = list(csv.DictReader(open(args.soluprot)))
    first = rows[0]
    seq_col = "seqs" if "seqs" in first else "sequence"
    keep = []
    for i, row in enumerate(rows):
        seq = (row.get(seq_col) or "").upper().strip()
        if not seq or len(seq) < 10:
            continue
        if hashlib.md5(seq.encode()).hexdigest()[:12] in clean_ids:
            keep.append(i)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("\n".join(str(i) for i in keep) + "\n")
    print(f"[regen] {len(keep)} clean row indices (of {len(rows)}) -> {args.out} "
          f"(matched {len(keep)}/{len(clean_ids)} committed clean ids)")
    if len(keep) != len(clean_ids):
        print(f"[regen] WARNING: {len(clean_ids) - len(keep)} committed ids had no matching row — "
              f"check that --soluprot is the same CSV the predictions were generated from.")


if __name__ == "__main__":
    main()
