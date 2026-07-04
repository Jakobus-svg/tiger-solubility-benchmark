#!/usr/bin/env python3
"""
validate_esol_holdout.py — end-to-end validation of the TiGer pipeline against eSOL.

eSOL (Niwa et al. 2009, PNAS 106:4201) measures the chaperone-free solubility of
~3,000 native E. coli K-12 ORFs in a reconstituted (PURE) cell-free system. It is a
QUASI-INDEPENDENT holdout for the production solubility model: eSOL is native E. coli
ORFs, whereas RP3Net was trained on SGC structural-genomics constructs and reports its
headline AUROC (0.83) on an independent AstraZeneca prospective set. They are different
data sources; overlap has been MEASURED via scripts/dedup_esol_vs_rp3net.py (~0.13% near-identical
to the public SGC training data). Run that script and pass --allowlist data/esol_clean_bnumbers.txt
here for a leak-free holdout; the AstraZeneca training portion is proprietary/unavailable.

What it measures (the END-TO-END pipeline output, not RP3Net in isolation):
  * AUROC  — predicted solubility score vs eSOL solubility binarised at a threshold
  * Spearman — predicted score vs continuous eSOL solubility (%)
reported overall, and split by eSOL cell location so the in-domain (cytoplasmic) figure
is separated from the out-of-domain (membrane / periplasmic) one — the latter is where
RP3Net is explicitly out of domain and the scope banner should fire.

The eSOL table has NO sequences (only b-numbers), so sequences are mapped from the
UniProt E. coli K-12 reference proteome (UP000000625) by ordered-locus-name (b-number);
the mapping is downloaded once and cached next to the eSOL file.

Whatever solubility path the engine uses is reported via sol_source: with RP3Net live
(sol_source=rp3net*) this is the PRODUCTION number; in a plain environment it falls back
to the leakage-free heuristic and the number is the FALLBACK floor, not production.

Usage:
    ANTHROPIC_KEY=x ADMIN_KEY=x PREMIUM_CODE=x CONTACT_EMAIL=x \
        python scripts/validate_esol_holdout.py /path/to/esol.csv \
            [--threshold 70] [--limit N] [--compartment cytoplasmic|all]
"""
import sys, csv, os, argparse, urllib.parse, urllib.request, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))   # for _val_stats
import _val_stats as st
import swi as _swi
import heuristic as _heur

# main.py (the production tool, needed ONLY for the RP3Net 'tool' engine) is imported lazily so
# this script runs standalone in the paper repo for the heuristic/swi/external engines.
m = None
def _load_main():
    global m
    if m is None:
        import main as _m
        m = _m
    return m


def _load_external_csv(path):
    """Read an external tool's prediction CSV into {id: solubility_score} (higher = more soluble).

    Auto-detects the id column (sid/id/name/header/accession/query/b number) and the solubility
    score column (prefers names containing 'solub'; else score/prob/prediction/usability). Raises
    with a clear message if columns can't be found, so it fails loud rather than silently mismatching.
    """
    import csv as _csv
    rows = list(_csv.DictReader(open(path, encoding="utf-8-sig")))
    if not rows:
        raise SystemExit(f"[external] FATAL: {path} is empty.")
    cols = list(rows[0].keys())
    low = {c.lower().strip(): c for c in cols}
    id_keys = ("sid", "id", "name", "header", "accession", "query", "b number", "bnumber", "fasta_id")
    id_col = next((low[k] for k in id_keys if k in low), None)
    # score: prefer an explicit solubility column, else a generic score/probability column
    score_col = next((c for c in cols if "solub" in c.lower()), None)
    if score_col is None:
        score_col = next((low[k] for k in ("score", "prob", "probability", "prediction", "predicted", "usability") if k in low), None)
    if id_col is None or score_col is None:
        raise SystemExit(f"[external] FATAL: could not auto-detect id/score columns in {path}. "
                         f"Columns present: {cols}. Rename to include an id (e.g. 'sid') and a "
                         f"solubility score (e.g. 'predicted_solubility').")
    out = {}
    for r in rows:
        rid = (r.get(id_col) or "").strip()
        try:
            out[rid] = float(r.get(score_col))
        except (TypeError, ValueError):
            continue
    if not out:
        raise SystemExit(f"[external] FATAL: parsed 0 numeric scores from {path} "
                         f"(id col '{id_col}', score col '{score_col}').")
    print(f"[external] parsed {path}: id col '{id_col}', score col '{score_col}', {len(out)} rows")
    return out

PROTEOME_URL = ("https://rest.uniprot.org/uniprotkb/stream?"
                + urllib.parse.urlencode({
                    "query": "proteome:UP000000625",
                    "fields": "accession,gene_oln,sequence",
                    "format": "tsv"}))


def _proteome_looks_valid(path: Path) -> bool:
    """A good UniProt TSV starts with the 'Entry\\t...' header and has thousands of rows.
    Guards against a download that silently wrote an HTML error page or a truncated file."""
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            if not fh.readline().lower().startswith("entry\t"):
                return False
            return sum(1 for _ in fh) > 1000
    except Exception:
        return False


def load_proteome_map(cache: Path) -> dict:
    """b-number (ordered locus name, e.g. b0003) -> protein sequence.
    NOTE: this is a SEQUENCE LOOKUP table (eSOL ships b-numbers, not sequences) — nothing is
    trained on it. Downloaded once from UniProt and validated; a corrupt cache is re-fetched."""
    if cache.exists() and not _proteome_looks_valid(cache):
        print(f"[esol] cached {cache.name} looks invalid (empty/HTML/truncated) — re-downloading", flush=True)
        cache.unlink()
    if not cache.exists():
        print(f"[esol] downloading E. coli K-12 proteome → {cache.name} ...", flush=True)
        urllib.request.urlretrieve(PROTEOME_URL, cache)
        if not _proteome_looks_valid(cache):
            raise SystemExit(
                f"[esol] FATAL: proteome download to {cache} is invalid "
                f"({cache.stat().st_size if cache.exists() else 0} bytes). Check network/proxy, then "
                f"retry — or download this URL manually to that path:\n  {PROTEOME_URL}")
    bmap = {}
    with open(cache, encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for row in rd:
            if len(row) < 3:
                continue
            _, oln, seq = row[0], row[1], row[2]
            seq = (seq or "").strip().upper()
            if not seq:
                continue
            for tok in oln.split():           # "b0003 JW..." → take the b-number token
                tok = tok.strip()
                if tok.startswith("b") and tok[1:].split(".")[0].isdigit():
                    bmap.setdefault(tok.split(".")[0], seq)
    if len(bmap) < 1000:
        raise SystemExit(
            f"[esol] FATAL: only {len(bmap)} b-number→sequence entries parsed from {cache} "
            f"(expected ~4400). The file is likely truncated or in an unexpected format. "
            f"Delete it and re-run to re-download.")
    return bmap


def auroc(scores, labels):
    """Tie-aware AUROC via the Mann-Whitney U statistic. labels in {0,1}."""
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return float("nan")
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    sum_pos = sum(r for r, y in zip(ranks, labels) if y == 1)
    return (sum_pos - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def spearman(a, b):
    n = len(a)
    if n < 3:
        return float("nan")
    def ranks(x):
        order = sorted(range(n), key=lambda i: x[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg = (i + j) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = ranks(a), ranks(b)
    ma, mb = sum(ra) / n, sum(rb) / n
    cov = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    va = sum((x - ma) ** 2 for x in ra) ** 0.5
    vb = sum((x - mb) ** 2 for x in rb) ** 0.5
    return cov / (va * vb) if va and vb else float("nan")


def report(tag, pred, sol_pct, thr, n_boot=2000):
    if not pred:
        print(f"  {tag:<26} n=   0  (no rows in this subset)")
        return {"tag": tag, "n": 0}
    labels = [1 if s >= thr else 0 for s in sol_pct]
    a, lo, hi = st.bootstrap_auroc_ci(pred, labels, n_boot=n_boot)
    r = spearman(pred, sol_pct)
    npos = sum(labels)
    print(f"  {tag:<26} n={len(pred):>4}  pos={npos:>4} ({npos/len(pred)*100:4.1f}%)  "
          f"AUROC={a:.3f} [95% CI {lo:.3f}-{hi:.3f}]  Spearman={r:.3f}")
    return {"tag": tag, "n": len(pred), "pos": npos,
            "auroc": round(a, 4), "auroc_ci95": [round(lo, 4), round(hi, 4)],
            "spearman": round(r, 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("esol_csv", nargs="?", default="data/esol.csv",
                    help="eSOL table (default: data/esol.csv, bundled)")
    ap.add_argument("--threshold", type=float, default=70.0,
                    help="eSOL solubility%% cutoff for the soluble label (default 70)")
    ap.add_argument("--limit", type=int, default=0, help="cap #proteins (0 = all)")
    ap.add_argument("--allowlist", default="", help="optional file of b-numbers (one per line) to "
                    "restrict the holdout to — e.g. the deduped list from dedup_esol_vs_rp3net.py")
    ap.add_argument("--out", default="", help="optional path to write a JSON summary of the results "
                                              "(by default the script only prints to stdout)")
    ap.add_argument("--bootstrap", type=int, default=2000,
                    help="bootstrap resamples for the AUROC 95%% CI (default 2000; 0 = skip CI speedups)")
    ap.add_argument("--threshold-sweep", default="",
                    help="comma-separated solubility%% cutoffs to also report (e.g. 60,70,80) — "
                         "checks the finding is not an artefact of the 70%% binarisation")
    ap.add_argument("--save-preds", default="",
                    help="optional JSON file of per-protein {bnumber, pred, solubility_pct, compartment} "
                         "— feed two such files (heuristic vs RP3Net) to compare_aurocs.py")
    ap.add_argument("--engine", choices=["tool", "swi", "heuristic", "external"], default="tool",
                    help="'tool' = TiGer solubility path (RP3Net/heuristic via main); "
                         "'swi' = Solubility-Weighted Index (Bhandari 2020); "
                         "'external' = read predictions from --pred-csv (e.g. NetSolP, ProteinSol output)")
    ap.add_argument("--pred-csv", default="",
                    help="for --engine external: CSV from an external tool; auto-detects an id column "
                         "(sid/id/name/header, matched to eSOL b-number) and a solubility score column "
                         "(higher = more soluble)")
    ap.add_argument("--heuristic-naive", action="store_true",
                    help="for --engine heuristic: use the eSOL-naive variant (no MW/pI terms)")
    ap.add_argument("--heuristic-continuous", action="store_true",
                    help="for --engine heuristic: use the un-rounded float score (fewer AUROC ties; "
                         "regenerate all numbers if you enable this — default reproduces the paper)")
    ap.add_argument("--pred-label", default="external predictor",
                    help="name of the external predictor (for reporting / pred-file provenance)")
    ap.add_argument("--export-fasta", default="",
                    help="write the (allowlist-filtered) eSOL sequences as FASTA (header = b-number) and exit; "
                         "this is the input you feed to NetSolP/ProteinSol etc.")
    args = ap.parse_args()

    allow = None
    if args.allowlist:
        allow = {ln.strip() for ln in open(args.allowlist) if ln.strip()}
        print(f"[esol] allowlist: restricting to {len(allow)} b-numbers from {args.allowlist}")

    esol_path = Path(args.esol_csv)
    bmap = load_proteome_map(esol_path.parent / "ecoli_k12_proteome.tsv")
    print(f"[esol] proteome map: {len(bmap)} b-number→sequence entries")

    rows = list(csv.DictReader(open(esol_path, encoding="latin-1")))  # eSOL supplement is latin-1
    recs = []  # (b, compartment, sol_pct, sequence)
    seen_b = {}  # b-number -> index in recs, to drop many-to-one b-number collisions
    n_dup, n_conflict = 0, 0
    for r in rows:
        b = (r.get("B number") or "").strip()
        if allow is not None and b not in allow:
            continue
        sol = (r.get("Solubility (%)") or "").strip()
        seq = bmap.get(b)
        if not seq or not sol:
            continue
        try:
            solf = float(sol)
        except ValueError:
            continue
        comp = (r.get("Cell location") or "").strip().lower()
        if b in seen_b:
            # eSOL->UniProt b-number mapping is not always 1:1 (a few b-numbers appear on multiple
            # eSOL rows). Keep the FIRST occurrence and flag conflicting solubility so it is never
            # silently double-counted with contradictory ground truth (e.g. b2094: 72% vs 14%).
            n_dup += 1
            prev = recs[seen_b[b]]
            if abs(prev[2] - solf) > 1e-9:
                n_conflict += 1
                print(f"[esol] WARNING: b-number {b} maps to multiple eSOL rows with CONFLICTING "
                      f"solubility ({prev[2]}% vs {solf}%); keeping first, dropping the rest.")
            continue
        seen_b[b] = len(recs)
        recs.append((b, comp, solf, seq))
    if n_dup:
        print(f"[esol] de-duplicated {n_dup} many-to-one b-number collision(s) "
              f"({n_conflict} with conflicting solubility).")
    if args.limit:
        recs = recs[:args.limit]
    print(f"[esol] {len(recs)} proteins with sequence + numeric solubility")
    if not recs:
        raise SystemExit(
            "[esol] FATAL: 0 proteins matched. Likely causes: (1) the proteome map is empty/bad "
            "(see messages above), or (2) the eSOL column names differ from what this script expects "
            "('B number', 'Solubility (%)', 'Cell location') — open the CSV header and check.")

    if args.export_fasta:
        with open(args.export_fasta, "w") as fh:
            for b, comp, solf, seq in recs:
                fh.write(f">{b}\n{seq}\n")
        print(f"[esol] wrote {len(recs)} sequences -> {args.export_fasta}")
        print("[esol] feed this FASTA to the external tool, then re-run with "
              "--engine external --pred-csv <tool_output.csv>")
        return

    ext_map = None
    if args.engine == "external":
        if not args.pred_csv:
            raise SystemExit("[esol] FATAL: --engine external requires --pred-csv <tool output CSV>.")
        ext_map = _load_external_csv(args.pred_csv)
        print(f"[esol] external predictor '{args.pred_label}': {len(ext_map)} id→score rows from {args.pred_csv}")

    pred, sol_pct, comp_all, bnum_all, srcs = [], [], [], [], set()
    n_skip = 0
    for k, (bnum, comp, solf, seq) in enumerate(recs):
        if args.engine == "swi":
            sc = _swi.swi_score(seq)
            if sc != sc:        # NaN
                n_skip += 1; continue
            pred.append(float(sc)); sol_pct.append(solf); comp_all.append(comp); bnum_all.append(bnum)
            srcs.add("SWI (Bhandari 2020, eSOL-independent)")
            if (k + 1) % 500 == 0: print(f"  ...{k+1}/{len(recs)}")
            continue
        if args.engine == "heuristic":
            sc = _heur.tiger_heuristic_score(seq, include_mw_pi=not args.heuristic_naive,
                                             continuous=args.heuristic_continuous)
            pred.append(float(sc)); sol_pct.append(solf); comp_all.append(comp); bnum_all.append(bnum)
            srcs.add("TiGer heuristic (naive)" if args.heuristic_naive else "TiGer heuristic (full)")
            continue
        if args.engine == "external":
            sc = ext_map.get(bnum)
            if sc is None:
                n_skip += 1; continue
            pred.append(float(sc)); sol_pct.append(solf); comp_all.append(comp); bnum_all.append(bnum)
            srcs.add(args.pred_label)
            continue
        try:
            a = _load_main().analyze_sequence(seq)   # RP3Net used here if models are live
        except Exception:
            n_skip += 1
            continue
        sc = (a.get("scores") or {}).get("solubility")
        if sc is None:
            n_skip += 1
            continue
        pred.append(float(sc)); sol_pct.append(solf); comp_all.append(comp); bnum_all.append(bnum)
        srcs.add(a.get("sol_source", "?"))
        if (k + 1) % 500 == 0:
            print(f"  ...{k+1}/{len(recs)}", flush=True)

    if not pred:
        raise SystemExit("[esol] FATAL: every sequence failed analysis (0 scored). Check that the "
                         "engine imports and that ANTHROPIC_KEY/ADMIN_KEY/etc. are set.")
    print(f"\n[esol] scored {len(pred)} proteins ({n_skip} skipped)")
    print(f"[esol] solubility path(s) exercised: {sorted(srcs)}")
    print(f"[esol] AUROC binarised at eSOL solubility >= {args.threshold:.0f}%, "
          f"Spearman on continuous solubility(%)\n")

    def sub(pred, sol_pct, comp_all, keep):
        idx = [i for i, c in enumerate(comp_all) if keep(c)]
        return [pred[i] for i in idx], [sol_pct[i] for i in idx]

    results = [report("ALL", pred, sol_pct, args.threshold, args.bootstrap)]
    p, s = sub(pred, sol_pct, comp_all, lambda c: "cytoplasmic" in c)
    if p: results.append(report("cytoplasmic (in-domain)", p, s, args.threshold, args.bootstrap))
    p, s = sub(pred, sol_pct, comp_all, lambda c: "membrane" in c)
    if p: results.append(report("membrane (out-of-domain)", p, s, args.threshold, args.bootstrap))
    p, s = sub(pred, sol_pct, comp_all, lambda c: "periplasm" in c)
    if p: results.append(report("periplasmic (out-of-domain)", p, s, args.threshold, args.bootstrap))

    # Threshold robustness: re-report cytoplasmic + ALL at alternative cutoffs so the finding can be
    # shown to be stable (not an artefact of the 70% binarisation).
    sweep = {}
    if args.threshold_sweep:
        cut = []
        for tok in args.threshold_sweep.split(","):
            tok = tok.strip()
            if tok:
                cut.append(float(tok))
        if cut:
            print(f"\n[esol] threshold sweep ({', '.join(f'{c:.0f}%' for c in cut)}):")
            cyto_p, cyto_s = sub(pred, sol_pct, comp_all, lambda c: "cytoplasmic" in c)
            for thr in cut:
                rr = {
                    "ALL": report(f"ALL @{thr:.0f}%", pred, sol_pct, thr, args.bootstrap),
                    "cytoplasmic": report(f"cytoplasmic @{thr:.0f}%", cyto_p, cyto_s, thr, args.bootstrap),
                }
                sweep[str(thr)] = rr

    if args.save_preds:
        Path(args.save_preds).parent.mkdir(parents=True, exist_ok=True)
        Path(args.save_preds).write_text(json.dumps({
            "engine": sorted(srcs),
            "threshold_pct": args.threshold,
            "preds": [{"bnumber": b, "pred": pv, "solubility_pct": sv, "compartment": cv}
                      for b, pv, sv, cv in zip(bnum_all, pred, sol_pct, comp_all)],
        }, indent=2))
        print(f"[esol] per-protein predictions → {args.save_preds}")

    if args.out:
        summary = {
            "dataset": "eSOL (Niwa 2009, native E. coli K-12 ORFs)",
            "threshold_pct": args.threshold,
            "n_scored": len(pred), "n_skipped": n_skip,
            "sol_source": sorted(srcs),
            "subsets": [r for r in results if r],
            "threshold_sweep": sweep,
        }
        Path(args.out).write_text(json.dumps(summary, indent=2))
        print(f"[esol] summary written -> {args.out}")

    if allow is not None:
        print("\n[dedup / leakage note] Running on the DEDUPLICATED eSOL set (--allowlist): near-identical")
        print("  matches to RP3Net's public SGC training data have been removed (measured ~0.13% leak,")
        print("  >=90% identity / >=80% coverage; see scripts/dedup_esol_vs_rp3net.py). The cytoplasmic")
        print("  AUROC is therefore a leak-free holdout vs the available training data. SCOPE: the")
        print("  AstraZeneca training portion is proprietary/unavailable, but AZ are drug-discovery")
        print("  targets unlikely to overlap the native E. coli K-12 proteome.")
        print("  Membrane/periplasmic rows are OUT OF DOMAIN — expect low/inverted AUROC there; that")
        print("  is the scope-banner regime, not a model failure.")
    else:
        print("\n[dedup / leakage note] eSOL = native E. coli K-12 ORFs (PURE cell-free); RP3Net")
        print("  trained on SGC structural-genomics constructs. Different sources, but this run is NOT")
        print("  deduplicated. Before citing the cytoplasmic AUROC as a holdout, run")
        print("  scripts/dedup_esol_vs_rp3net.py and re-run with --allowlist data/esol_clean_bnumbers.txt")
        print("  (measured leak is ~0.13%, so the deduped number is essentially identical).")
        print("  Membrane/periplasmic rows are OUT OF DOMAIN — expect low/inverted AUROC there; that")
        print("  is the scope-banner regime, not a model failure.")


if __name__ == "__main__":
    main()
