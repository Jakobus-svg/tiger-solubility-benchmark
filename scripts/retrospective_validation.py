"""
TiGer Biotech – Retrospektive Validierungsstudie

Vergleicht TiGer-Vorhersagen gegen bekannte Outcomes im
SoluProt Test-Set (3.100 Proteine).

Nutzung:  python scripts/retrospective_validation.py
Output:   validation_results.json, retrospective_paper_results.txt, validation_roc.png
          (does NOT touch a hand-curated paper_results.txt)
"""
import sys, os, csv, json, hashlib, argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))   # for _val_stats

# Mock only deps that are genuinely ABSENT — never shadow a real install. Mocking pydantic
# breaks RP3Net's import chain (pydantic.alias_generators), so try the real import first and
# only fall back to a stub when the module truly isn't installed. In the normal repo env
# (fastapi/pydantic present) nothing is stubbed and RP3Net loads as the production model.
import importlib, unittest.mock as _m
for _k in ['fastapi', 'fastapi.responses', 'fastapi.staticfiles', 'pydantic', 'anthropic', 'uvicorn']:
    if _k in sys.modules:
        continue
    try:
        importlib.import_module(_k)
    except Exception:
        sys.modules[_k] = _m.MagicMock()

os.environ.setdefault('ANTHROPIC_KEY', 'x')
os.environ.setdefault('PREMIUM_CODE', '1')

import _val_stats as st
import swi as _swi
import heuristic as _heur

# main.py (production tool) is loaded lazily, ONLY for the RP3Net 'tool' engine, so this script
# runs standalone in the paper repo for the heuristic/swi/external engines.
mod = None
def _load_main():
    global mod
    if mod is None:
        import importlib, unittest.mock as _mk
        for _k in ['fastapi','fastapi.responses','fastapi.staticfiles','pydantic','anthropic','uvicorn']:
            if _k in sys.modules: continue
            try: importlib.import_module(_k)
            except Exception: sys.modules[_k] = _mk.MagicMock()
        os.environ.setdefault('ANTHROPIC_KEY','x'); os.environ.setdefault('PREMIUM_CODE','1')
        import importlib.util
        spec = importlib.util.spec_from_file_location('main','main.py')
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def _load_external_csv(path):
    """Read an external tool's CSV into {id: solubility_score} (higher = more soluble)."""
    import csv as _csv
    rows = list(_csv.DictReader(open(path, encoding="utf-8-sig")))
    if not rows:
        raise SystemExit(f"[external] FATAL: {path} is empty.")
    cols = list(rows[0].keys()); low = {c.lower().strip(): c for c in cols}
    id_keys = ("sid","id","name","header","accession","query","row","fasta_id")
    id_col = next((low[k] for k in id_keys if k in low), None)
    score_col = next((c for c in cols if "solub" in c.lower()), None)
    if score_col is None:
        score_col = next((low[k] for k in ("score","prob","probability","prediction","predicted","usability") if k in low), None)
    if id_col is None or score_col is None:
        raise SystemExit(f"[external] FATAL: could not auto-detect id/score columns in {path}. Columns: {cols}")
    out = {}
    for r in rows:
        rid = (r.get(id_col) or "").strip()
        try: out[rid] = float(r.get(score_col))
        except (TypeError, ValueError): continue
    print(f"[external] parsed {path}: id col '{id_col}', score col '{score_col}', {len(out)} rows")
    return out

def calc_auc(y_true, y_score):
    pairs = sorted(zip(y_score, y_true), reverse=True)
    n_pos = sum(y_true); n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0: return 0.5
    tp = fp = auc = prev_fp = prev_tp = 0; prev_score = None
    for score, label in pairs:
        if score != prev_score:
            auc += (fp - prev_fp) * (tp + prev_tp) / 2
            prev_fp = fp; prev_tp = tp; prev_score = score
        if label == 1: tp += 1
        else: fp += 1
    auc += (fp - prev_fp) * (tp + prev_tp) / 2
    return auc / (n_pos * n_neg)

def run(out_json='validation_results.json', save_preds='', n_boot=2000, allowlist='', engine='tool', pred_csv='', pred_label='external predictor', export_fasta='', heuristic_naive=False, heuristic_continuous=False):
    print("TiGer Biotech – Retrospektive Validierung")
    print("=" * 55)

    data_file = Path('data/soluprot_test.csv')
    with open(data_file) as f:
        rows = list(csv.DictReader(f))

    first    = rows[0]
    seq_col  = 'seqs'   if 'seqs'   in first else 'sequence'
    lbl_col  = 'labels' if 'labels' in first else 'label'

    allow = None
    if allowlist:
        allow = {int(x) for x in Path(allowlist).read_text().split() if x.strip().isdigit()}
        print(f"[retro] allowlist: restricting to {len(allow)} non-leaked row indices from {allowlist}")
    print(f"Test-Set: {len(rows):,} Proteine | Spalten: {seq_col}/{lbl_col}\n")

    if export_fasta:
        n = 0
        with open(export_fasta, "w") as fh:
            for i, row in enumerate(rows):
                if allow is not None and i not in allow: continue
                seq = row.get(seq_col,'').upper().strip()
                if not seq or len(seq) < 10: continue
                fh.write(f">{i}\n{seq}\n"); n += 1
        print(f"[retro] wrote {n} sequences -> {export_fasta} (header = row index)")
        print("[retro] feed to external tool, then re-run --engine external --pred-csv <out.csv>")
        return

    ext_map = None
    if engine == "external":
        if not pred_csv: raise SystemExit("[retro] FATAL: --engine external requires --pred-csv")
        ext_map = _load_external_csv(pred_csv)

    results = []; errors = 0; srcs = set()
    for i, row in enumerate(rows):
        if allow is not None and i not in allow: continue
        seq   = row.get(seq_col,'').upper().strip()
        label = int(row.get(lbl_col, 0))
        if not seq or len(seq) < 10: continue
        if engine == "heuristic":
            sc = _heur.tiger_heuristic_score(seq, include_mw_pi=not heuristic_naive,
                                             continuous=heuristic_continuous)
            srcs.add("TiGer heuristic (naive)" if heuristic_naive else "TiGer heuristic (full)")
            results.append({'id': hashlib.md5(seq.encode()).hexdigest()[:12],
                            'label': label, 'sol_score': sc, 'sol_prob': None})
            continue
        if engine == "external":
            sc = ext_map.get(str(i))
            if sc is None: errors += 1; continue
            srcs.add(pred_label)
            results.append({'id': hashlib.md5(seq.encode()).hexdigest()[:12],
                            'label': label, 'sol_score': sc, 'sol_prob': None})
            continue
        if engine == "swi":
            sc = _swi.swi_score(seq)
            if sc != sc: errors += 1; continue
            srcs.add("SWI (Bhandari 2020, eSOL-independent)")
            results.append({'id': hashlib.md5(seq.encode()).hexdigest()[:12],
                            'label': label, 'sol_score': sc, 'sol_prob': _swi.swi_prob(seq)})
            continue
        try:
            a = _load_main().analyze_sequence(seq)
            srcs.add(a.get('sol_source', ''))
            results.append({
                'id':       hashlib.md5(seq.encode()).hexdigest()[:12],
                'label':    label,
                'sol_score':a['scores']['solubility'],
                'sol_prob': a.get('sol_prob'),
            })
        except: errors += 1
        if (i+1) % 500 == 0: print(f"  {i+1}/{len(rows)}...")

    print(f"\n{len(results):,} analysiert, {errors} Fehler")
    # Which solubility model actually drove scores['solubility'] this run? (Determines labels:
    # with RP3Net live, scores['solubility'] and sol_prob are RP3Net's — NOT 'heuristic'.)
    active = ('external' if engine == 'external'
              else 'swi' if engine == 'swi'
              else 'rp3net' if any('rp3net' in s for s in srcs)
              else 'esm2' if any('esm2' in s for s in srcs)
              else 'heuristic')
    active_label = {'rp3net': 'RP3Net (production)', 'esm2': 'TiGer ESM-2 35M',
                    'heuristic': 'TiGer heuristic', 'swi': 'SWI (Bhandari 2020)',
                    'external': pred_label}[active]
    print(f"Active solubility source: {active}\n")

    labels = [r['label'] for r in results]
    scores = [r['sol_score'] for r in results]
    probs  = [(r['sol_prob'], r['label']) for r in results if r['sol_prob'] is not None]

    auc_h = calc_auc(labels, scores)
    auc_ml = calc_auc([l for p,l in probs], [p for p,l in probs]) if probs else None

    import statistics as _stats
    # NOTE: AUROC is the ONLY metric reported in the paper. accuracy/precision/recall/F1 below use a
    # median score split, which on the balanced SoluProt test set is ~uninformative by construction
    # (it forces ~50/50 predictions); they are kept for diagnostics only and are NOT reported.
    # scores scale: 0-100 for RP3Net/heuristic (threshold 50), 0-1 for SWI -> use median split
    threshold = _stats.median(scores) if engine in ("swi", "external", "heuristic") else 50
    preds     = [1 if s >= threshold else 0 for s in scores]
    accuracy  = sum(p==l for p,l in zip(preds,labels)) / len(labels)
    tp = sum(p==1 and l==1 for p,l in zip(preds,labels))
    fp = sum(p==1 and l==0 for p,l in zip(preds,labels))
    fn = sum(p==0 and l==1 for p,l in zip(preds,labels))
    prec = tp/(tp+fp) if (tp+fp) > 0 else 0
    rec  = tp/(tp+fn) if (tp+fn) > 0 else 0
    f1   = 2*prec*rec/(prec+rec) if (prec+rec) > 0 else 0

    # auc_h = AUC of scores['solubility']; auc_ml = AUC of sol_prob. With RP3Net live BOTH are
    # RP3Net-driven (so they nearly coincide and the headline number IS RP3Net measured on this set).
    measured_auc = auc_ml if (auc_ml is not None and active in ('rp3net', 'esm2')) else auc_h

    # Bootstrap 95% CI on the headline (measured) AUROC, on the SAME quantity reported.
    if measured_auc is auc_ml and auc_ml is not None:
        _bs_scores = [p for p, l in probs]; _bs_labels = [l for p, l in probs]
    else:
        _bs_scores = scores; _bs_labels = labels
    _, ci_lo, ci_hi = st.bootstrap_auroc_ci(_bs_scores, _bs_labels, n_boot=n_boot)

    print(f"MEASURED on this run ({active_label}, n={len(results):,}):")
    print(f"  AUC={measured_auc:.3f} [95% CI {ci_lo:.3f}-{ci_hi:.3f}]  "
          f"Acc={accuracy:.3f}  P={prec:.3f}  R={rec:.3f}  F1={f1:.3f}")

    if save_preds:
        Path(save_preds).parent.mkdir(parents=True, exist_ok=True)
        # Save BOTH scores so compare_aurocs.py can align heuristic-vs-RP3Net runs by sequence id.
        Path(save_preds).write_text(json.dumps({
            "engine": sorted(srcs), "active": active, "threshold": threshold,
            "preds": [{"id": r["id"], "pred": r["sol_score"],
                       "sol_prob": r["sol_prob"], "label": r["label"]} for r in results],
        }, indent=2))
        print(f"  per-protein predictions -> {save_preds}")
    if auc_ml is not None and abs(auc_ml - auc_h) > 0.005:
        print(f"  (scores['solubility'] AUC={auc_h:.3f} | sol_prob AUC={auc_ml:.3f})")

    print("\nReference AUCs (literature / paper figures — NOT measured on this run):")
    refs = [('CamSol (struct.)',          0.73),
            ('ProteinSol',               0.68),
            ('ccSol',                    0.65),
            ('RP3Net paper (AstraZeneca)', 0.83)]
    rows_cmp = [(f'>> {active_label} (this set)', measured_auc)] + refs
    for name, auc in sorted(rows_cmp, key=lambda x: -x[1]):
        bar = '█' * int(auc * 25)
        print(f"  {name:<30} AUC={auc:.3f}  {bar}")
    if active == 'rp3net':
        print("  NOTE: 'RP3Net paper 0.83' is the AstraZeneca prospective figure, a DIFFERENT set;")
        print("        the '>> RP3Net (this set)' row is the one measured on SoluProt just now.")

    out = {'n': len(results), 'active_source': active,
           'auc_measured': round(measured_auc,3),
           'auc_measured_ci95': [round(ci_lo,3), round(ci_hi,3)],
           'auc_scores_solubility': round(auc_h,3),
           'auc_sol_prob': round(auc_ml,3) if auc_ml is not None else None,
           'accuracy': round(accuracy,3), 'precision': round(prec,3),
           'recall': round(rec,3), 'f1': round(f1,3)}
    json.dump(out, open(out_json,'w'), indent=2)

    import datetime
    paper = f"""TiGer measured AUC={measured_auc:.3f} ({active_label}) on the held-out SoluProt test
set (n={len(results):,}), threshold={threshold}/100, accuracy={accuracy*100:.1f}%.

Reference AUCs (literature / paper figures, NOT measured on this run):
  CamSol structure-based:            AUC=0.73
  ProteinSol (Hebditch et al. 2017): AUC=0.68
  ccSol (Agostini et al. 2012):      AUC=0.65
  RP3Net paper (AstraZeneca):        AUC=0.83  (different, prospective set)

Active solubility source this run: {active}. SoluProt test is independent of RP3Net's SGC
training source, so with active=rp3net this AUC is a genuine independent-holdout figure;
before external use, identity-dedup SoluProt against RP3Net's training sequences.

Evaluated: {datetime.date.today()} | Dataset: SoluProt (GleghornLab)
"""
    # NOTE: written to a dedicated file so it never clobbers a hand-curated paper_results.txt
    open('retrospective_paper_results.txt','w').write(paper)
    print(f"\nGespeichert: {out_json}, retrospective_paper_results.txt")

    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        n_pos = sum(labels); n_neg = len(labels)-n_pos
        thresholds = list(range(0,101,2))
        tprs = [sum(s>=t and l==1 for s,l in zip(scores,labels))/n_pos for t in thresholds]
        fprs = [sum(s>=t and l==0 for s,l in zip(scores,labels))/n_neg for t in thresholds]
        plt.figure(figsize=(7,6))
        plt.plot(fprs,tprs,'teal',lw=2,label=f'{active_label} (AUC={measured_auc:.3f})')
        plt.plot([0,1],[0,1],'k--',alpha=0.4)
        plt.xlabel('FPR'); plt.ylabel('TPR')
        plt.title('TiGer Biotech – ROC Curve\nSoluProt Test Set')
        plt.legend(); plt.grid(alpha=0.3)
        plt.savefig('validation_roc.png',dpi=150,bbox_inches='tight')
        print("ROC-Kurve: validation_roc.png")
    except: pass

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="SoluProt held-out validation (RP3Net live or, with "
                                             "SOL_FORCE_HEURISTIC=1, the leakage-free heuristic).")
    ap.add_argument("--out", default="validation_results.json", help="JSON summary path")
    ap.add_argument("--save-preds", default="", help="optional per-protein predictions JSON "
                    "(feed heuristic + RP3Net runs to compare_aurocs.py)")
    ap.add_argument("--bootstrap", type=int, default=2000, help="bootstrap resamples for the 95%% CI")
    ap.add_argument("--allowlist", default="", help="row-index allowlist from dedup_soluprot_vs_rp3net.py (leak-free re-run)")
    ap.add_argument("--engine", choices=["tool","swi","heuristic","external"], default="tool", help="'swi'=SWI; 'heuristic'=TiGer composition; 'external'=read --pred-csv")
    ap.add_argument("--heuristic-naive", action="store_true", help="for --engine heuristic: eSOL-naive variant (no MW/pI)")
    ap.add_argument("--heuristic-continuous", action="store_true", help="for --engine heuristic: un-rounded float score (fewer AUROC ties; regenerate numbers if enabled)")
    ap.add_argument("--pred-csv", default="", help="for --engine external: external tool's prediction CSV")
    ap.add_argument("--pred-label", default="external predictor", help="name of the external predictor")
    ap.add_argument("--export-fasta", default="", help="write (allowlist-filtered) SoluProt sequences as FASTA (header=row index) and exit")
    _a = ap.parse_args()
    run(out_json=_a.out, save_preds=_a.save_preds, n_boot=_a.bootstrap, allowlist=_a.allowlist, engine=_a.engine,
        pred_csv=_a.pred_csv, pred_label=_a.pred_label, export_fasta=_a.export_fasta, heuristic_naive=_a.heuristic_naive,
        heuristic_continuous=_a.heuristic_continuous)
