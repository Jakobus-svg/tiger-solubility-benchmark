"""
swi.py — Solubility-Weighted Index (SWI), an eSOL-INDEPENDENT composition baseline.

Exact weights and formula from the authors' published implementation:
  Bhandari BK, Gardner PP, Lim CS (2020) Bioinformatics 36(18):4691-4698,
  doi:10.1093/bioinformatics/btaa578
  source: github.com/Gardner-BinfLab/SoDoPE_paper_2020/SWI/swi.py

SWI weights were optimised on the PSI:Biology dataset and VALIDATED on eSOL as an independent
set, so SWI is a clean (eSOL-naive) second baseline for the eSOL comparison — unlike the in-house
heuristic, whose MW/pI terms were fit to eSOL. Higher SWI = more soluble; AUROC uses the raw SWI
(monotonic with the probability, so the ROC is identical).
"""
import math

# Exact per-residue weights (verbatim from the authors' swi.py)
SWI_WEIGHTS = {
    'A': 0.8356471476582918, 'C': 0.5208088354857734, 'E': 0.9876987431418378,
    'D': 0.9079044671339564, 'G': 0.7997168496420723, 'F': 0.5849790194237692,
    'I': 0.6784124413866582, 'H': 0.8947913996466419, 'K': 0.9267104557513497,
    'M': 0.6296623675420369, 'L': 0.6554221515081433, 'N': 0.8597433107431216,
    'Q': 0.789434648348208,  'P': 0.8235328714705341, 'S': 0.7440908318492778,
    'R': 0.7712466317693457, 'T': 0.8096922697856334, 'W': 0.6374678690957594,
    'V': 0.7357837119163659, 'Y': 0.6112801822947587,
}
# Logistic constants from the authors' fit: prob = 1 / (1 + exp(-(A*SWI + B)))
_A, _B = 81.0581, -62.7775


def swi_score(seq: str) -> float:
    """Mean per-residue SWI weight (U treated as C, per the authors' reader). NaN if no scorable residues."""
    seq = (seq or "").upper().replace("U", "C")
    vals = [SWI_WEIGHTS[a] for a in seq if a in SWI_WEIGHTS]
    return sum(vals) / len(vals) if vals else float("nan")


def swi_prob(seq: str) -> float:
    """Authors' logistic probability of solubility from SWI (monotonic with swi_score)."""
    s = swi_score(seq)
    if s != s:  # NaN
        return float("nan")
    return 1.0 / (1.0 + math.exp(-(_A * s + _B)))
