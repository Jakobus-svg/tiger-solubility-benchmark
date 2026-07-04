"""
heuristic.py — the TiGer composition solubility heuristic, as a standalone function.

This is the deliberately simple, deterministic composition score used as a baseline in the
generalisation benchmark. It is extracted VERBATIM from the production tool (TiGer Biotech) so
that the baseline is fully reproducible and inspectable, independent of the rest of the tool.

Two variants (selected by `include_mw_pi`):
  * include_mw_pi=False  -> the eSOL-NAIVE variant: only literature-derived composition terms
    (Wilkinson & Harrison 1991 charge term; Idicula-Thomas & Balaji 2005 propensity; aliphatic /
    aromatic / cysteine terms). No parameter is fit to any benchmark. This is the leakage-free
    baseline reported in the paper (eSOL cytoplasmic AUROC 0.723).
  * include_mw_pi=True   -> the FULL variant: additionally adds molecular-weight and isoelectric-
    point terms whose sign/magnitude were informed by an OLS fit to eSOL itself. Because it is
    partly fit to eSOL, its eSOL number (0.787) is partly circular and is NOT claimed; it is kept
    only to quantify, by ablation, how much an eSOL-fit prior inflates apparent eSOL performance.

    Reported eSOL cytoplasmic AUROC for the naive variant: 0.723 on the RP3Net-clean set, 0.725 on
    the doubly-clean (RP3Net ∩ PSI:Biology) intersection used for the paper's native comparisons.

The score is a clipped linear combination on an 8-95 scale (higher = more soluble). The default
integer form reproduces the published 0.723 / 0.787 AUROCs exactly, but the rounding collapses the
score onto ~46 distinct values and so introduces avoidable ranking ties; pass continuous=True for
the un-rounded float, which is the fairer form for a paired AUROC comparison (see tiger_heuristic_score).

References for the composition terms:
  Wilkinson DL, Harrison RG (1991) Bio/Technology 9(5):443-448.
  Idicula-Thomas S, Balaji PV (2005) Protein Sci 14:582-592.
  eSOL OLS reference for the MW/pI directions: Niwa et al. (2009) PNAS 106:4201.
"""
from collections import Counter

# ── Per-residue solubility propensity (Idicula-Thomas & Balaji 2005, as used by TiGer) ──
_SOL_PROP = {
    'A': 0.0, 'C': -0.3, 'D': 0.5, 'E': 0.5, 'F': -0.6, 'G': 0.2, 'H': 0.1,
    'I': -0.6, 'K': 0.5, 'L': -0.6, 'M': -0.1, 'N': 0.3, 'P': 0.0, 'Q': 0.3,
    'R': 0.5, 'S': 0.3, 'T': 0.2, 'V': -0.4, 'W': -0.8, 'Y': -0.5,
}

# ── Residue monoisotopic-ish MW (Da) for the optional MW term ──
AA_MW = {'A': 89.0935, 'R': 174.2017, 'N': 132.1184, 'D': 133.1032, 'C': 121.1590,
         'Q': 146.1451, 'E': 147.1299, 'G': 75.0669, 'H': 155.1552, 'I': 131.1736,
         'L': 131.1736, 'K': 146.1882, 'M': 149.2124, 'F': 165.1900, 'P': 115.1310,
         'S': 105.0930, 'T': 119.1197, 'W': 204.2262, 'Y': 181.1894, 'V': 117.1469}

# ── Bjellqvist/ExPASy pKa set for the optional pI term (ProtParam-equivalent) ──
PKA = {'D': 4.05, 'E': 4.45, 'H': 5.98, 'C': 9.00, 'Y': 10.00, 'K': 10.00, 'R': 12.00}
PKA_NTERM, PKA_CTERM = 7.5, 3.55
PKA_NTERM_RES = {'A': 7.59, 'M': 7.00, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.70}
PKA_CTERM_RES = {'D': 4.55, 'E': 4.75}


def molecular_weight_kda(seq: str) -> float:
    """Average MW in kDa (water-corrected), matching the tool's `mw`."""
    return round((18.0153 + sum(AA_MW.get(aa, 111) - 18.0153 for aa in seq)) / 1000, 2)


def isoelectric_point(seq: str) -> float:
    """pI via Henderson-Hasselbalch + bisection (Bjellqvist/ProtParam), matching the tool's `pi`."""
    if not seq:
        return 7.0
    comp = Counter(seq)
    nterm = PKA_NTERM_RES.get(seq[0], PKA_NTERM)
    cterm = PKA_CTERM_RES.get(seq[-1], PKA_CTERM)

    def charge_at(pH):
        ch = 1.0 / (1.0 + 10.0 ** (pH - nterm))
        ch -= 1.0 / (1.0 + 10.0 ** (cterm - pH))
        for aa, pka in PKA.items():
            c = comp.get(aa, 0)
            if not c:
                continue
            if aa in ('D', 'E', 'C', 'Y'):
                ch -= c / (1.0 + 10.0 ** (pka - pH))
            else:
                ch += c / (1.0 + 10.0 ** (pH - pka))
        return ch

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if charge_at(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def tiger_heuristic_score(seq: str, include_mw_pi: bool = True, continuous: bool = False):
    """
    TiGer composition solubility score on an 8-95 scale (higher = more soluble).

    include_mw_pi=False -> eSOL-naive baseline (literature composition terms only).
    include_mw_pi=True  -> full variant (adds the eSOL-informed MW/pI terms).
    continuous=False (default) -> integer score (verbatim production output; reproduces the
        published AUROCs 0.723/0.725/0.787 exactly). The integer scale collapses the score onto
        ~46 distinct values, so for AUROC/ranking it introduces avoidable mid-rank ties that
        handicap the baseline relative to continuous predictors (SWI, NetSolP).
    continuous=True -> return the raw float sol_f (clipped to 8-95, NOT rounded). For AUROC only
        the ranking matters; the continuous score removes the discretisation ties and is the
        recommended form for a fair paired comparison. Point estimates shift only marginally, but
        callers that switch to this MUST regenerate all prediction files and paper numbers, since
        the committed results/ JSONs were produced with the integer (default) form.

    Verbatim port of the production formula (TiGer main.py); see module docstring.
    """
    seq = (seq or "").upper()
    n = len(seq)
    if n == 0:
        return 8.0 if continuous else 8
    comp = Counter(seq)
    cys = comp.get('C', 0)

    charge_frac = (comp.get('D', 0) + comp.get('E', 0) + comp.get('K', 0) + comp.get('R', 0)) / n
    hydrophob_frac = (comp.get('I', 0) + comp.get('L', 0) + comp.get('V', 0) + comp.get('M', 0)) / n
    aliphatic_frac = (comp.get('I', 0) + comp.get('L', 0) + comp.get('V', 0)) / n
    aromatic_frac = (comp.get('W', 0) + comp.get('F', 0) + comp.get('Y', 0)) / n
    aa_sol_score = sum(_SOL_PROP.get(aa, 0.0) for aa in seq) / n
    wh_insol = (hydrophob_frac + aromatic_frac) > 0.45 or aliphatic_frac > 0.40
    wh_soluble = charge_frac >= 0.24 and aliphatic_frac < 0.35

    sol_f = 50.0
    sol_f += charge_frac * 55       # Wilkinson & Harrison 1991 (charged residues -> soluble)
    sol_f -= hydrophob_frac * 45    # aliphatic hydrophobic (I/L/V/M)
    sol_f -= aromatic_frac * 30     # aromatic (W/F/Y) - pi-stacking aggregation
    sol_f += aa_sol_score * 20      # Idicula-Thomas & Balaji 2005 propensity
    sol_f -= (cys > 3) * 7          # disulfide complexity
    sol_f -= (cys > 6) * 5
    sol_f -= 12 if wh_insol else 0
    sol_f += 8 if wh_soluble else 0
    if include_mw_pi:
        # eSOL-OLS-informed terms (omit for the leakage-free naive baseline)
        mw = molecular_weight_kda(seq)
        pi = isoelectric_point(seq)
        sol_f -= max(-8, min(15, (mw - 25) * 0.30))   # MW: larger -> less soluble (softened -0.30)
        sol_f += max(-6, min(8, (6.5 - pi) * 1.5))    # pI: away from 6.5 -> less soluble (softened)
    sol_f -= 15 if n < 50 else 0    # very short peptides -> inclusion-body risk
    sol_f = min(95.0, max(8.0, sol_f))
    return sol_f if continuous else int(sol_f)


if __name__ == "__main__":
    # quick self-check on a couple of sequences
    demo = "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQ"
    print("full :", tiger_heuristic_score(demo, include_mw_pi=True))
    print("naive:", tiger_heuristic_score(demo, include_mw_pi=False))
