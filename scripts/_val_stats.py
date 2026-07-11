"""
_val_stats.py — small, dependency-free statistics for the validation scripts.

Provides:
  * auroc(scores, labels)                         -- Mann-Whitney AUROC
  * bootstrap_auroc_ci(scores, labels, ...)       -- AUROC + percentile 95% CI
  * paired_bootstrap_diff(a, b, labels, ...)      -- AUROC(a)-AUROC(b) on the SAME proteins,
                                                     95% CI of the difference + two-sided p
  * delong_test(a, b, labels)                     -- DeLong (1988) paired AUROC comparison
                                                     (auc_a, auc_b, z, p); independent cross-check

Pure stdlib (random, math). No numpy/scipy required, so it runs anywhere the validation scripts do.
The paired bootstrap is the primary, robust comparison; DeLong is included as a second opinion.
"""
import math
import random


# ── AUROC (Mann-Whitney U) ───────────────────────────────────────────────────
def auroc(scores, labels):
    """AUROC via the rank-sum identity; handles ties with mid-ranks. labels in {0,1}."""
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return float("nan")
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    n = len(scores)
    s_sorted = [scores[order[k]] for k in range(n)]
    while i < n:
        j = i
        while j < n and s_sorted[j] == s_sorted[i]:
            j += 1
        avg = (i + j - 1) / 2.0 + 1.0          # average rank (1-based) for the tie block
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    sum_pos = sum(r for r, y in zip(ranks, labels) if y == 1)
    n_pos, n_neg = len(pos), len(neg)
    return (sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


# ── Spearman rank correlation (threshold-free complement to AUROC) ───────────
def spearman(scores, target):
    """Spearman's rho between predicted scores and a CONTINUOUS target (e.g. eSOL % solubility).

    A threshold-free complement to AUROC: it uses the whole solubility scale rather than a binarised
    label, so it checks that the AUROC ranking is not an artefact of the chosen cutoff. Ties are
    handled with mid-ranks. Returns NaN for degenerate input. Pure stdlib."""
    if len(scores) != len(target) or len(scores) < 2:
        return float("nan")

    def _rank(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i, n = 0, len(x)
        xs = [x[order[k]] for k in range(n)]
        while i < n:
            j = i
            while j < n and xs[j] == xs[i]:
                j += 1
            avg = (i + j - 1) / 2.0
            for k in range(i, j):
                r[order[k]] = avg
            i = j
        return r

    a, b = _rank(scores), _rank(target)
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    da = sum((x - ma) ** 2 for x in a) ** 0.5
    db = sum((x - mb) ** 2 for x in b) ** 0.5
    return num / (da * db) if da > 0 and db > 0 else float("nan")


# ── bootstrap CI for a single AUROC ──────────────────────────────────────────
def bootstrap_auroc_ci(scores, labels, n_boot=2000, alpha=0.05, seed=0):
    """Return (auroc, lo, hi) — point estimate and percentile (1-alpha) CI by case resampling."""
    a = auroc(scores, labels)
    if math.isnan(a):
        return a, float("nan"), float("nan")
    rng = random.Random(seed)
    n = len(scores)
    idx = range(n)
    boots = []
    for _ in range(n_boot):
        bi = [rng.randint(0, n - 1) for _ in idx]
        bs = [scores[k] for k in bi]
        bl = [labels[k] for k in bi]
        v = auroc(bs, bl)
        if not math.isnan(v):
            boots.append(v)
    if not boots:
        return a, float("nan"), float("nan")
    boots.sort()
    lo = boots[int((alpha / 2) * len(boots))]
    hi = boots[min(len(boots) - 1, int((1 - alpha / 2) * len(boots)))]
    return a, lo, hi


# ── paired bootstrap difference of two AUROCs on the same samples ────────────
def paired_bootstrap_diff(scores_a, scores_b, labels, n_boot=2000, alpha=0.05, seed=0):
    """AUROC(a) - AUROC(b) on the SAME proteins (paired). Returns
    (diff, lo, hi, p_two_sided).

    NOTE ON THE P-VALUE. This is a bootstrap achieved-significance-level (ASL), read off the
    resampling distribution of the difference (which is centred on the observed diff, not on the
    null). It is a ROBUSTNESS CROSS-CHECK, not a null-hypothesis test; delong_test() is the primary,
    parametric significance test and should carry the headline p-values. The ASL cannot resolve a
    p-value below 1/n_boot: when no replicate crosses zero the returned p is 0.0, which means only
    "p < 1/n_boot" (e.g. < 5e-4 at n_boot=2000) — it does NOT justify reporting "p < 0.0001". Use
    DeLong (or raise n_boot) for small p-values. seed is fixed for reproducibility; the CI/p are one
    Monte-Carlo realisation, so use n_boot >= 10000 for final reported figures.

    p = 2*min(P(diff*<=0), P(diff*>=0)) over bootstrap replicates."""
    da = auroc(scores_a, labels)
    db = auroc(scores_b, labels)
    diff = da - db
    rng = random.Random(seed)
    n = len(labels)
    diffs = []
    for _ in range(n_boot):
        bi = [rng.randint(0, n - 1) for _ in range(n)]
        la = [labels[k] for k in bi]
        va = auroc([scores_a[k] for k in bi], la)
        vb = auroc([scores_b[k] for k in bi], la)
        if not (math.isnan(va) or math.isnan(vb)):
            diffs.append(va - vb)
    if not diffs:
        return diff, float("nan"), float("nan"), float("nan")
    diffs.sort()
    lo = diffs[int((alpha / 2) * len(diffs))]
    hi = diffs[min(len(diffs) - 1, int((1 - alpha / 2) * len(diffs)))]
    n_le = sum(1 for d in diffs if d <= 0)
    n_ge = sum(1 for d in diffs if d >= 0)
    p = 2.0 * min(n_le, n_ge) / len(diffs)
    # p == 0 only means "below the resolution of this bootstrap" (1/len(diffs)); never report it as
    # an exact tiny number — defer to DeLong for small p. See docstring.
    return diff, lo, hi, min(1.0, p)


# ── DeLong test (fast, midrank) — independent cross-check ────────────────────
def _midrank(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    rank = [0.0] * len(x)
    i = 0
    n = len(x)
    xs = [x[order[k]] for k in range(n)]
    while i < n:
        j = i
        while j < n and xs[j] == xs[i]:
            j += 1
        avg = (i + j - 1) / 2.0 + 1.0
        for k in range(i, j):
            rank[order[k]] = avg
        i = j
    return rank


def delong_test(scores_a, scores_b, labels):
    """DeLong (1988) paired comparison of two AUROCs on the same samples.
    Returns (auc_a, auc_b, z, p_two_sided). NaN p if variance is degenerate."""
    pos = [i for i, y in enumerate(labels) if y == 1]
    neg = [i for i, y in enumerate(labels) if y == 0]
    m, n = len(pos), len(neg)
    if m == 0 or n == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")

    def structural(scores):
        xp = [scores[i] for i in pos]
        xn = [scores[i] for i in neg]
        tx = _midrank(xp)
        ty = _midrank(xn)
        tz = _midrank(xp + xn)
        auc = (sum(tz[:m]) - m * (m + 1) / 2.0) / (m * n)
        v01 = [(tz[k] - tx[k]) / n for k in range(m)]          # over positives
        v10 = [1.0 - (tz[m + k] - ty[k]) / m for k in range(n)]  # over negatives
        return auc, v01, v10

    auc_a, a01, a10 = structural(scores_a)
    auc_b, b01, b10 = structural(scores_b)

    def cov(p, q, k):
        mp = sum(p) / k
        mq = sum(q) / k
        return sum((p[i] - mp) * (q[i] - mq) for i in range(k)) / (k - 1) if k > 1 else 0.0

    s01_aa, s01_bb, s01_ab = cov(a01, a01, m), cov(b01, b01, m), cov(a01, b01, m)
    s10_aa, s10_bb, s10_ab = cov(a10, a10, n), cov(b10, b10, n), cov(a10, b10, n)
    var_a = s01_aa / m + s10_aa / n
    var_b = s01_bb / m + s10_bb / n
    cov_ab = s01_ab / m + s10_ab / n
    se = math.sqrt(max(var_a + var_b - 2 * cov_ab, 0.0))
    if se == 0:
        return auc_a, auc_b, float("nan"), float("nan")
    z = (auc_a - auc_b) / se
    p = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))
    return auc_a, auc_b, z, p
