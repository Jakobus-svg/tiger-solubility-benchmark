# Input datasets

The committed files here are the **deduplicated id allowlists** used for leak-free evaluation:

- `esol_clean_bnumbers.txt` — eSOL b-numbers with near-identical matches to RP3Net's public SGC
  training data removed (4/3,131 = 0.13% removed; the **RP3Net-clean** set).
- `esol_clean_both.txt` — the **doubly-clean** set (RP3Net-clean ∩ PSI:Biology-clean), 3,000
  b-numbers; this defines the native-protein table (2,154 cytoplasmic). Derived from the leak-free
  prediction files; use it to reproduce the exact table values without a NetSolP-dedup re-run.
- `soluprot_clean_idx.txt` — SoluProt test-set row indices after removing 40/3,100 (1.29%)
  overlapping RP3Net's training data. **Not committed** (row indices depend on the raw CSV);
  regenerate offline with `scripts/regen_soluprot_clean_idx.py` (from the committed clean predictions),
  or from scratch with `scripts/dedup_soluprot_vs_rp3net.py`.
- `esol_clean_psi_homology.txt` — eSOL b-numbers surviving the **homology-level** screen (MMseqs2
  clustering with PSI:Biology training at ≤30% identity; 256/3,000 = 8.5% removed), 2,744 b-numbers.
  Used to show the native ranking is homology-robust (see `results/esol_psi_homology_dedup.json` and
  `scripts/subset_by_allowlist.py`).

**Note on eSOL→UniProt mapping.** The mapping is not strictly 1:1: two b-numbers (b0621, b2094) appear
on multiple eSOL rows. They are de-duplicated on load (keep first); b2094 has *conflicting* solubility
across its source rows (72% vs 14%) and is flagged by the loader — verify against the raw eSOL
supplement before treating it as ground truth. The near-duplicate screen is at ≥90% identity; for a
homology-level control use `scripts/dedup_homology.py` (≤30% identity clustering).

The raw input datasets are **not redistributed** (third-party licenses). Obtain and place here:

| file | columns expected | source |
|---|---|---|
| `esol.csv` | `B number`, `Solubility (%)`, `Cell location` | Niwa et al. 2009 PNAS 106:4201 (supplementary) |
| `ecoli_k12_proteome.tsv` | UniProt TSV incl. `Sequence`, `Gene Names (ordered locus)` | UniProt UP000000625 |
| `soluprot_test.csv` | `seqs`, `labels` | Hon et al. 2021 Bioinformatics 37(1):23 (SoluProt test set) |
