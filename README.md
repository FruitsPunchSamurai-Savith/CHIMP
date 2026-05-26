# CHIMP

**CHIpr Many-body Polynomial fitter** — a full-dimensional potential energy surface (PES) fitter for triatomic systems, combining the **Aguado–Paniagua (AP)** many-body expansion with key improvements from **CHIPR** (Combined-Hyperbolic-Inverse-Power-Representation). Implemented as a self-contained Google Colab notebook.

> **Data not included.** Ab initio datasets are proprietary and cannot be distributed. See [Bring Your Own Data](#bring-your-own-data) for the expected format.

---

## Overview

CHIMP constructs a full-dimensional PES as follows:

$$V(R_{12}, R_{13}, R_{23}) = V_2^{AB} + V_2^{AC} + V_2^{BC} + V_3 + V_\text{lr} + E_0$$

where 2-body terms are fit with either the AP or CHIPR diatomic basis, and the 3-body term $V_3$ is expanded in symmetry-adapted polynomials with nonlinear range parameters $\beta$ optimized by a grid scan + Levenberg–Marquardt refinement.

### Key features

- **Dual diatomic basis**: classic AP (`rho = R·exp(−β·R)`) or CHIPR sech/csech basis — the latter strongly preferred for ionic/polar systems (LiF, NaH)
- **REPDAMP**: CHIPR-inspired short-range damping that forces $V_3 \to 0$ at compressed bonds, preventing polynomial waste on the repulsive wall
- **Energy weighting**: cubic down-weighting of high-energy points via `(E_ref / (E_ref + E_rel))³`, tunable per system
- **Angular coverage boost**: sparse-angle reweighting in the potential well (used for LiFH)
- **Long-range corrections**: ion–induced-dipole (`−α/2R⁴`) for ionic systems; switched dipole–dipole (`−C₆/R⁶`) for polar systems
- **Stratified RMSD**: error breakdown by energy stratum (0–5, 5–10, … kcal/mol above minimum)
- **Symmetry-adapted 3-body monomials**: exact gfit3c index tables for AAA (`indice=1`), ABB (`indice=2`), and ABC (`indice=3`) symmetries

---

## Repository Structure

```
CHIMP/
├── PES_fit_CHIPR_2_.ipynb   ← main notebook (14 cells, run top-to-bottom)
└── README.md
```

---

## Bring Your Own Data

Data files are not included. You need to supply them and point Cell 4's `CONFIGS` dict to the right paths. Expected formats:

**Diatomic** — two-column whitespace-delimited, R in bohr, E in hartree (absolute energy):
```
1.200   -0.987654
1.400   -0.991234
...
```

**Triatomic** — four-column, distances in bohr, energy in hartree. Column order follows the gfit3c convention (`R12, R23, R13, E`); CHIMP re-orders to `R12, R13, R23` internally:
```
1.80   2.40   2.10   -1.499123
...
```

For NaH2-style systems where diatomic and triatomic data are embedded in a `gfit3c.inp` file, set `triat_loader = 'inp'` in the config and use Cell 3 to extract the diatomic blocks.

---

## Quick Start (Colab)

```python
# 1. Place your data files in Google Drive under a folder (e.g. MyDrive/PES_data/)

# 2. Open the notebook in Colab and run Cell 1 to mount Drive

# 3. In Cell 4, update BASE and file paths, then set:
SYSTEM = 'H3'   # or 'FH2' | 'LiFH' | 'NaH2' | your own key

# 4. Run all cells (Runtime → Run all)
```

Output includes:
- Diatomic fit parameters (β₁, β₂ or CHIPR γ/ζ) and RMS
- 3-body fit β values, RMS, E_max
- Stratified RMSD table
- Four plots per system: error vs index, error vs energy, PES contours, ab initio vs fitted overlays

---

## Benchmark Results

Results on the development systems (data not distributed):

| System | Type | Basis | N points | RMS (kcal/mol) | E_max (kcal/mol) |
|--------|------|-------|----------|----------------|-----------------|
| H3     | AAA  | AP    | 25 804   | 0.447          | 10.31           |
| FH2    | ABB  | AP    | 11 160   | 0.787          | 6.02            |
| LiFH   | ABC  | CHIPR | 16 329   | 4.017 (3b only)| 95.12           |

### Stratified RMSD — H3 (SS-MRPT2)

| Stratum (kcal) | N pts | RMS (kcal/mol) | E_max (kcal/mol) |
|---------------|-------|----------------|-----------------|
| 0–5           | 2034  | 0.610          | 1.866           |
| 5–10          | 1205  | 0.548          | 1.865           |
| 10–20         | 2186  | 0.603          | 2.105           |
| 20–50         | 3441  | 0.607          | 6.775           |
| 50–100        | 8196  | 0.470          | 10.309          |
| 100–200       | 8742  | 0.142          | 6.930           |

### Stratified RMSD — FH2

| Stratum (kcal) | N pts | RMS (kcal/mol) | E_max (kcal/mol) |
|---------------|-------|----------------|-----------------|
| 0–5           | 472   | 0.358          | 0.829           |
| 5–10          | 526   | 0.354          | 1.099           |
| 10–20         | 1015  | 0.337          | 1.436           |
| 20–50         | 4759  | 0.696          | 3.753           |
| 50–100        | 3683  | 0.735          | 2.949           |
| 100–200       | 705   | 1.836          | 6.021           |

### Stratified RMSD — LiFH (100 kcal/mol cutoff; 3-body RMS)

| Stratum (kcal) | N pts | RMS (kcal/mol) | E_max (kcal/mol) |
|---------------|-------|----------------|-----------------|
| 0–5           | 107   | 2.735          | 7.071           |
| 5–10          | 1687  | 2.530          | 17.855          |
| 10–20         | 2239  | 4.017          | 31.258          |
| 20–50         | 4670  | 8.046          | 111.060         |
| 50–100        | 7626  | 11.024         | 239.477         |

---

## Notebook Walkthrough

| Cell | Contents |
|------|----------|
| 1    | Mount Google Drive |
| 2    | Imports (`numpy`, `scipy`, `matplotlib`) |
| 3    | Extract NaH2 embedded diatomic data from `gfit3c.inp` |
| 4    | **Config** — set `SYSTEM` and file paths here |
| 5    | Data loaders (`load_diatomic`, `load_triatomic`, `load_triatomic_from_inp`) |
| 6    | Diatomic basis: AP (`build_ap_basis`, `fit_ap`) and CHIPR (`build_chipr_basis`, `fit_chipr`) |
| 7    | Long-range corrections (`vlr_ion_induced_dipole`, `vlr_dipole_induced_dipole`) |
| 8    | 3-body monomial index tables (gfit3c IPA31/IPA32/IPA33) + `repdamp` + `build_3body_basis` |
| 9    | `fit_3body`: grid scan → LM refinement on subset → final weighted lstsq on full data |
| 10   | `print_stratified_rmsd` |
| 11   | `eval_pes`: full PES evaluator (2-body + 3-body + long-range) |
| 12   | `run_fit`: main pipeline |
| 13   | `plot_all`: contour plots, error-vs-index, error-vs-energy, ab initio vs fitted overlays |
| 14   | **Run cell** — executes `run_fit` + `plot_all` for selected system |

---

## Configuration Reference

All system parameters live in the `CONFIGS` dict (Cell 4):

| Field | Description |
|-------|-------------|
| `indice` | Symmetry: 1=AAA, 2=ABB, 3=ABC |
| `diat_basis` | `'AP'` or `'CHIPR'` for the AB diatomic |
| `chipr_AB` | Dict with `M`, `gama`, `R_ref0`, `zeta`, `ZA`, `ZB` (CHIPR only) |
| `MT_diat` / `MT_triat` | Polynomial degree cutoff index (maps to n_terms via NFAS tables) |
| `e0_AB/BC/AC/ABC` | Atomic/diatomic reference energies in hartree |
| `energy_cutoff_kcal` | Discard points above this threshold (kcal above minimum) |
| `energy_weight_ref_kcal` | Reference energy for cubic weighting; `None` = uniform |
| `angular_boost` | `True` to upweight sparsely-sampled angles in the well |
| `polar` / `C6_AB` | Enable switched dipole–dipole long-range correction |
| `charged` / `alpha_neutral` | Enable ion–induced-dipole long-range correction |
| `init_vex_triat` | Initial β guesses for 3-body grid scan |

---

## Method Notes

### Diatomic fitting (AP)

$$V_2(R) = \frac{c_0 e^{-\beta_2 R}}{R} + \sum_{i=1}^{M-1} c_i \left(R e^{-\beta_1 R}\right)^i$$

Grid scan over ±0.5 in β space (9 points), then TRF refinement with bounds [0.05, 6.0].

### Diatomic fitting (CHIPR)

$$V_2(R) = \frac{Z_A Z_B}{R} \sum_{j=0}^{M-1} c_j \,\phi_j(R), \quad \phi_j = \operatorname{sech}\!\left(\gamma_j (R - \zeta R_0^j)\right)$$

with a long-range `csech` function for $j = M-1$. Nonlinear parameters (γ, R₀, ζ) are fixed from a separate optimisation; linear coefficients by lstsq.

### 3-body fitting

Range-scaling coordinates: $\rho_{ij} = R_{ij} e^{-\beta_{ij} R_{ij}}$

Fit proceeds in three phases:
1. **Grid scan** on a random 8000-pt subset (±0.5 steps in each β)
2. **LM refinement** on the same subset
3. **Final lstsq** on full dataset with optional energy weights

REPDAMP multiplies every row of the design matrix by $\prod_i \left[\frac{1+\tanh(\kappa(R_i - R_0))}{2}\right]^\xi$ (default κ=100, ξ=10, R₀=0.5 bohr).

---

## Dependencies

```
numpy
scipy
matplotlib
```

Runs on standard Colab CPU runtime, no additional installs needed.

---

## References

- Aguado, A. & Paniagua, M. *J. Chem. Phys.* **96**, 1265 (1992) — AP many-body expansion
- Varandas, A. J. C. *J. Chem. Phys.* **138**, 054120 (2013) — CHIPR method
- Werner, H.-J. et al. — gfit3c monomial index convention
