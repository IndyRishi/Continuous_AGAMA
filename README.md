# Constraining the Inner Dark-Matter Slope of Sculptor

A chemo-dynamical pipeline that measures the **dark-matter inner density slope** of the
Sculptor (and Fornax) dwarf spheroidal galaxies across several independent modeling
frameworks, and tests a central methodological question:

> Is the standard **two-population split** of a dwarf's stars (metal-rich vs metal-poor)
> justified, or is the galaxy a **continuous metallicity–kinematics sequence** whose
> arbitrary splitting biases the inferred DM slope?

The pipeline reproduces the **Walker & Peñarrubia (2011)** mass-slope method and compares
against the **Arroyo-Polonio et al. (2025, "AP25")** action-based analysis on real
spectroscopic data. It also applies a continuous **f(J, [Fe/H])** action-distribution-function
model in which metallicity is a coordinate inside a single DF rather than a label for two
populations. **This continuous model is an adaptation of the extended distribution-function
framework of Sanders & Binney (2015) — developed for the Milky Way disc — to a
pressure-supported dwarf spheroidal; it is not a new DF construction.** The contribution is the
application to the dSph inner-slope problem and the explicit test of the two-population bias.

---

## Key results (real Sculptor data, Tolstoy et al. 2023; 1339 members)

| Framework | Method | Inner slope |
|---|---|---|
| Spherical Jeans (`--dm5`) | gNFW, σ_los only | γ ≈ 0.78 (+0.29 / −0.39) |
| GravSphere (`--gravsphere`) | Jeans + Virial Shape Parameters + free β(r) | γ ≈ 0.48 (+0.41 / −0.29) |
| Walker & Peñarrubia 2011 (`--wp11`) | Two-subcomponent mass slope | Γ = 2.50 (+0.23 / −0.20); NFW excluded 99.6% |
| Continuous f(J,[Fe/H]) (`--continuous`) | Single-DF, metallicity as coordinate | mock-validated; real-data gradient detected (k_J < 0) |

All frameworks point to a **core-like** (shallow) inner slope, consistent within
uncertainties with AP25's published γ = 0.39 (+0.23 / −0.26).

**Status note on the continuous result:** the mock recovery is validated (recovers injected γ
and k_J). On real data it detects the metallicity gradient (k_J < 0), independently corroborated
by an action-space analysis using Gaia proper motions. The real-data continuous fit is currently
run on a representative subsample and is convergence-limited; the full-sample converged run
requires cluster-scale compute.

### Evidence that Sculptor is one continuous population

- **Metallicity is statistically unimodal** by three independent tests: bimodality coefficient
  (BC = 0.30 < 0.555), ΔBIC = −8 (favours a single Gaussian), and Hartigan's dip test (p = 0.79).
- **Kinematics vary smoothly** with any imposed metallicity split (no plateau).
- **σ_los alone is degenerate** in γ (a core-through-cusp range fits equally well), which is why
  slopes are reported as MCMC posteriors, not point estimates.
- **On mocks, a two-population split biases γ high** (+0.32 to +0.42), while the continuous
  treatment is far closer to truth.
- **Membership robustness:** the WP11 slope is stable across Gaia proper-motion membership
  thresholds, so the result is not driven by the membership selection.

---

## Installation

Requires Python 3.11+. The dynamical modeling uses
[AGAMA](https://github.com/GalacticDynamics-Oxford/Agama).

```bash
pip install -r requirements.txt
```

AGAMA needs a C++ toolchain and GSL; if `pip install agama` fails, build from source (on a
cluster: `module load gsl`, then `pip install agama`). Real-data commands query
[VizieR](https://vizier.cds.unistra.fr/) and Gaia DR3, so they require network access; mock
commands run offline.

---

## Usage

All analyses are flags on the single entry-point script. Add `--galaxy fornax` to any real-data
command to switch targets (outputs are automatically prefixed `fornax_`).

### Data presentation and evidence (fast)
```bash
python sculptor_full_pipeline.py --figure1        # two-galaxy data presentation (Sculptor + Fornax)
python sculptor_full_pipeline.py --overview       # data overview + [Fe/H] unimodality stats
python sculptor_full_pipeline.py --slide          # sliding-threshold σ_los + γ-degeneracy curve
python sculptor_full_pipeline.py --biasgate       # two-population split biases γ (mocks)
python sculptor_full_pipeline.py --biasconv       # mean bias vs number of realizations
```

### Dark-matter slope measurements (MCMC; resume by re-running)
```bash
python sculptor_full_pipeline.py --dm5            # spherical-Jeans gNFW inner slope
python sculptor_full_pipeline.py --gravsphere     # GravSphere: Jeans + VSPs + free β(r)
python sculptor_full_pipeline.py --wp11           # Walker & Peñarrubia 2011 (reproduces their Fig. 10)
python sculptor_full_pipeline.py --continuous     # continuous f(J,[Fe/H]) model
python sculptor_full_pipeline.py --chain          # full 25-parameter AP25 action-DF (cluster-scale)
```

### Robustness, comparison, and Gaia-based figures
```bash
python sculptor_full_pipeline.py --robustness     # WP11 slope vs Gaia PM membership threshold
python sculptor_full_pipeline.py --compare        # γ posteriors across frameworks
python sculptor_full_pipeline.py --fig4all        # DM density ρ(r) from all chains + AP25's curve
python sculptor_full_pipeline.py --dispersion     # radial σ_los profile (add --gaia for PM dispersion)
python sculptor_full_pipeline.py --skymap         # Gaia PM sky map with membership-cut panels
python sculptor_full_pipeline.py --actions        # action-space chemodynamics ([Fe/H] vs J_r)
```

### Validation / cross-checks
```bash
python sculptor_full_pipeline.py --continuous --mock     # recover known γ, k_J (validates the method)
python sculptor_full_pipeline.py --wp11 --mock           # WP11 recovery on a known-Γ mock
python sculptor_full_pipeline.py --crosscheck --repo <path/to/gravsphere>   # vs reference GravSphere
```

### Common flags
`--steps N` (MCMC steps to add), `--walkers N`, `--nproc N` (parallel processes), `--K N`
(metallicity DF nodes for `--continuous`), `--nsub N` (fit a random subsample; omit for the full
sample), `--backend FILE.h5` (checkpoint), `--no-resume`, `--galaxy {sculptor,fornax}`, `--mock`.

**MCMC runs checkpoint every step** to an HDF5 backend and **resume** on re-launch — run the same
command until the convergence report stops reporting `NOT CONVERGED`. Check progress from a
second terminal:
```bash
python -c "import emcee; b=emcee.backends.HDFBackend('wp11.h5'); print('steps:', b.iteration)"
```

---

## Outputs

### Chains (`.npy`, posterior samples)
`dm5_chain.npy`, `gravsphere_chain.npy`, `wp11_chain.npy`, `cont_chain.npy`,
`ap25_chain.npy` (Fornax runs prefix `fornax_`).

### Figures
| File | Content |
|---|---|
| `figure1_two_galaxy.png` | Two-galaxy data presentation (Sculptor + Fornax) |
| `figure_data_overview.png` | Sample overview + [Fe/H] unimodality (BC, ΔBIC, dip test) |
| `figure_sliding_metallicity.png` | σ_los continuum + σ_los γ-degeneracy curve |
| `figure_bias_gate.png` | Two-population split biases γ (mocks) |
| `figure_bias_vs_realizations.png` | Bias convergence vs number of realizations |
| `figure_gravsphere_beta.png` | GravSphere anisotropy profile β(r) |
| `figure_wp11.png` | WP11 Fig. 10 reproduction (Γ, NFW exclusion) |
| `figure_membership_robustness.png` | WP11 slope stable across membership thresholds |
| `figure_continuous_corner.png` | Continuous f(J,[Fe/H]) posterior |
| `figure_action_space.png` | Action-space chemodynamics ([Fe/H] vs J_r) |
| `figure_dispersion_profile.png` | Radial velocity-dispersion profile (LOS + PM) |
| `figure_gaia_skymap.png` | Gaia PM sky map with membership-cut panels |
| `figure_framework_comparison.png` | γ posteriors across frameworks |
| `figure_fig4_all_chains.png` | DM density ρ(r) from all frameworks + AP25's curve |
| `fornax_figure_wp11.png` | Fornax WP11 (failure-case contrast) |

---

## Methods

**Data.** Sculptor: Tolstoy et al. (2023) VizieR catalog `J/A+A/675/A49`, members with reliable
[Fe/H], in the galaxy rest frame, on elliptical (semi-major-axis) radii (Muñoz+2018
centre/ellipticity/PA, D = 84 kpc). Fornax: Walker et al. (2009) MMFS (`J/AJ/137/3100`), a
multi-galaxy catalog filtered to Fornax members (Mmb > 0.5), using the Mg spectral index W′ as
the metallicity separator (as WP11 did), D = 147 kpc.

**Frameworks.**
- *Spherical Jeans* — gNFW halo (α=1, η=3), Osipkov–Merritt anisotropy, fit to the binned
  σ_los(R) of the two metallicity halves.
- *GravSphere* — from-scratch implementation of Read & Steger (2017): spherical Jeans with two
  Virial Shape Parameters and a free Baes–van Hese anisotropy profile β(r), cross-checked against
  the reference `gravsphere` code (σ_los to 0.3%, VSPs to 0.14%).
- *Walker & Peñarrubia (2011)* — two chemo-dynamically distinct stellar subcomponents, each a
  Plummer sphere with Gaussian velocity and metallicity distributions; the estimator
  M(r_h) = 5 r_h σ² / (2G) at two half-light radii gives the slope Γ ≡ Δlog M / Δlog r. Γ > 2
  excludes an NFW cusp.
- *AP25 action-DF* — the full 25-parameter model with the paper's Eq. 4 gNFW-with-cutoff
  potential.
- *Continuous f(J,[Fe/H])* — a single DoublePowerLaw DF whose scale action varies smoothly with
  metallicity, `log₁₀ J₀(z) = logJ₀ + k_J·(z − ⟨z⟩)`; each star's likelihood marginalises over
  its true metallicity. A detected `k_J < 0` means metal-rich stars are centrally concentrated —
  a gradient captured without splitting. **Adapted from the extended-DF framework of Sanders &
  Binney (2015).**

---

## Scope & caveats

- Reported slopes are **MCMC posteriors, not point estimates**, because the σ_los-only
  likelihood is degenerate in γ.
- The continuous method is an **adaptation of extended DFs (Sanders & Binney 2015)**, not a new
  construction; its contribution is the dSph application and the two-population-bias test.
- The continuous real-data result is a **representative subsample** and convergence-limited; the
  full-sample converged run requires cluster-scale compute.
- The full 25-parameter AP25 action-DF (`--chain`) is **cluster-scale**; the reduced frameworks
  serve as the discrete baselines.
- WP11's absolute masses differ from their Table 4 by a constant estimator offset that **cancels
  in Γ**; the slope reproduces their value.
- The AP25 curve overplotted in the Fig. 4 figures is computed from their published parameters
  via their Eq. 4 (not by re-running their method).
- **Fornax is a failure-case contrast, not a measurement:** the two-population method does not
  constrain Fornax's slope (populations do not separate); the sample is small (452 stars) and
  Mg-selection-biased.
- Gaia PM figures (sky map, PM dispersion, action space) are illustrative; the core measurements
  use the spectroscopic sample.

---

## Repository layout

- `sculptor_full_pipeline.py` — the complete pipeline (all frameworks and figures).
- `requirements.txt` — dependencies.
- `crosscheck_gravsphere.py`, `run_phase5.py`, `selection_function.py` — thin
  backward-compatibility shims around the main module.

---

*Developed as part of a research project on dwarf-spheroidal dark matter.*
