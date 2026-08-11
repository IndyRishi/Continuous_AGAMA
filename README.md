# Constraining the Inner Dark-Matter Slope of Sculptor

A chemo-dynamical pipeline that measures the **dark-matter inner density slope** of the
Sculptor dwarf spheroidal across four methodologically distinct frameworks, and asks a
methodological question that cuts across all of them:

> The standard **two-population split** of a dwarf's stars (metal-rich vs metal-poor) is
> used to break the mass-anisotropy degeneracy. Does the decomposition itself bias the
> recovered slope -- and if so, is the bias chemical or geometric?

The pipeline reproduces the **Walker & Penarrubia (2011)** mass-slope method on both
Sculptor and Fornax, implements **GravSphere** (Read & Steger 2017) from scratch with a
cross-validation against the reference code, and adds a continuous **f(J, [Fe/H])**
action-distribution-function model in which metallicity is a coordinate inside a single DF
rather than a label for two populations. **That continuous model is an adaptation of the
extended-DF framework of Sanders & Binney (2015) -- developed for the Milky Way disc -- to a
pressure-supported dwarf spheroidal; it is not a new DF construction.** Its contribution is
the coupling parameter `k_J`, which no split-based method can produce.

---

## Key results

Real Sculptor data, Tolstoy et al. (2023), 1339 stars with reliable [Fe/H].

| Framework | Method | Inner slope |
|---|---|---|
| Spherical Jeans (`--dm5`) | gNFW, sigma_los only | gamma = 0.97 (+0.21 / -0.33) |
| GravSphere (`--gravsphere`) | Jeans + Virial Shape Parameters + free beta(r) | gamma = 0.46 (+0.40 / -0.28) |
| Walker & Penarrubia 2011 (`--wp11`) | Two-subcomponent mass slope | Gamma = 2.50 (+0.22 / -0.20) -> gamma = 0.50; 99.6% of the posterior above Gamma = 2 |
| Continuous f(J,[Fe/H]) (`--continuous`) | Single DF, metallicity as coordinate | chain not yet converged; action-space gradient detected independently |

GravSphere and the two-population estimator agree on a **shallow (core-like)** slope,
consistent with the published gamma = 0.39 (+0.23 / -0.26) of Arroyo-Polonio et al. (2025)
and the gamma = 0.5 +/- 0.3 of Zhu et al. (2016) from independent data. The spherical Jeans
fit sits higher at 0.97; a prior-sensitivity test shows this is **not** a prior artifact
(below), so the disagreement between methods is real rather than an inference failure.

### What the mocks actually show

The two-population split is often assumed to bias the slope. Controlled experiments
separate the two candidate causes:

- **Null control** (one Plummer tracer, no metallicity structure): split and all-star fits
  agree within their standard errors (+0.230 vs +0.231 at gamma_true = 0.4).
- **Label-scramble control** (two Plummers, metallicity permuted): likewise (-0.091 vs -0.084).
- Introducing the **two-Plummer geometry alone** moves the recovery by -0.32 to -0.37;
  restoring the true metallicity link moves it back by only +0.04 to +0.11.

So **>=75% of the apparent decomposition bias is tracer geometry, not the decomposition** --
which means mock-based claims about population splitting must control for geometry, this
work's own included.

Which model recovers best depends on which is right: under continuous truth the split
biases high (+0.41, +0.33) while the continuous fit does better (+0.15, +0.07); under
**discrete** truth the ordering reverses (split -0.04, -0.30; continuous -0.12, -0.43).
For Sculptor, which independent photometry and chemodynamics indicate is discrete, the
split is the better-matched model.

### Negative results worth recording

- **Marginal metallicity unimodality does not settle anything.** Sculptor's [Fe/H] is
  unimodal by three tests (BC = 0.30, dBIC = -8, Hartigan dip p = 0.79) -- but two-population
  mocks matched to the literature decomposition return the same, because the components
  overlap at 1.35 sigma. The test cannot discriminate a continuum from two overlapping
  populations.
- **A two-component proper-motion mixture is not identifiable** on a spectroscopically
  preselected sample. The fit returns P_mem = 1 for every star (sd = 8e-13); bounding the
  foreground width relocates the pathology rather than removing it. Membership is therefore
  characterised by direct offset from the systemic PM.
- **The sigma_los-only likelihood is nearly flat in gamma** (dchi2 < 1 across 0 <= gamma
  <= 1.5), which is why slopes are reported as posterior medians. Refitting under a prior
  with median 0.75 rather than 0.95 moves the Jeans posterior by **0.02** -- so the
  posterior is set by the data, not the prior, despite the flat likelihood.

### Fornax as a cross-galaxy validation

Applying the same implementation to the Walker et al. (2009) MMFS Fornax sample recovers
**Gamma = 2.97 (+0.71 / -0.56)** against the published 2.61 (+0.43 / -0.37), with f_mem,
f_sub, the chemical separation and both velocity dispersions agreeing within their
uncertainties.

Reaching that required rebuilding the per-star sample from the per-exposure table (the
averaged Mg-index column is populated for only 459 of 2633 rows, against the 2603 stars
WP11 modelled), retaining non-members with a fitted member fraction, and reconstructing the
spectroscopic completeness w(R) against a Gaia DR3 photometric parent. The residual
difference in half-light radii (~20% larger than theirs) traces to that reconstruction: the
target catalogue behind their w(R) is unpublished.

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

All analyses are flags on the single entry-point script. Add `--galaxy fornax` to any
real-data command to switch targets (outputs are automatically prefixed `fornax_`).

### Data presentation and evidence (fast)

```bash
python continuous_agama.py --figure1        # two-galaxy data presentation (Sculptor + Fornax)
python continuous_agama.py --overview       # data overview + [Fe/H] unimodality stats
python continuous_agama.py --slide          # sliding-threshold sigma_los + gamma-degeneracy curve
python continuous_agama.py --biasgate       # continuous-truth mocks: split vs continuous
python continuous_agama.py --revgate        # discrete-truth mocks: the reverse comparison
python continuous_agama.py --gatediag       # null and label-scramble controls (geometry vs chemistry)
python continuous_agama.py --shrinkage      # prior shrinkage vs chain length
python continuous_agama.py --biasconv       # mean bias vs number of mock realisations
```

`--overview` writes its figure and exits; the others can be combined with `--galaxy`.

### Dark-matter slope measurements (MCMC; resume by re-running)

```bash
python continuous_agama.py --dm5            # spherical-Jeans gNFW inner slope
python continuous_agama.py --gravsphere     # GravSphere: Jeans + VSPs + free beta(r)
python continuous_agama.py --wp11           # Walker & Penarrubia 2011 two-population estimator
python continuous_agama.py --continuous     # continuous f(J,[Fe/H]) model
python continuous_agama.py --chain          # full 25-parameter AP25 action-DF (cluster-scale)
```

### Robustness, comparison, and Gaia-based figures

```bash
python continuous_agama.py --robustness     # WP11 slope vs Gaia PM membership threshold
python continuous_agama.py --pop3           # WP11 slope vs very-metal-poor [Fe/H] floor
python continuous_agama.py --compare        # gamma posteriors across frameworks
python continuous_agama.py --fig4all        # DM density rho(r) from all chains + AP25's curve
python continuous_agama.py --dispersion     # radial sigma_los profile (add --gaia for PM dispersion)
python continuous_agama.py --skymap         # Gaia PM sky map with membership-cut panels
python continuous_agama.py --actions        # action-space chemodynamics ([Fe/H] vs J_r)
python continuous_agama.py --gmmdiag        # PM mixture identifiability diagnostic
```

### Validation / cross-checks

```bash
python continuous_agama.py --continuous --mock     # recover known gamma, k_J
python continuous_agama.py --wp11 --mock           # WP11 recovery on a known-Gamma mock
python continuous_agama.py --crosscheck --repo <path/to/gravsphere>   # vs reference GravSphere
python continuous_agama.py --crosscheck --repo <path> --mcmc          # + end-to-end posterior equivalence
```

### Common flags

`--steps N` (MCMC steps to add), `--walkers N`, `--nproc N` (parallel processes), `--K N`
(metallicity DF nodes for `--continuous`), `--nsub N` (fit a random subsample; omit for the
full sample), `--backend FILE.h5` (checkpoint), `--no-resume`, `--galaxy {sculptor,fornax}`,
`--mock`, `--seed N` (walker initialisation), `--mockseed N` (mock draw), `--gammaprior LO HI`
(Jeans prior bounds, for the prior-sensitivity test), `--astar KPC` (GravSphere tracer scale),
`--selection {none,radial,2d}` (AP24-style completeness for the continuous DF).

**MCMC runs checkpoint every step** to an HDF5 backend and **resume** on re-launch -- run the
same command until the convergence report stops reporting `NOT CONVERGED`. To check progress
from a second terminal, **copy the backend first**; opening the live file will lock it and
crash the running chain:

```bash
cp wp11.h5 peek.h5
python -c "import emcee; b=emcee.backends.HDFBackend('peek.h5',read_only=True); print('steps:', b.iteration)"
```

---

## Outputs

### Chains

HDF5 backends named by `--backend`, plus `.npy` posterior samples: `dm5_chain.npy`,
`gravsphere_chain.npy`, `wp11_chain.npy`, `cont_chain.npy`, `ap25_chain.npy` (Fornax runs
prefix `fornax_`).

### Figures

| File | Content |
|---|---|
| `figure1_two_galaxy.png` | Two-galaxy data presentation (Sculptor + Fornax) |
| `figure_data_overview.png` | Sample overview + [Fe/H] unimodality (BC, dBIC, dip test) |
| `figure_sliding_metallicity.png` | sigma_los continuum + gamma-degeneracy curve |
| `figure_bias_gate.png` | Continuous-truth mocks: split vs continuous |
| `figure_reverse_bias_gate.png` | Discrete-truth mocks: the reverse comparison |
| `figure_gate_diagnostics.png` | Null and scramble controls: geometry vs chemistry |
| `figure_shrinkage_test.png` | Prior shrinkage, flat across chain length |
| `figure_gravsphere_beta.png` | GravSphere anisotropy profile beta(r) |
| `figure_wp11.png` | Sculptor two-population mass slope and Gamma posterior |
| `figure_membership_robustness.png` | Gamma stable across PM membership thresholds |
| `figure_pop3_robustness.png` | Gamma vs very-metal-poor [Fe/H] floor |
| `figure_continuous_corner.png` | Continuous f(J,[Fe/H]) posterior |
| `figure_action_space.png` | Action-space chemodynamics ([Fe/H] vs J_r) |
| `figure_dispersion_profile.png` | Radial velocity-dispersion profile (LOS + PM) |
| `figure_gaia_skymap.png` | Gaia PM sky map with membership-cut panels |
| `figure_fig4_all_chains.png` | DM density rho(r) from all frameworks + AP25's curve |
| `fornax_figure_wp11.png` | Fornax two-population fit (cross-galaxy validation) |

---

## Methods

**Data.** Sculptor: Tolstoy et al. (2023) VizieR catalog `J/A+A/675/A49`, members with
reliable [Fe/H], in the galaxy rest frame, on elliptical (semi-major-axis) radii
(Munoz+2018 centre/ellipticity/PA, D = 83.9 kpc). Fornax: Walker et al. (2009) MMFS
(`J/AJ/137/3100`), per-star sample rebuilt by inverse-variance averaging the per-exposure
table, using the Mg spectral index W' as the chemical discriminant, D = 147 kpc.

**Frameworks.**

- *Spherical Jeans* -- gNFW halo (alpha=1, eta=3), Osipkov-Merritt anisotropy, fit to the
  binned sigma_los(R) of the two metallicity halves.
- *GravSphere* -- from-scratch implementation of Read & Steger (2017): spherical Jeans with
  two Virial Shape Parameters and a free Baes-van Hese anisotropy profile beta(r),
  cross-checked against the reference code (sigma_los to 0.32%, VSPs to 0.14%).
- *Walker & Penarrubia (2011)* -- two chemo-dynamically distinct stellar subcomponents, each
  a Plummer sphere with Gaussian velocity and metallicity distributions; the estimator
  M(r_h) = 5 r_h sigma^2 / (2G) at two half-light radii gives Gamma = dlogM/dlogr, with
  Gamma > 2 excluding an NFW cusp. Sculptor is fitted members-only in nine parameters;
  Fornax retains non-members and adds a fitted member fraction plus a reconstructed
  selection function, giving ten.
- *Continuous f(J,[Fe/H])* -- a single DoublePowerLaw DF whose scale action varies smoothly
  with metallicity, `log10 J0(z) = logJ0 + k_J*(z - <z>)`; each star's likelihood
  marginalises over its true metallicity. `k_J < 0` places metal-rich stars on more tightly
  bound, less radially excursive orbits. **Adapted from Sanders & Binney (2015).**

---

## Scope & caveats

- Reported slopes are **MCMC posteriors, not point estimates**, because the sigma_los-only
  likelihood is degenerate in gamma. The prior sensitivity of that choice is measured, not
  assumed.
- The continuous method is an **adaptation of extended DFs (Sanders & Binney 2015)**, not a
  new construction; its contribution is `k_J` and the dSph application.
- **The continuous real-data chain has not converged.** The action-space gradient it encodes
  is corroborated independently by a direct computation on 574 stars with valid actions,
  which does not depend on chain convergence.
- All four methods assume **spherical symmetry** while using elliptical projected radii.
  Sculptor's intrinsic 3D shape is unobservable from one viewing angle, so the geometric
  systematic of Genina et al. (2018) cannot be tested here.
- Massari et al. (2018) measure Sculptor's anisotropy directly from HST+Gaia and prefer
  strongly radial orbits, which they note favours a cusp. Their 15-star sample gives a broad
  posterior that does not formally exclude this work's near-isotropic inner beta, but the
  direction of the preference runs against the shallow-slope result.
- For Fornax, the selection function is **reconstructed** from a Gaia DR3 photometric parent
  rather than the unpublished target catalogue, and the WP11 perspective term is omitted
  (quantified at 0.006 in Gamma, negligible). The selection reconstruction dominates,
  bounded at ~0.2 in Gamma by a colour-cut variant.
- WP11's absolute masses differ from their Table 4 by a constant estimator offset that
  **cancels in Gamma**.
- The AP25 curve overplotted in the rho(r) figure is computed from their published
  parameters via their Eq. 4, not by re-running their method.
- Gaia PM figures are illustrative; the core measurements use the spectroscopic sample.

---

## Repository layout

- `continuous_agama.py` -- the complete pipeline: all four frameworks, the mock
  experiments, the cross-validation against the reference GravSphere, and every figure.
- `requirements.txt` -- pinned dependencies.
- `README.md` -- this file.

---

*Developed as part of a research project on dwarf-spheroidal dark matter.*
