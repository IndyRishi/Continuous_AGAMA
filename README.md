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

This repository accompanies the manuscript *Constraining the Inner Dark-Matter Slope of
Sculptor: A Comparative Analysis of Dynamical Methods*. Every number below matches the
submitted version; where the two could drift, the manuscript is authoritative.

---

## Key results

Real Sculptor data, Tolstoy et al. (2023), 1339 stars with reliable [Fe/H].

| Framework | Method | Inner slope |
|---|---|---|
| Spherical Jeans (`--dm5`) | gNFW, sigma_los only | gamma = 0.97 (+0.21 / -0.33) |
| GravSphere (`--gravsphere`) | Jeans + Virial Shape Parameters + free beta(r) | gamma = 0.46 (+0.40 / -0.28) |
| Walker & Penarrubia 2011 (`--wp11`) | Two-subcomponent mass slope | Gamma = 2.49 (+0.23 / -0.20) -> gamma = 0.50; 99.6% of the posterior above Gamma = 2 |
| Continuous f(J,[Fe/H]) (`--continuous`) | Single DF, metallicity as coordinate | non-converging; action-space gradient detected independently |

GravSphere and the two-population estimator agree on a **shallow (core-like) asymptotic**
slope, consistent with the published gamma = 0.39 (+0.23 / -0.26) of Arroyo-Polonio et al.
(2025) and the gamma = 0.5 +/- 0.3 of Zhu et al. (2016) from independent data. The spherical
Jeans fit sits higher at 0.97; a prior-sensitivity test shows this is **not** a prior
artifact (below), so the disagreement between methods is real rather than an inference
failure.

### The central-density cross-check

The same posteriors, evaluated on the `rho_DM(150 pc)` diagnostic of Read, Walker & Steger
(2019) rather than on the asymptotic slope:

```bash
python continuous_agama.py --rho150
```

| Source | rho_DM(150 pc)  [1e8 Msun/kpc^3] | local slope at 150 pc |
|---|---|---|
| GravSphere (this work) | 1.77 (+0.33 / -0.26) | -1.00 (+0.19 / -0.23) |
| Spherical Jeans (this work) | 1.70 +/- 0.11 | -1.11 (+0.16 / -0.13) |
| Read, Walker & Steger (2019) | 1.49 (+0.28 / -0.23) | -0.83 (+0.30 / -0.25) |

Agreement at **0.74 sigma** against an independent GravSphere implementation on a different
photometric catalog. Every posterior sample of both chains lies above their 1e8 threshold,
placing Sculptor in their **cusp-like** class.

This does not contradict the shallow gamma: the two are different quantities. `gamma` is the
asymptotic slope as r -> 0, while the slope *attained* at 150 pc is consistent with NFW. The
two fits here agree on the density to 4% while their asymptotic slopes differ by 0.51 --
the gamma-r_s-M_DM degeneracy runs almost exactly along the direction that preserves the
density near the constrained radius. The core-like gamma is an extrapolation inward of the
innermost kinematic constraint.

### What the mocks actually show

The two-population split is often assumed to bias the slope. Controlled experiments
separate the two candidate causes:

- **Null control** (one Plummer tracer, no metallicity structure): split and all-star fits
  agree within their standard errors (+0.174 vs +0.162 at gamma_true = 0.4).
- **Label-scramble control** (two Plummers, metallicity permuted): likewise
  (-0.082 vs -0.085).
- Introducing the **two-Plummer geometry alone** moves the recovery by -0.26 at
  gamma_true = 0.4 and -0.44 at 1.0; restoring the true metallicity link moves it back by
  only +0.05 and +0.13.

So **77-83% of the apparent decomposition bias is tracer geometry, not the decomposition** --
which means mock-based claims about population splitting must control for geometry, this
work's own included.

Which model recovers best depends on which is right: under continuous truth the split
biases high (+0.459, +0.296) while the continuous fit does better (+0.165, -0.009); under
**discrete** truth the ordering reverses (split -0.029, -0.292; continuous -0.094, -0.420).
For Sculptor, which independent photometry, chemodynamics and orbit-based modeling indicate
is discrete, the split is the better-matched model.

### Negative results worth recording

- **Marginal metallicity unimodality does not settle anything.** Sculptor's [Fe/H] is
  unimodal by three tests (BC = 0.30, dBIC = -8, Hartigan dip p = 0.79) -- but two-population
  mocks matched to the literature decomposition return the same (dip p ~ 0.96), because the
  components overlap at 1.35 sigma. The test cannot discriminate a continuum from two
  overlapping populations.
- **A two-component proper-motion mixture is not identifiable** on a spectroscopically
  preselected sample. The fit returns P_mem = 1 for every star (sd = 8e-13); bounding the
  foreground width relocates the pathology rather than removing it. Membership is therefore
  characterised by direct offset from the systemic PM.
- **The sigma_los-only likelihood is weakly constraining in gamma** (dchi2 < 1 for
  gamma <~ 0.4, dchi2 < 4 for gamma <~ 0.9, with the minimum at or below the lower scan
  limit), which is why slopes are reported as posterior medians. Refitting under a prior
  with median 0.75 rather than 0.95 moves the Jeans posterior by **0.02** -- so the
  posterior is set by the data, not the prior, despite the flat likelihood.

### Fornax as a cross-galaxy validation

Applying the same implementation to the Walker et al. (2009) MMFS Fornax sample recovers
**Gamma = 3.16 (+0.83 / -0.63)** against the published 2.61 (+0.43 / -0.37) -- agreement at
0.8 sigma -- with f_mem, f_sub, the chemical separation and both velocity dispersions
agreeing within their uncertainties.

Reaching that required rebuilding the per-star sample from the per-exposure table (the
averaged Mg-index column is populated for only 459 of 2633 rows, against the 2603 stars
WP11 modelled), retaining non-members with a fitted member fraction, and reconstructing the
spectroscopic completeness w(R) against a Gaia DR3 photometric parent. The residual
difference in half-light radii (~20% larger than theirs) traces to that reconstruction: the
target catalogue behind their w(R) is unpublished. A magnitude-only parent, which admits
main-sequence foreground, gives Gamma = 2.97 (+0.71 / -0.56); the difference bounds the
selection systematic at ~0.2 in Gamma.

Note that the Fornax posterior sits largely above Gamma = 3, where the local-power-law
assumption behind `Gamma = 3 - gamma` fails. It is reported as a constraint on Gamma and on
the exclusion of an NFW cusp, and **no gamma is inferred from it**.

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

VizieR mirrors occasionally fall out of sync and return an empty result rather than an
error. Each loader prints its row count (1339 of 1701 for Sculptor, 2603 rebuilt from 3150
Fornax exposures); if those disagree, set another mirror through the `VIZIER_SERVER`
environment variable.

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

### Derived quantities

```bash
python continuous_agama.py --rho150         # rho_DM(150 pc) and local slope vs Read+2019
python continuous_agama.py --menc           # enclosed mass at the fitted half-light radii
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

### Recovering lost outputs

Chains survive in their HDF5 backends even when the derived `.npy` or figure is lost or
overwritten. Both can be rebuilt without re-sampling:

```bash
python continuous_agama.py --rebuild-npy dm5.h5      # flattened chain from its backend
python continuous_agama.py --cont-corner contfull.h5 # continuous corner plot, seconds not hours
```

### Common flags

`--steps N` (MCMC steps to add), `--walkers N`, `--nproc N` (parallel processes), `--K N`
(metallicity DF nodes for `--continuous`), `--nsub N` (fit a random subsample; omit for the
full sample), `--backend FILE.h5` (checkpoint), `--no-resume`, `--galaxy {sculptor,fornax}`,
`--mock`, `--seed N` (walker initialisation), `--mockseed N` (mock draw), `--gammaprior LO HI`
(Jeans prior bounds, for the prior-sensitivity test), `--nbins N` (sigma_los bins per
population, default 6), `--ra KPC` (Osipkov-Merritt anisotropy radius), `--astar KPC`
(GravSphere tracer scale), `--vclip SIGMA` (velocity outlier clip),
`--selection {none,radial,2d}` (AP24-style completeness for the continuous DF).

**Sensitivity runs are isolated by construction.** Any flag that changes the likelihood
redirects the run to its own backend, `.npy` and figures -- `--dm5 --nbins 8` writes
`dm5_nb8.h5`, `dm5_nb8_chain.npy` and `figure_dm5_corner_nb8.png`, leaving the headline
files untouched. `_backend_config_guard()` additionally refuses to resume a chain whose
stored configuration differs from the current one, so a variant can never silently append
to a headline chain.

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

HDF5 backends named by `--backend`, plus `.npy` posterior samples. The archived deposit
contains the following, mapped to where each is used in the manuscript:

| Backend | Configuration | Used in |
|---|---|---|
| `dm5.h5` | headline spherical Jeans, 6 bins, r_a = 1.5 kpc, seed 42 | Table 9, Fig. E3 |
| `dm5_seed1.h5`, `dm5_seed7.h5` | alternative walker initialisations | Table E2 (seed stability) |
| `dm5_gp0-1.5.h5` | narrower gamma prior, median 0.75 | Table E3 (prior sensitivity) |
| `dm5_ra100.h5` | r_a = 100 kpc, effectively isotropic | Table E4 (anisotropy systematic) |
| `dm5_nb8.h5` | 8 sigma_los bins per population | Sec. 2.2 (binning check) |
| `gravsphere.h5` | headline GravSphere, 7 bins | Table 9, Figs. 11, E2 |
| `wp11.h5` | headline two-population, Sculptor | Table 9, Figs. 12, E1 |
| `fornax_wp11.h5` | Fornax cross-galaxy validation | Sec. 3.9, Figs. 15, E5 |
| `contfull.h5` | continuous DF, full 1339-star sample | Table 5, Fig. E4, App. E |
| `cont500_good.h5` | seeded restart, 500-star subsample, 56 walkers | App. E (boundary diagnosis) |

The two continuous-DF chains are **non-converging by design of the model family, not by
under-sampling** -- see Appendix E of the manuscript. They are archived because the
diagnosis rests on them: the walker split at ln P = -6753 vs -6695, and the pile-up of
`g_z`, `h_z` against their 1.5 prior ceiling.

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
| `figure_dm5_corner.png` | Spherical Jeans posterior |
| `figure_gravsphere_corner.png` | GravSphere posterior |
| `figure_wp11_corner.png` | Two-population posterior (Sculptor) |
| `figure_continuous_corner.png` | Continuous f(J,[Fe/H]) posterior |
| `figure_action_space.png` | Action-space chemodynamics ([Fe/H] vs J_r) |
| `figure_dispersion_profile.png` | Radial velocity-dispersion profile (LOS + PM) |
| `figure_gaia_skymap.png` | Gaia PM sky map with membership-cut panels |
| `figure_fig4_all_chains.png` | DM density rho(r) from all frameworks + AP25's curve |
| `fornax_figure_wp11.png` | Fornax two-population fit (cross-galaxy validation) |
| `fornax_figure_wp11_corner.png` | Fornax two-population posterior |

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
  cross-checked against the reference code (sigma_los to 0.32%, theory VSPs to 0.10%,
  observed-VSP estimators to 1 part in 1e12).
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
- **The continuous real-data chain does not converge.** The action-space gradient it encodes
  is corroborated independently by a direct computation on 575 stars with valid actions,
  which does not depend on chain convergence.
- All four methods assume **spherical symmetry** while using elliptical projected radii.
  Sculptor's intrinsic 3D shape is unobservable from one viewing angle, so the geometric
  systematic of Genina et al. (2018) cannot be fully excluded. The *misalignment mechanism*
  specifically is testable and is tested: the two subcomponents' major axes agree to
  -4.1 +/- 4.3 deg, placing Sculptor with the simulated galaxies for which the cusp is
  correctly recovered rather than with the misaligned case that produces a spurious core.
- Massari et al. (2018) measure Sculptor's anisotropy directly from HST+Gaia and prefer
  strongly radial orbits. Strigari et al. (2018) reanalyse the same measurements and find
  cored and cuspy potentials predict almost identical transverse dispersions at that radius,
  so they do not discriminate; two orbit-based analyses (Breddels et al. 2013; Breddels &
  Helmi 2013) independently support the near-isotropic-to-tangential inner beta found here.
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
- `figs/` -- generated figures.
- `README.md` -- this file.

---

## Citation

If you use this code, please cite the manuscript and this archive. Large-language-model
assistance was used in developing and debugging the pipeline; every quantitative result is
produced by code independently validated as described in Appendix B of the manuscript and
in the cross-check commands above.

Licensed CC-BY.
