# Wandering bumps in stochastic neural fields

MATLAB code that reproduces every figure of the paper: the three analytic bump results and the seven simulation figures.

Spatially localized bumps of persistent activity in a neural field can sit anywhere on a ring when the synaptic weights are translation invariant, so they are only marginally stable to a shift. In the model, firing is a threshold function of the synaptic input, the weight is a cosine kernel, and additive noise with a cosine spatial correlation drives the field. Noise then makes the bump wander as pure diffusion, and the code compares the measured diffusion coefficient with the small-noise prediction. It also covers bump extinction near the saddle-node, pinning of the bump by a stationary input or by weak periodic modulation of the weights, switching between attractors, and the reduced effective diffusion that a high-frequency modulation of the weights gives.

## Citation

Zachary P. Kilpatrick and Bard Ermentrout. Wandering bumps in stochastic neural fields. SIAM Journal on Applied Dynamical Systems 12(1):61–94 (2013). [doi:10.1137/120877106](https://doi.org/10.1137/120877106)

```bibtex
@article{kilpatrick2013wandering,
  author  = {Kilpatrick, Zachary P. and Ermentrout, Bard},
  title   = {Wandering Bumps in Stochastic Neural Fields},
  journal = {SIAM Journal on Applied Dynamical Systems},
  volume  = {12},
  number  = {1},
  pages   = {61--94},
  year    = {2013},
  doi     = {10.1137/120877106}
}
```

## Requirements

MATLAB R2021a or later. No toolboxes are needed. The figures are written with `exportgraphics`.

## Quick start

Clone the repository, `wandering-bumps`, from https://github.com/zpkilpat/wandering-bumps. Then start MATLAB in the repository root.

```matlab
run_all                      % all ten figures with the grid solver (paper scale)
run_all('fast', true)        % exact 2D reduction wherever there is no input: minutes, not hours
figure5                      % a single figure
figure5('quick', true)       % fewer realizations and shorter runs, a quick smoke test
```

From a shell:

```bash
matlab -batch "run_all('fast',true)"
```

`figureN` draws Fig. N of the paper and writes it to `output/figs/figN.png` and `output/figs/figN.pdf`. Options (`'quick'`, `'fast'`, `'M'`, `'seed'`, `'outdir'`) are name-value pairs, described in `src/plotting/wb_opts.m`. Figures 1 to 4 are analytic and take seconds.

## Figures

| file | paper figure | contents |
| --- | --- | --- |
| `figure1.m` | 1 | input-locked bumps: half-width against input strength, stable and unstable profiles |
| `figure2.m` | 2 | the 2n bump locations under synaptic modulation, and which of them are stable |
| `figure3.m` | 3 | bumps with modulation n = 1: half-width against heterogeneity, profiles |
| `figure4.m` | 4 | eigenvalue of shift perturbations against heterogeneity and against modulation frequency |
| `figure5.m` | 5 | free diffusion: single realizations, variance against time, D against threshold |
| `figure6.m` | 6 | extinction near the saddle-node: single run, mean extinction time and its exponential fit |
| `figure7.m` | 7 | bumps pinned by a stationary input: saturating variance for n = 1 and 2 |
| `figure8.m` | 8 | switching between the two attractors of a bimodal input, and the mean waiting time |
| `figure9.m` | 9 | pinning by low-frequency synaptic modulation, n = 2 and 3 |
| `figure10.m` | 10 | reduced effective diffusion under high-frequency modulation |

## Solvers and runtime

The weight and the noise only ever drive the cos x and sin x modes, so the integral in the model is evaluated exactly in O(N) per realization. `nf_sim` is the grid solver the paper describes: Euler–Maruyama with Δx ≈ 0.01, Δt = 0.01, trapezoid rule, and 1000 realizations for the variance estimates.

`nf_sim_fast` is an exact reduction for the cases with no input. The field is then U = a(t) cos x + b(t) sin x, the active set is a single arc, and the integrals have closed forms, so each realization is a two-variable SDE with the same Euler step. With the same seed it reproduces the grid solver's statistics to three digits, about 40 times faster. `fast_vs_grid_check` in `src/checks/` shows the comparison. Figures 7 and 8 have an input and always use the grid solver. The default is the grid solver, which is the paper's method.

| | runtime |
| --- | --- |
| `run_all` with the grid solver | about 1.5 hours |
| `run_all('fast', true)` | about 10 minutes, mostly Figs. 7 and 8 |
| `run_all('quick', true)` | a few minutes |

Times are for MATLAB R2023b on a laptop.

## What the paper does not state

* **Fig. 7 input strength.** I0 = 0.1, backed out of the saturation levels printed in the figure.
* **Fig. 8 switch criterion.** A switch is recorded when the bump centre first reaches |x| ≥ 0.9π. The Fig. 8b runs use Δt = 0.05 because they are very long, and a smaller number of realizations than the 1000 used elsewhere. Panel (a) shows a run with a hop, chosen by scanning seeds.
* **Fig. 6a.** The single run is chosen by scanning seeds for one that goes extinct near t = 65.
* **Fig. 6b and 8b fits** are least squares on the mean event time, so expect scatter in the fitted exponent.
* **Fig. 10d simulation points** use the long-time slope of the variance over the second half of a T = 200 run.
* The initial condition is the stable bump A cos x, relaxed deterministically for 30 time units so that it sits on the attractor of the modulated or input-driven network. The bump centre is the centre of the active set, unwrapped in time.

## Two formulas in the paper that do not match its figures

1. **Eq. (2.19), the shift eigenvalue for n = 1.** The closed form as printed, −σa / (2 sin a + σa + σ sin 2a), is slightly off. Evaluating (2.18) directly gives σ(sin a cos a − a)/A with A = 2 sin a + σ(a + sin a cos a), and a deterministic simulation of a shifted bump decays at exactly that rate (−0.0997 against −0.1001 at σ = 0.2). The paper's Fig. 4 follows the printed form, so `figure4` draws the printed curve (solid) and the direct evaluation (dashed and ×). Everywhere else the code uses (2.18), which for n > 1 agrees with the closed form (2.20) to machine precision.
2. **Eq. (B.9), the potential amplitude.** The printed expression equals −2κ/n², with κ = −λ_o from (4.5). The cosine well that follows from the SDE (3.26) with K(nΔ) = (κ/n) sin(nΔ) has amplitude κ/n², and that reproduces the paper's Fig. 10d curve. The code uses κ/n²; `wb_theory().Veff` returns both.

## Repository layout

```
.
├── run_all.m                   regenerate every figure
├── figure1.m ... figure10.m    one file per figure
├── src/
│   ├── model/
│   │   ├── nf_sim.m            grid Euler–Maruyama solver, vectorized over realizations
│   │   ├── nf_sim_fast.m       exact two-variable reduction (no input)
│   │   └── wb_theory.m         every analytic formula: bump existence, stability, D, kappa, V, D_eff
│   ├── plotting/               options, panel labels, field plots, export
│   └── checks/                 fast_vs_grid_check
└── output/figs/                figN.png and figN.pdf from a full run
```

## License

MIT. See `LICENSE`.
