# h93 — Hopenhayn & Rogerson (1993)

A replication of **Hopenhayn & Rogerson (1993, JPE)**, "Job Turnover and Policy
Evaluation: A General Equilibrium Analysis" — the Hopenhayn (1992) industry
model in general equilibrium, used to evaluate a tax on job destruction
(firing costs). The paper is [`h93.pdf`](h93.pdf).

The model has two regimes:

- **`τ = 0` benchmark** — with no firing cost the firm's problem is static in
  employment, so the state is just productivity `s` (**1-D**). This is where
  the model is calibrated to the data (Table 1/2).
- **`τ > 0` experiments** — a firing cost `τ·max(0, n₋₁−n)` makes employment
  sticky, so the state becomes `(s, n₋₁)` (**2-D**) and the policy is an
  (S,s)-style inaction band. These are the policy experiments (Table 3).

## Model

Firm state `(s, n₋₁)`, employment choice `n'`. Prices are explicit (`p`, `w`),
with `p` fixed as the numeraire and `w` solved from free entry. Bellman
(paper p.7), in dollars:

```
W(s, n₋₁) = max_{n'≥0} { p·s·n'^θ − w·n' − p·cf − w·τ·max(0, n₋₁−n')     ← return_fn
                         + β·max[ E_{s'|s} W(s', n'),  −w·τ·n' ] }
                                └ stay next period ┘  └ exit: fire all ┘
```

- **Exit is decided at the start of a period, before the new shock** (HR's
  timing): a firm with shock `s` and chosen `n'` exits next period iff
  `E_sW(s',n') < −w·τ·n'`. This is why the `max` sits on the *continuation*.
- Costs: `cf`, `ce` in **goods** (dollar cost `p·cf`, `p·ce`); the firing cost
  is in **wages** (`τ` = fraction of the period wage; `τ=0.2` ≈ one year's pay
  for a 5-year period).
- **Entrants** draw `s ~ v`, have `n₋₁ = 0`, and pay **no `cf`** (footnote 5).
  Free entry: `Σ_s v(s)·W_e(s) = p·ce`.
- **Household** (Hansen indivisible labor, `u=ln c`, `v=A·N`) pins the scale
  (mass of entrants) via labor-market clearing and delivers the employment
  and welfare numbers. `A` is set so employment `N = 0.6` at `τ = 0`.

Production is `s·n^θ` with `s` in **levels** (`s = exp(log s)`); the shock is
`log s' = a + ρ·log s + ε`.

## Calibration (paper values + calibrated to Table 1/2 targets)

Structural (from the paper): `β=0.8`, `θ=0.64`, `ρ=0.93`,
`σ_ε=(1−θ)√0.53=0.262`, 5-year period. Calibrated at `p=w=1` to the Table 1
targets (mean size 61.7, 5-yr exit 0.37, avg entrant size 7.5): `a=0.061`
(mean log s = 0.87), `cf=15.15`, `ce=14.91`, entrants uniform over the bottom
74% of the `s`-grid. A fixed log-`s` grid has its top state at `n*=5000`.
(This is the same calibration used by [`../h92`](../h92).)

## Results

**Table 2 (benchmark, `τ=0`)** — [`calibrate_benchmark.m`](calibrate_benchmark.m)
matches average size (61.7 vs 61.2), exit rate (0.40 vs 0.39), job turnover
(0.31 vs 0.30), growth variance (0.54 vs 0.55), and the firm/employment size
distributions; the one soft moment is the survivor serial correlation (see
Caveats).

**Table 3 (policy, `τ = 0, 0.1, 0.2`)** — [`hr1993_2d.m`](hr1993_2d.m), at the
`nz=100, na=250` grid:

| | τ=0 | τ=.1 | τ=.2 | HR (0/.1/.2) |
|---|---|---|---|---|
| Price `z=p/w` (rel.) | 1.000 | 1.022 | 1.039 | 1.000 / 1.026 / 1.048 |
| Consumption (output) | 100 | 97.9 | 96.2 | 100 / 97.5 / 95.4 |
| Average productivity | 100 | 99.1 | 97.9 | 100 / 99.2 / 97.9 |
| Total employment | 100 | 98.6 | 97.9 | 100 / 98.3 / 97.5 |
| Utility-adj. consumption | 100 | 98.9 | 97.8 | 100 / 98.7 / 97.2 |
| Average firm size | 61.7 | 62.8 | 64.2 | 61.2 / 61.8 / 65.1 |
| Layoff costs / wage bill | 0 | .025 | .043 | 0 / .026 / .044 |
| Job turnover rate | .308 | .253 | .216 | .30 / .26 / .22 |
| Serial corr. log n | .79 | .83 | .84 | .92 / .94 / .94 |
| Var. of growth | .54 | .43 | .37 | .55 / .45 / .39 |

The model reproduces HR's headline findings: a firing tax of one year's wages
(`τ=0.2`) **cuts job turnover ~30%**, **lowers average productivity ~2%**,
**reduces employment ~2%**, and **costs ~2.2% of consumption in welfare** — the
productivity/welfare channel that is the paper's main message. At `τ=0` the 2-D
model collapses exactly to the 1-D benchmark (`w=1`, `z=1`, avg size 61.7).

## Files

```
h93/
├── calibrate_benchmark.m   tau=0 benchmark: calibrate + report Table 2 (1-D)
├── hr1993_2d.m             tau>0 driver: solve, close with household, Table 3
├── return_fn.m             per-period return p·s·n'^θ − w·n' − p·cf − w·g
├── labor_adjustment.m      firing cost g(n',n₋₁) = τ·max(0, n₋₁−n')
├── vfn_2d.m                2-D value iteration (nested exit max)
├── stationary_2d.m         2-D stationary distribution + Table-3 aggregates
├── solve_wage_2d.m         free entry Σ v·W_e = p·ce, solved for w
├── household.m             representative-household closure (employment, welfare)
├── tauchen_fixed.m         Tauchen transition on a fixed (top = n*=5000) grid
├── README.md
└── h93.pdf                 the paper
```

## How to run

```matlab
% From the h93/ directory:
calibrate_benchmark   % tau=0 benchmark, reports Table 2
hr1993_2d             % tau = 0,.1,.2 policy experiments, reports Table 3
```

`hr1993_2d` defaults to `nz=100, na=250` (tight match; a few minutes); set
`P.nz=50, P.na=150` at the top for fast iteration.

## Figures

[`makefigs.m`](makefigs.m) generates the result plots (saved to `figures/`):

- `fig1_decision_bands.png` — the employment decision-rule bands
  `[n_l(s), n_u(s)]` for `τ = 0, 0.1, 0.2` (the visual form of HR Table 4),
  plus the fractional band width showing the distortion grow with `τ`.
- `fig2_table3_effects.png` — the Table-3 policy effects vs `τ`, model vs HR.
- `fig3_size_distribution.png` — firm and employment size distributions.
- `fig4_cohort_hazard.png` — cohort exit-hazard by age, model vs HR data.
- `fig5_stationary_measure.png` — the stationary measure over `(s, n₋₁)`.

## Caveats

- **Survivor serial correlation** (0.79–0.84 vs HR 0.92–0.94): the *level* is
  low because the exit threshold sits mid-distribution, so the surviving
  sample is truncated (the *unconditional* autocorrelation is 0.93, exactly ρ).
  It rises with `τ` as HR report. Closing the level would need HR's exact
  (under-specified) entrant distribution / exit age-profile. The calibration
  can be refined later.
- **Formulation vs the VFI Toolkit**: this uses HR's own grid construction and
  the paper's `σ_ε`; it is not bit-identical to other replications.

## Reference

Hopenhayn, H. and Rogerson, R. (1993). "Job Turnover and Policy Evaluation: A
General Equilibrium Analysis." *Journal of Political Economy* 101(5), 915–938.
