# h92 — single-type Hopenhayn (1992) model

A lean, single-type **Hopenhayn (1992)** stationary firm-dynamics model,
obtained from the explicit-`p`/`w` recode in [`../code`](../code) (the HNS22
benchmark) by stripping it down to the classic setup. It keeps the explicit
output price `p` and wage `w`, and the "fix `p = 1`, solve for `w` given a
hardcoded entry cost" default.

## What changed relative to `../code`

| Change | Detail |
|---|---|
| **Single firm type** | Only the "low" type is kept — no low/high mixture, no `omega`, one transition matrix `F`, one operating cost `cf`. |
| **No labor-force growth** | `g` and every `1+g` are removed. The stationary distribution is `mu = (I − Fnew)⁻¹ · (mstar·Gnew)`. |
| **Costs in goods** | The entry cost `ce` and operating cost `cf` are denominated in **goods** (output units), so their dollar costs are `p·ce` and `p·cf` — not `w·ce`, `w·cf`. |
| **No transition** | `feed_data.m` (the labor-force-series transition simulation) is removed; the model is stationary only. |
| **No lifecycle** | The age / survival / concentration-by-age block is removed from `stationary.m`. |

## Model

Firm's static problem (dollars), with `cf` in goods:

```
max_n  p·exp(s)·n^alpha − w·n − p·cf
FOC:   n*(s) = ( alpha·p·exp(s) / w )^(1/(1−alpha))                 [n.m]
Profit: Pi(s) = p·exp(s)·n*^alpha − w·n* − p·cf                     [prof_fn.m]
```

Incumbent value (dollars) and exit:

```
V(s) = max{ 0, Pi(s) + beta·E[V(s')|s] }          exit where V(s) = 0   [vfn.m]
```

Free entry, with `ce` in goods (dollar entry cost `p·ce`):

```
sum_s V(s;p,w)·G(s) = p·ce                          [entry_residual.m]
```

Since `cf` and `ce` scale with `p`, firm value is homogeneous of degree 1 in
`(p,w)`, so free entry pins only the ratio `z = p/w`. The default fixes
`p = 1` and solves for `w` ([`solve_wage.m`](solve_wage.m)); `z = p/w` is
reported as a diagnostic.

Stationary firm-size distribution (no growth), employment, and aggregates are
computed in [`stationary.m`](stationary.m). With `cf` and `ce` in goods, the
fixed cost and entry cost employ no labor, so firm employment is the
production labor `n*(s)` and total labor is `N = Σ n*(s)·mu(s)`.

## Files

```
h92/
├── main.m            DEFAULT: fix p=1, hardcode ce & cf (goods), solve for w
├── stationary.m      stationary equilibrium (single type, no growth, no lifecycle)
├── solve_wage.m      solve free entry E[V(s0)] = p·ce for w
├── entry_residual.m  free-entry residual E[V(s0)] − p·ce
├── vfn.m             value function iteration (dollars)
├── prof_fn.m         dollar profit p·eˢ·n*^α − w·n* − p·cf   (cf in goods)
├── n.m               labor demand n*(s) from the FOC
├── tauchen.m         single-AR(1) Tauchen discretization
├── check_h92.m       verification (free entry clears + homogeneity)
└── README.md
```

No data files are needed (the transition machinery that read them is removed).

## How to run

```matlab
% From the h92/ directory:
main         % solve and print the equilibrium; saves h92_benchmark.mat
check_h92    % verification suite
```

`main` prints:

```
=== Hopenhayn (1992): equilibrium (p fixed, solve for w) ===
  p (output price, numeraire) = 1
  w (wage, solved)            = 0.936762657168836
  z = p/w (relative price)    = 1.06750625929335
  ce (hardcoded, goods)       = 14.9087
  cf (hardcoded, goods)       = 15.1537
  E[V(s0)] at solved w        = 14.9087  (free entry: should equal p*ce)
  entry-condition residual    = 1.776e-15
  --
  exit-threshold productivity = 1.29025 (grid index 60)
  total employment N          = 201.136
  average firm size           = 114.012
  average entrant size        = 45.5648
  startup rate                = 0.1149
  exit rate                   = 0.1149
```

## Calibration (shared with h93)

The parameters hardcoded at the top of [`main.m`](main.m) are the
**Hopenhayn & Rogerson (1993)** calibration, taken from [`../h93`](../h93), so
this model and the HR1993 model use the *same* parameterization and grid:
`beta = 0.8`, `theta = alpha = 0.64`, `rho = 0.93`,
`sigma_eps = (1-theta)*sqrt(0.53) = 0.262`, `a = 0.061` (mean log s = 0.87),
`cf = 15.15`, `ce = 14.91` (goods), a fixed log-`s` grid with the top state at
`n* = 5000`, and entrants uniform over the bottom 74% of the grid.

The only difference between this folder and `../h93` is the **model**, not the
numbers: `h92` uses the exit-AFTER-shock timing `V(s) = max(0, prof + beta*E[V])`
(a firm exits the moment its value would go negative), while `h93` uses HR's
exit-BEFORE-shock timing. On the shared calibration this yields, e.g., average
firm size ≈ 114 and exit rate ≈ 0.11 here, versus ≈ 62 and ≈ 0.40 in `../h93`
— a clean illustration of how much the exit-timing assumption matters.

## Figures

[`makefigs.m`](makefigs.m) generates result plots (saved to `figures/`):
`fig1_size_distribution.png` (firm/employment size distribution vs HR data),
`fig2_value_policy.png` (value function `V(s)`, labor demand `n*(s)`, and the
exit threshold), and `fig3_cohort_survival.png` (cohort survival and exit
hazard by age).

## Verification

[`check_h92.m`](check_h92.m) confirms (all pass): free entry clears
(`E[V(s0)] = p·ce` to machine precision), the exit threshold is interior,
entry equals exit in steady state, and **homogeneity** — rescaling `p` leaves
`z`, the exit threshold, `N`, firm sizes, and the entry/exit rate invariant.

## Value-function iteration cap

As in `../code`, `vfn.m` caps the iteration at `params.maxiter` (default
`100`). Raise it (e.g. `5000`) for a fully converged solve when running
counterfactuals.

## Reference

Hopenhayn, H. (1992). "Entry, Exit, and Firm Dynamics in Long Run
Equilibrium." *Econometrica* 60(5), 1127–1150.
