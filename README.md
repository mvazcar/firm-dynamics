# firm-dynamics

A recoding of the benchmark firm-dynamics model of **Hopenhayn, Neira and
Singhania (2022)** ("Firm Dynamics and the Declining Startup Rate", HNS22) in
which the **output price `p` and the wage `w` appear explicitly**, instead of
the single composite price `z = p/w` used in the original replication code.

The original code is preserved unchanged under [`hns22_raw/`](hns22_raw/) for
reference. The recoded model lives in [`code/`](code/).

## What changed and why

In the HNS22 firm problem a firm with productivity state `s` chooses labor `n`:

```
max_n   p·exp(s)·n^alpha  −  w·n  −  w·cf
```

The first-order condition is derived explicitly in
[`code/n.m`](code/n.m) and [`code/prof_fn.m`](code/prof_fn.m):

```
FOC:   alpha·p·exp(s)·n^(alpha−1) − w = 0
   =>  n*(s) = ( alpha · p · exp(s) / w )^(1/(1−alpha))          [n.m]

Profit (dollars), at n*:
   Pi(s) = p·exp(s)·n*^alpha − w·n* − w·cf                       [prof_fn.m]
         = revenue − wage bill − fixed-cost bill
```

Both the output price `p` (on revenue) and the wage `w` (on the labor and
fixed-cost bills) appear **explicitly** in their natural roles — not folded
into a single ratio. Because the firm's problem is homogeneous, its real
decisions still depend on prices only through the ratio `z = p/w` (the real,
wage-deflated output price); the original code exploited this and carried the
single variable `z` (named `zstar`), fixing `z = 1` and *computing* the entry
cost `ce` that makes `z = 1` an equilibrium. This recode instead carries `p`
and `w` as separate objects and reports `z = p/w` as a derived diagnostic.

This recode keeps `p` and `w` as **separate, explicit objects** everywhere, and
inverts the default experiment:

| | Original (`hns22_raw`) | Recode (`code`) |
|---|---|---|
| Price object | single `z = zstar` | explicit `p` **and** `w` |
| Numeraire | `z = 1` fixed | `p = 1` fixed |
| `ce` | **output** (computed at `z=1`) | **input** (hardcoded) |
| Solved for | — | **`w`** (free entry) |

### The default experiment

[`code/main.m`](code/main.m):

1. Fixes the output price `p = 1` (numeraire).
2. Hardcodes the entry cost `ce` (in units of labor) — constant `CE` at the top
   of `main.m`.
3. Solves the free-entry condition for the wage `w`
   ([`code/solve_wage.m`](code/solve_wage.m)).

Free entry requires the expected **dollar** value of an entrant to equal the
**dollar** entry cost. The entry cost is `ce` units of labor, i.e. `w·ce`
dollars, so ([`code/entry_residual.m`](code/entry_residual.m)):

```
omega_h · Σ_s V_h(s; p,w)·G(s)  +  omega_l · Σ_s V_l(s; p,w)·G(s)  =  w·ce
```

Because firm value `V` is homogeneous of degree 1 in `(p,w)` and `ce` is in
labor units, this condition pins down only the **ratio** `z = p/w`. Fixing
`p = 1` turns "solve for the entry-clearing price" into "**solve for the wage
`w`**", and `z = p/w` is recovered afterwards.

The hardcoded `CE = 0.012815402758375272` is the entry cost the original code
produced at `z = 1`. With `p = 1` the solver therefore returns **`w = 1`**
(hence `z = 1`), and the recode reproduces the published benchmark exactly (see
Verification). Change `CE` in `main.m` to study how the equilibrium wage and
firm dynamics respond to a different entry cost.

## Units convention

Profit and the value function are in **dollars**, with `p` and `w` explicit:
`prof_fn` returns `Pi(s) = p·exp(s)·n*^alpha − w·n* − w·cf` and `vfn` returns the
dollar value `V`. Free entry is `E[V(s0)] = w·ce`. Because the firm problem is
homogeneous, the real allocation depends on prices only through `z = p/w`, so
every physical quantity — employment, firm-size and age distributions, total
labor `N`, the aggregate labor share — is **invariant to the overall `(p,w)`
scale**; only the ratio matters, which is why fixing `p = 1` is a pure
numeraire choice. (`N` and the entry cost `ce` are denominated in units of
labor; the aggregate labor share is the dollar wage bill `w·n` over dollar
revenue `p·exp(s)·n^alpha`, in which the wage cancels.)

## Value-function iteration cap

`vfn.m` caps the iteration at `params.maxiter`, **default 100** — the value
HNS22 uses. At 100 iterations the VFI does **not** reach `tol = 1e-8` (the
operator is a `beta`-contraction needing several hundred iterations), so the
solved value is a mild approximation. This is preserved as the default so the
recode reproduces the **calibrated** benchmark exactly (the published `ce` and
moments were themselves produced with the 100-iteration operator). For a fully
converged solve — recommended when running counterfactuals with a different
`CE` — set `params.maxiter` higher (e.g. `5000`); note this departs slightly
from the published calibration.

## Repository layout

```
firm-dynamics/
├── code/                      recoded model (explicit p, w)
│   ├── main.m                 DEFAULT: fix p=1, hardcode ce, solve for w
│   ├── solve_wage.m           solve free entry for w given ce and p
│   ├── entry_residual.m       free-entry residual E[V(s0;p,w)] − w·ce
│   ├── stationary.m           stationary equilibrium given (p, w)
│   ├── stationary_alt.m       alt-experiments variant (also solves for w)
│   ├── vfn.m                  value function iteration (dollars), args (F, p, w, cf)
│   ├── prof_fn.m              dollar profit p·eˢ·n*^α − w·n* − w·cf, args (params, s, p, w)
│   ├── n.m                    labor demand n*(s) from the FOC, args (params, s, p, w)
│   ├── feed_data.m            feeds labor-force series through the model
│   ├── tauchen2.m             Tauchen discretization (unchanged)
│   ├── makefigs.m             figures (unchanged; loads benchmark.mat)
│   └── maketables.m           tables (unchanged; see note below)
├── data_summary_stats/        data series (copied from hns22_raw)
├── tests/
│   └── check_benchmark.m      regression + homogeneity verification
└── hns22_raw/                 original HNS22 code, untouched (reference)
```

### File-by-file mapping from the original

- `n.m`, `prof_fn.m`: signature `(...,z)` → `(..., p, w)`, rewritten from the
  **first-order conditions** with `p` and `w` explicit. `n.m` is the FOC labor
  demand `(alpha·p·exp(s)/w)^(1/(1−alpha))`; `prof_fn.m` is dollar profit
  `p·exp(s)·n*^alpha − w·n* − w·cf` (was the labor-unit `Pi/w`).
- `vfn.m`: `vfn(F, zstar, cf, params)` → `vfn(F, p, w, cf, params)`; returns the
  **dollar** value `V`; the iteration cap is now `params.maxiter`.
- `stationary.m`: reads `p, w` instead of `zstar`; takes `ce` as an **input**
  rather than computing it; reports `p`, `w`, `z` and an entry-condition
  diagnostic `entry_value_check` (= `w·ce` at the solved wage); warns on a
  degenerate no-exit equilibrium.
- `solve_wage.m` + `entry_residual.m`: **new**. Together they replace the old
  `zstar_fun.m` and do the "solve for `w` given `ce`" step, with a post-solve
  check that free entry actually clears. (`zstar_fun.m` is not carried over.)
- `stationary_alt.m`: the old "solve for `zstar` given `ce`" now fixes `p` and
  solves for `w` via `solve_wage`.
- `feed_data.m`: only the labor-share block changed — `rev_firms` is written
  with `p` and `w` explicit (kept in labor units so `w` cancels in the share).
- `tauchen2.m`, `makefigs.m`, `maketables.m`: unchanged (no price variable).

## How to run

Requires MATLAB (developed on R2024b; uses `normcdf` from the Statistics and
Machine Learning Toolbox and base `fzero`).

```matlab
% From the code/ directory — default experiment (p=1, solve for w):
main

% From the tests/ directory — verification:
check_benchmark
```

`main.m` prints the solved prices and saves `benchmark.mat`, e.g.:

```
=== Equilibrium prices (default: p fixed, solve for w) ===
  p (output price, numeraire) = 1
  w (wage, solved)            = 1
  z = p/w (relative price)    = 1
  ce (hardcoded, labor units) = 0.0128154027583753
  E[V(s0)] at solved w        = 0.0128154027583753  (free entry: should equal w*ce)
  w*ce (dollar entry cost)    = 0.0128154027583753
  entry-condition residual    = 1.041e-17
```

## Verification

[`tests/check_benchmark.m`](tests/check_benchmark.m) runs two suites (all pass):

1. **Benchmark regression.** With `p = 1` and the hardcoded `ce`, the solved
   wage is `w = 1`, free entry clears (`E[V(s0)] − w·ce ≈ 1e-17`), and 25
   equilibrium and transition moments (`ce`, `N`, exit thresholds, startup/exit
   rates, average size, concentration, survival, labor share, ...) match the
   original z-based code to relative error `0`.
2. **Homogeneity in `(p,w)`.** Rescaling `p` by `lambda = 2` and re-solving
   yields `w = 2` (so `z = 1` unchanged) and leaves every physical moment,
   including the aggregate labor share, invariant (to solver tolerance) —
   confirming the recode never silently introduces a wage-level dependence.

## Note on `maketables.m`

`maketables.m` reproduces the paper's Table 6, which requires the alternative-
experiments output files (`alteqbm_*`). The driver that produces those files is
not part of `hns22_raw`, so Table 6 will not run out of the box; Tables 3–5 (and
`makefigs.m`) work from the `benchmark.mat` that `main.m` saves.

## Reference

Hopenhayn, H., Neira, J., and Singhania, R. (2022). "Firm Dynamics and the
Declining Startup Rate." *Econometrica*.
