# hopenhayn_neira_singhania_2022 — HNS22 recoded with explicit `p` and `w`

A recoding of the benchmark firm-dynamics model of **Hopenhayn, Neira and
Singhania (2022)** ("Firm Dynamics and the Declining Startup Rate", HNS22) in
which the **output price `p` and the wage `w` appear explicitly**, instead of
the single composite price `z = p/w` used in the original replication code.
This folder is self-contained: the model, the data (`data_summary_stats/`) and
the verification script (`check_benchmark.m`) all live here.

## What changed and why

In the HNS22 firm problem a firm with productivity state `s` chooses labor `n`:

```
max_n   p·exp(s)·n^alpha  −  w·n  −  w·cf
```

The first-order condition is derived explicitly in [`n.m`](n.m) and
[`prof_fn.m`](prof_fn.m):

```
FOC:   alpha·p·exp(s)·n^(alpha−1) − w = 0
   =>  n*(s) = ( alpha · p · exp(s) / w )^(1/(1−alpha))          [n.m]

Profit (dollars), at n*:
   Pi(s) = p·exp(s)·n*^alpha − w·n* − w·cf                       [prof_fn.m]
         = revenue − wage bill − fixed-cost bill
```

Both prices appear **explicitly** in their natural roles, not folded into a
single ratio. Because the firm's problem is homogeneous, real decisions still
depend only on `z = p/w`; the original code exploited this and carried `z`
(named `zstar`), fixing `z = 1` and *computing* the entry cost `ce`. This
recode carries `p` and `w` separately and reports `z = p/w` as a diagnostic.

| | Original HNS22 | This recode |
|---|---|---|
| Price object | single `z = zstar` | explicit `p` **and** `w` |
| Numeraire | `z = 1` fixed | `p = 1` fixed |
| `ce` | **output** (computed at `z=1`) | **input** (hardcoded) |
| Solved for | — | **`w`** (free entry) |

## Default experiment ([`main.m`](main.m))

Fix `p = 1`, hardcode the entry cost `ce`, and solve free entry for the wage
`w` ([`solve_wage.m`](solve_wage.m)). With `ce` in units of labor the dollar
entry cost is `w·ce`, so ([`entry_residual.m`](entry_residual.m)):

```
omega_h · Σ_s V_h(s;p,w)·G(s) + omega_l · Σ_s V_l(s;p,w)·G(s) = w·ce
```

The hardcoded `CE = 0.012815402758375272` is the entry cost the original code
produced at `z = 1`, so with `p = 1` the solver returns **`w = 1`** and the
recode reproduces the published benchmark exactly. Change `CE` to study how the
equilibrium responds to a different entry cost.

## Units and the `maxiter` cap

Profit and the value function are in **dollars** with `p, w` explicit; every
physical quantity (employment, distributions, total labor `N`, the labor
share) is invariant to the overall `(p,w)` scale, so fixing `p = 1` is a pure
numeraire choice. `vfn.m` caps value-function iteration at `params.maxiter`
(default `100`, the HNS22 value, which does not reach `tol`); it is kept so the
recode reproduces the *calibrated* benchmark. Raise it (e.g. `5000`) for a
fully converged solve.

## File mapping from the original

- `n.m`, `prof_fn.m`: signature `(...,z)` → `(..., p, w)`, rewritten from the
  first-order conditions with `p, w` explicit (dollar profit).
- `vfn.m`: `vfn(F, zstar, cf, params)` → `vfn(F, p, w, cf, params)`; dollar
  value `V`; iteration cap `params.maxiter`.
- `stationary.m`: reads `p, w` and takes `ce` as an **input**; reports `p, w, z`.
- `solve_wage.m` + `entry_residual.m`: the wage analogue of the original
  `zstar_fun.m` + `fzero` (solve for `w` given `ce`).
- `stationary_alt.m`: the old "solve for `zstar`" now fixes `p` and solves `w`.
- `feed_data.m`: only the labor-share block changed (`p, w` explicit).
- `tauchen2.m`, `makefigs.m`, `maketables.m`: unchanged.

## How to run

```matlab
% From this folder:
main            % default experiment (p=1, solve for w)
check_benchmark % verification (regression + homogeneity, all pass)
```

`check_benchmark` confirms the recode reproduces the original z-based
benchmark to relative error `0`, and that rescaling `p` leaves every physical
moment invariant.

`maketables.m` reproduces the paper's Table 6, which needs alternative-
experiment output files not shipped with the original; Tables 3–5 and `makefigs.m`
work from the `benchmark.mat` that `main.m` saves.

## Reference

Hopenhayn, H., Neira, J., and Singhania, R. (2022). "Firm Dynamics and the
Declining Startup Rate." *Econometrica*.
