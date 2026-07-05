# h93 — Hopenhayn & Rogerson (1993)

The **Hopenhayn & Rogerson (1993, JPE)** firm-dynamics model at its baseline
firing tax `tau = 0`, at which it reduces to the frictionless single-type
Hopenhayn (1992) model. This folder reuses that model engine — identical to
[`../h92`](../h92) — and supplies the HR1993 calibration in [`main.m`](main.m).
The paper is included as [`h93.pdf`](h93.pdf).

Since `tau = 0` removes the labor-adjustment (firing) cost, the firm's state
is just its productivity (no lagged-employment state), and the model is
exactly the goods-cost, single-type, no-growth stationary Hopenhayn model of
`../h92`. Only the calibration differs.

## Calibration

Parameter values follow the [VFI Toolkit replication of
HR1993](https://github.com/vfitoolkit/vfitoolkit-matlab-replication/tree/master/HopenhaynRogerson1993):

| Parameter | Value |
|---|---|
| `beta` (discount factor) | `0.8` |
| `alpha` (curvature, HR's `theta`) | `0.64` |
| `cf` (fixed cost, goods) | `12` |
| `ce` (entry cost, goods) | `40` |
| log-AR(1) `log z' = a + rho·log z + eps` | `a = 0.078`, `rho = 0.93`, `sigma_eps = sqrt((1-rho)·0.53)` |
| Tauchen grid | `n_z = 20`, half-width `q = 4` |
| entrant distribution | uniform over the bottom 65% of productivity states |

Productivity enters production in levels (`z = exp(svec)`, `svec` the log-z
grid), so `prof_fn`'s `p·exp(s)·n^alpha` equals HR's `p·z·n^alpha`. Costs are
in goods (`p·cf`, `p·ce`), matching HR.

## Model files

Same engine as `../h92` (copied so the folder is self-contained):
`n.m`, `prof_fn.m`, `vfn.m`, `entry_residual.m`, `solve_wage.m`,
`stationary.m`, `tauchen.m`. The only folder-specific file is
[`main.m`](main.m) (the HR calibration) plus [`check_h93.m`](check_h93.m).

## How to run

```matlab
% From the h93/ directory:
main         % solve and print the HR1993 equilibrium; saves h93_benchmark.mat
check_h93    % verification (free entry clears + homogeneity)
```

`main` prints, for example:

```
=== Hopenhayn-Rogerson (1993) calibration, tau = 0 (p fixed, solve for w) ===
  p (output price, numeraire) = 1
  w (wage, solved)            = 0.483448724196798
  z = p/w (relative price)    = 2.06847169089421
  ce (hardcoded, goods)       = 40
  cf (hardcoded, goods)       = 12
  E[V(z0)] at solved w        = 40  (free entry: should equal p*ce)
  entry-condition residual    = 7.105e-15
  --
  exit-threshold productivity (log z) = 0.783317 (grid index 9 of 20)
  average firm size (employment)      = 235.302
  average entrant size                = 92.7481
  startup rate                        = 0.0300432
  exit rate                           = 0.0300432
```

## Normalization and comparison caveats

- **Normalization.** HR1993 fix the wage (`w = 1`) and let the output price
  clear free entry. This model instead fixes `p = 1` (numeraire) and solves
  free entry `E[V(z0)] = p·ce` for `w`. Only `z = p/w` matters, so the real
  allocation is identical; the reported `z = p/w` equals HR's equilibrium
  output price under their `w = 1` normalization.
- **Continuous vs grid employment.** This model uses the continuous
  first-order-condition labor demand `n*(s)`, whereas the VFI Toolkit
  discretizes employment on a grid capped at 5000. Upper-tail firm sizes (and
  hence the average) therefore differ from the Toolkit's numbers even at
  identical parameters.
- **Entrant mass.** `mstar` is a free scale normalization here (the Toolkit
  pins it via the representative household), so total employment `N` is not
  directly comparable; firm sizes and turnover rates are.

## Reference

Hopenhayn, H. and Rogerson, R. (1993). "Job Turnover and Policy Evaluation: A
General Equilibrium Analysis." *Journal of Political Economy* 101(5), 915–938.
