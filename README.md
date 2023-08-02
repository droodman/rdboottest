# rdboottest
Stata package for bias correction and robust variance estimation in regression discontinuity design using the wild bootstrap ([He and Bartalotti 2020](https://doi.org/10.1093/ectj/utaa002)).

This Stata command wraps [rdrobust](https://github.com/rdpackages/rdrobust). To use it, replace `rdrobust` in your command lines with
`rdboottest`. It will run `rdrobust`, show results, and append bootstrap-based estimates of the bias-corrected coefficent, p value,
and confidence interval. It will also add these results to e().

## Installation
When more mature, this package will be posted on SSC. For now, install it from Stata with
```
net install rdboottest, replace from(https://raw.github.com/droodman/rdboottest/v[X.Y.Z])
```
where "[X.Y.Z]" represents the [latest release version number](https://github.com/droodman/rdboottest/releases).

## Documentation
Install and type `help rdboottest` in Stata.

## Example
```
. use https://github.com/rdpackages/rdrobust/raw/master/stata/rdrobust_senate

. rdboottest vote margin, seed(1438)

Sharp RD estimates using local polynomial regression.

      Cutoff c = 0 | Left of c  Right of c            Number of obs =       1297
-------------------+----------------------            BW type       =      mserd
     Number of obs |       595         702            Kernel        = Triangular
Eff. Number of obs |       360         323            VCE method    =         NN
    Order est. (p) |         1           1
    Order bias (q) |         2           2
       BW est. (h) |    17.754      17.754
       BW bias (b) |    28.028      28.028
         rho (h/b) |     0.633       0.633

Outcome: vote. Running variable: margin.
--------------------------------------------------------------------------------
            Method |   Coef.    Std. Err.    z     P>|z|    [95% Conf. Interval]
-------------------+------------------------------------------------------------
      Conventional |  7.4141     1.4587   5.0826   0.000     4.5551      10.2732
            Robust |     -          -     4.3110   0.000     4.0937      10.9193
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
    Wild bootstrap |  7.4383                       0.000    4.59704      10.2795
--------------------------------------------------------------------------------
Bias-corrected. Bootstrap method of He and Bartalotti (2020)
```
