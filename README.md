# rdboottest
Stata package for bias correction and robust variance estimation in regression discontinuity design using the wild bootstrap ([He and Bartalotti 2020](https://doi.org/10.1093/ectj/utaa002)).

This Stata command wraps [rdrobust](https://github.com/rdpackages/rdrobust). To use it, replace `rdrobust` in your command lines with
`rdboottest`. It will run `rdrobust`, show results, and append bootstrap-based estimates of the bias-corrected coefficent, p value,
and confidence interval. It will add these results to e(), along with the sample marker and other results normally missing from `rdrobust`
return values.

## Installation
When more mature, this package will be posted on SSC. For now, install it in Stata with
```
net install rdboottest, replace from(https://raw.github.com/droodman/rdboottest/vX.Y.Z)
```
where `X.Y.Z` represents the [latest release version number](https://github.com/droodman/rdboottest/releases).

## Documentation
Install and type `help rdboottest` in Stata.

## Example
As in He and Bartalotti (2020), this example [jackknifes residuals used in the bootstrap](https://onlinelibrary.wiley.com/doi/full/10.1002/jae.2969) and constructs an _equal-tailed_ 95% confidence interval, in which 2.5% of the distribution lies beyond either end. It defaults to Rademacher weights.
```
. use https://github.com/rdpackages/rdrobust/raw/master/stata/rdrobust_senate

. rdboottest vote margin, seed(71438) jk ptype(equaltail)

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
    Wild bootstrap |  7.5839                       0.000    4.13822      11.2495
--------------------------------------------------------------------------------
Bias-corrected. Bootstrap method of He and Bartalotti (2020)
```

# Simulation evidence
Emulating Tables 1 and 2 of He and Bartalotti (2020), this table shows the bias, standard deviation, root-mean-squared error,
and empirical coverage (size) of three estimators. The three are conventional (CL) and robust bias-corrected (RBC), both as produced
by `rdrobust` with the `bwselect(cerrd)` option; and the wild bootstrap, using the same bandwidths.


### Fuzzy RDD, non-clustered data, 1000 observations
```
  +-------------------------------------------------------+
  |   ρ   DGP   estimator    bias      SD    RMSE      EC |
  |-------------------------------------------------------|
  | -.9     1          CL   0.015   0.050   0.052   0.934 |
  | -.9     1         RBC   0.011   0.054   0.055   0.943 |
  | -.9     1         WBS   0.013   0.053   0.055   0.939 |
  | -.9     2          CL   0.065   0.074   0.098   0.659 |
  | -.9     2         RBC   0.017   0.066   0.068   0.904 |
  | -.9     2         WBS   0.019   0.064   0.067   0.955 |
  | -.9     3          CL   0.002   0.052   0.052   0.948 |
  | -.9     3         RBC   0.002   0.057   0.057   0.947 |
  | -.9     3         WBS   0.004   0.055   0.056   0.954 |
  |   0     1          CL   0.017   0.051   0.053   0.920 |
  |   0     1         RBC   0.014   0.054   0.056   0.931 |
  |   0     1         WBS   0.014   0.054   0.055   0.939 |
  |   0     2          CL   0.072   0.071   0.101   0.697 |
  |   0     2         RBC   0.024   0.063   0.067   0.926 |
  |   0     2         WBS   0.023   0.062   0.066   0.969 |
  |   0     3          CL   0.006   0.051   0.051   0.940 |
  |   0     3         RBC   0.007   0.055   0.056   0.940 |
  |   0     3         WBS   0.007   0.055   0.055   0.953 |
  |  .9     1          CL   0.018   0.051   0.055   0.945 |
  |  .9     1         RBC   0.015   0.055   0.058   0.955 |
  |  .9     1         WBS   0.013   0.054   0.055   0.949 |
  |  .9     2          CL   0.067   0.066   0.094   0.804 |
  |  .9     2         RBC   0.021   0.063   0.066   0.952 |
  |  .9     2         WBS   0.018   0.061   0.063   0.949 |
  |  .9     3          CL   0.006   0.051   0.051   0.947 |
  |  .9     3         RBC   0.007   0.055   0.056   0.950 |
  |  .9     3         WBS   0.005   0.054   0.054   0.939 |
  +-------------------------------------------------------+
```
Notes: ρ = correlation between first and second-stage errors. DGP is a data-generating process defined in He and Bartalotti (2020).
SD = standard deviation. RMSE = root-mean-square error. EC = empirical coverage, which is ideally 0.95. All bootstrap estimates use a
triangular kernel, jackknifing, equal-tailed confidence intervals, and Rademacher weights. 1000 simulations are run for each case (each row).

### Fuzzy RDD, clustered data, 1000 observations

```
  +------------------------------------------------------+
  |  G   DGP   estimator    bias      SD    RMSE      EC |
  |------------------------------------------------------|
  |  5     1          CL   0.016   0.048   0.050   0.927 |
  |  5     1         RBC   0.012   0.052   0.053   0.942 |
  |  5     1         WBS   0.012   0.051   0.053   0.935 |
  |  5     2          CL   0.077   0.070   0.104   0.662 |
  |  5     2         RBC   0.024   0.063   0.068   0.902 |
  |  5     2         WBS   0.022   0.063   0.066   0.937 |
  |  5     3          CL   0.009   0.051   0.052   0.938 |
  |  5     3         RBC   0.010   0.059   0.060   0.940 |
  |  5     3         WBS   0.011   0.069   0.070   0.932 |
  | 10     1          CL   0.019   0.049   0.053   0.893 |
  | 10     1         RBC   0.015   0.054   0.056   0.914 |
  | 10     1         WBS   0.015   0.053   0.055   0.901 |
  | 10     2          CL   0.081   0.065   0.104   0.628 |
  | 10     2         RBC   0.027   0.060   0.065   0.902 |
  | 10     2         WBS   0.025   0.059   0.065   0.951 |
  | 10     3          CL   0.006   0.047   0.048   0.936 |
  | 10     3         RBC   0.007   0.052   0.053   0.943 |
  | 10     3         WBS   0.007   0.052   0.052   0.932 |
  | 25     1          CL   0.020   0.047   0.051   0.913 |
  | 25     1         RBC   0.015   0.052   0.054   0.934 |
  | 25     1         WBS   0.015   0.051   0.054   0.918 |
  | 25     2          CL   0.084   0.069   0.109   0.613 |
  | 25     2         RBC   0.024   0.062   0.067   0.904 |
  | 25     2         WBS   0.023   0.061   0.066   0.917 |
  | 25     3          CL   0.006   0.049   0.050   0.931 |
  | 25     3         RBC   0.007   0.055   0.055   0.929 |
  | 25     3         WBS   0.007   0.055   0.055   0.914 |
  +------------------------------------------------------+
```
G = number of clusters. All simulations have ρ=0. Notes to previous table apply.
