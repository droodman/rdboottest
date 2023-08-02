{smcl}
{* *! rdboottest 0.1.0 1 August 2023}{...}
{help rdboottest:rdboottest}
{hline}{...}

{title:Title}

{pstd}
Supplements the {stata ssc describe rdrobust:rdrobust} command with wild bootstrap-based bias correction and confidence intervals using the method of 
He and Bartalotti (2020).

{title:Overview}
{pstd}
{cmd:rdboottest} wraps {stata ssc describe rdrobust:rdrobust}. Just change "rdrobust" to "rdboottest" in your command lines and add any optional
options to customize the wild bootstrapping. {cmd:rdboottest} will run {cmd:rdrobust}, show the output, and append bootstrap-based a bias-corrected estimate, p value,
and confidence interval. All results are saved in {cmd:e()} return macros.

{title:Syntax}

{phang}
{cmd:rdboottest} {it:depvar} {it:runvar} {ifin} [, {it:rdrobust_options} {it:rdboottest_options}]

{phang}
where {it:depvar}, {it:runvar}, and {it:rdrobust_options} follow the syntax defined in the {help rdrobust:rdobust help file} and
{it:rdboottest_options} are:

{synoptset 37 tabbed}{...}
{synopthdr:rdboottest_options}
{synoptline}
{synopt:{opt nobc}}suppress bias correction{p_end}
{synopt:{opt jack:knife} or {opt jk}}request jackknifing of residuals in bootstrap data construction{p_end}
{synopt:{opt rep:s(#)}}specifies number of replications for bootstrap estimation of confidence intervals; deafult is 999{p_end}
{synopt:{opt bcrep:s(#)}}specifies number of replications for estimation of bias; deafult is 500{p_end}
{synopt:{cmdab:weight:type(}{it:rademacher} {cmd:|} {it:mammen} {cmd:|} }specify weight type for bootstrapping; default is {it:rademacher}{p_end}
{synopt:{space 12} {it:webb} {cmd:|} {it:normal} {cmd:|} {it:gamma}{cmd:)}}{p_end}
{synopt:{cmdab:p:type(}{it:symmetric} {cmd:|} {it:equaltail} {cmd:|}}set p value type; {it:symmetric} is default{p_end}
{synopt:{space 12} {it:lower} {cmd:|} {it:upper}{cmd:)}}{p_end}
{synopt:{opt seed(#)}}initialize random number seed to {it:#}{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
Regression discontinuity design (RDD) offers hope of high-credibility causal identification in non-experimental settings. However, an intrinsic complication is that the assumptions needed
for "clean" identification hold only in the limit as one restricts to observational units closest to the identifying threshold. Real observational units
are not infinitessimally close, so all RDD estimates in principle contain bias. Worse, the bias depends in an unknown way on the kernel and bandwdith chosen
to operationalize nearness to the threshold.

{pstd}
This complication has given rise to a literature debating many proposals for kernel and bandwidth selection, bias estimation, and correction of variance estimates (e.g., Imbens and
Kalyanamaran 2012; Calonico, Cattaneo, and Tituniuk (CCT) 2014; Kolesar and Rothe 2018; 
Gelman and Imbens 2019; Pei et al. 2020; Kettlewell and Siminksi 2022 ({stata ssc describe pzms:pzms})). Bartalotti, Calhoun, and He (2017), along with He and
Bartalotti (2020), departs from this literature somewhat in exploiting the bootstrap to estimate both bias and variance. (I will refer to the latter paper
since it generalizes from the first.) Bootstrapping is appealing because it reduces reliance on asymptotic
theory to justify various complex bias and variance estimators. And the wild bootstrap is attractive in its ability to preserve traits of the original data, such as observation- and cluster-
specific variance structurs. The method of He and Bartalotti does not of course resolve the fundamental challenge of RDD noted above. In fact, as proposed and as implemented here, it
relies on the non-bootstrap CCT method for bandwidth selection.

{pstd}
The following description of the He and Bartalotti (2020) adaption of the wild bootstrap to RDD assumes familiarity with the application to linear models. (See 
Cameron, Gelbach, and Miller (2008) or Roodman et al. (2019). Stata implementations include the user-written {stata ssc describe boottest:boottest} and the built-in
{help wildbootstrap:wildbootstrap}.) From that starting point, He and Bartalotti introduce several modifications:

{p 4 7 0}1. It is a double bootstrap. A collection of simulated data sets in generated, and the researcher's estimator is applied to each, in order to empirically
estimate the variance. But before results such as {it:p} values and confidence intervals are extracted, each of these bootstrap estimates is bias-corrected
through a second layer of "wild" simulation (Horowitz 2001, p. 3172). For example, in the first iteration of the upper, "variance-simulation" level, a 
researcher's estimator could return 1.24. 1000 wild-bootstrap data sets are then constructed in which 1.24 is the true parameter value. If the researcher's
estimator returns 1.26 on average in these data sets, that points to a bias of +0.02, which is then subtracted from the original 1.24, yielding 1.22. The
run time of the double bootstrap is therefore proprortional to the product of the number of iterations at the variance-simulation and bias-correction levels. 

{p 4 7 0}2. Separately, and in general, as a paremetric bootstrap, the wild bootstrap is itself a two-stage process. An initial 
(data-generating process or DGP) regression is run as part of the process for constructing simulated data sets to which the researcher's specfication is applied
(replication regressions). In the usual wild bootstrap for linear models, the DGP regression specification is the same as the replication specification, or is restricted
by imposoition of the researcher's null hypothesis. In the He and Bartalotti adaption, roles are reversed. The DGP specification is that used by CCT to estimate bias,
typically with quadratic-polynomial controls in the running variable on each of the threshold and the CCT bandwidth. The researcher's null is not imposed. The replication
regressions are still necessarily the same as the researchers, which typically will mean use of linear controls and a distinct bandwidth.

{p 4 7 0}3. The object of the bootstrap is not a pivotal statistic such as a {it:t}, but the coefficient of interest. In other words, the method is bootstrap-c rather than
bootstrap-t.

{p 4 7 0}4. To reduce the influence of outliers, residuals from DGP regressions are jackknifed--individually or cluster-wise, as appropriate. (Jackknifing has been
added to {cmd:boottest} since Roodman et al. (2019); see MacKinnon, Nielsen, and Webb (2023))

{p 4 7 0}5. The method of Calonico, Cattaneo, and Titiiuk (2014) and Calonico, Cattaneo, and Farrell (2022) is used to select bandwidths. These
bandwidths, along with the researcher's kernel shape (triangular, Epanechnikov, rectangular), are incorporated into the DGP and replication regressions of the wild
bootstrap.

{pstd}
The method can be applied to sharp and fuzzy RDD. But since fuzzy RDD is an exactly-identified instrumental variables estimator, it
has no first moment. Its bias does not formally exist. It is ratio of the impact of crossing the threshold on the outcome to the impact on the treatment variable
and if there is substantial distributional mass for the latter near zero, the ratio is unstable. As a practical matter, 
bias-correction can be destabilizing when the discontinuity in treatment variable is not statistically clear.

{pstd}
Because the algorithm uses CCT bandwidths, and because {cmd:rdrobust} does not save all aspects of a user's specification in e() results, {cmd:rdboottest}
is wrapper for {cmd:rdrobust} rather than a post-estimation command. It starts by calling {cmd:rdrobust} and adds bootstrap results to the output
and the e() results.

{pstd}
{cmd:rdboottest} will adapt if the user overrides the {cmd:rdrobust} defaults for kernel type, bandwidth setting, or polynomial order of controls in the bias-correction
and estimation regressions. It also works when the left and right bandwidths differ. And it works for kink designs, requested through {cmd:rdrobust}'s {cmd:deriv()},
although this functionality has not been tested with simulations.


{marker options}{...}
{title:Options}

{phang}{opt nobc} prevents bias correction, which may be warranted in fuzzy RDD with a weak first stage.

{phang}{opt jack:knife} or {opt jk} requests jackknifing of the bootstrap data-generating process for OLS or linear IV estimates. At the DGP regression stage, the 
fit for each observation, or within each cluster, is computed using only the data from all
the other observations or clusters. He and Bartalotti (2020)'s simulation results incorporate jackknifing.

{phang}{opt rep:s(#)} sets the number of replications at the variance-sampling stage. The default is 999. Especially when clusters are few, increasing this number costs little in run 
time.

{phang}{opt bcrep:s(#)} sets the number of replications at the bias-correction stage. The default is 500. Since run time is proportional to the product of the {cmdab:rep:s(#)}
and {cmdab:bcrep:s(#)} options, and can be high in some cases, the question arises of which option to set higher. A general principle is that acturately estimating second moments
through simulation takes more replications than first moments. This argues in the direction of higher {cmdab:rep:s(#)} and lower {cmdab:bcrep:s(#)}

{phang}{opt p:type(symmetric | equaltail | lower | upper)} sets the p value type. The default, {it:symmetric}, has the p value derived from the
square of the {it:t}/{it:z} statistic, or, equivalently, the absolute value. {it:equaltail} performs a two-tailed test using the {it:t}/{it:z} statistic. For example, 
if the confidence level is 95, then the symmetric p value is less than 0.05 if the square of the test statistic is in the top 5 centiles of the corresponding bootstrapped 
distribution. The equal-tail p value is less than 0.05 if the test statistic is in the top or bottom 2.5 centiles. In addition, {it:lower} and {it:upper} allow
one-sided tests.

{phang}{opt weight:type(rademacher | mammen | webb | normal | gamma)} specifies the type of random weights to apply to the residuals or scores from the base regression, when
bootstrapping. The default is {it:rademacher}. However if the number of replications exceeds 2^(# of clusters)--that is, the number of possible
Rademacher draws--{cmd:boottest} will take each possible draw once. It will not do that with Mammen weights even though the same issue arises. In all such cases,
Webb weights are probably better.

{phang}{opt seed(#)} sets the initial state of the random number generator. See {help set seed}.


{title:Stored results}

{pstd}
{cmd:rdboottest} adds the following to {cmd:rdrobust}'s {cmd:e()} results:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(reps)}}number of bootstrap replications for variance simulation{p_end}
{synopt:{cmd:e(bcreps)}}number of bootstrap replications for bias correction{p_end}
{synopt:{cmd:e(p_bs)}}bootstrap p value that the parameter of interest is 0{p_end}
{synopt:{cmd:e(bc_wb)}}0/1 indicating whether bias correction was performed in bootstrapping{p_end}
{synopt:{cmd:e(jk_wb)}}0/1 indicating whether bootstrap used jackknifing{p_end}
{synopt:{cmd:e(tau_bc_wb)}}bootstrap-based bias-corrected parameter estimate{p_end}
{synopt:{cmd:e(ci_l_rb_wb)}}lower bound of bootstrap CI{p_end}
{synopt:{cmd:e(ci_r_rb_wb)}}upper bound of bootstrap CI{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(seed)}}value of {opt seed(#)} option, if any, or else c(seed) at command invocation{p_end}
{synopt:{cmd:e(weighttype_wb)}}bootstrapping weight type{p_end}
{synopt:{cmd:e(ptype_wb)}}bootstrapping p value type{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(dist)}}bootstrap distribution of the (bias-corrected) coefficient estimate{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}sample marker, which {cmd:rdrobust} omits{p_end}


{title:Author}

{p 4}David Roodman{p_end}
{p 4}david@davidroodman.com{p_end}


{title:Example}

{phang}. {stata `"use https://github.com/rdpackages/rdrobust/raw/master/stata/rdrobust_senate"'}{p_end}

{phang}. {stata rdboottest vote margin}{p_end}


{title:References}
{p 4 8 2}Bartalotti, Otávio, Gray Calhoun, and Yang He. 2017. Bootstrap Confidence Intervals for Sharp Regression Discontinuity Designs*. In Regression Discontinuity Designs, 38:421–53. Advances in Econometrics. Emerald Publishing Limited. https://doi.org/10.1108/S0731-905320170000038018.{p_end}
{p 4 8 2}Calonico, Sebastian, M.D. Cattaneo, and M.H. Farrell. 2022. Coverage Error Optimal Confidence Intervals for Local Polynomial Regression. Bernoulli 28 (4) 2998-3022, November. https://doi.org/10.3150/21-BEJ1445.{p_end}
{p 4 8 2}Calonico, Sebastian, Matias D. Cattaneo, and Rocio Titiunik. 2014. Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs. Econometrica: Journal of the Econometric Society 82 (6): 2295–2326. https://doi.org/10.3982/ecta11757.{p_end}
{p 4 8 2}Cameron, A.C., J.B. Gelbach, and D.L. Miller. 2008. Bootstrap-based improvements for inference with clustered errors.{p_end}
{p 4 8 2}Gelman, Andrew, and Guido Imbens. 2019. Why High-Order Polynomials Should Not Be Used in Regression Discontinuity Designs. Journal of Business & Economic Statistics: A Publication of the American Statistical Association 37 (3): 447–56. https://doi.org/10.1080/07350015.2017.1366909.{p_end}
{p 4 8 2}He, Y., and O. Bartalotti. 2020. Wild Bootstrap for Fuzzy Regression Discontinuity Designs: Obtaining Robust Bias-Corrected Confidence Intervals. The Econometrics Journal 23 (2): 211–31. https://doi.org/10.1093/ectj/utaa002.{p_end}
{p 4 8 2}Horowitz, J.L. 2001. Chapter 52 - The Bootstrap. In Handbook of Econometrics, edited by James J. Heckman and Edward Leamer, 5:3159–3228. Elsevier. https://doi.org/10.1016/S1573-4412(01)05005-X.{p_end}
{p 4 8 2}Imbens, Guido, and Karthik Kalyanaraman. 2012. Optimal Bandwidth Choice for the Regression Discontinuity Estimator. The Review of Economic Studies 79 (3): 933–59. https://doi.org/10.1093/restud/rdr043.{p_end}
{p 4 8 2}Kettlewell, N. and P. Siminksi. 2021. Optimal model selection in RDD and related settings using placebo zones. https://www.uts.edu.au/sites/default/files/2021-10/Optimal%20Model%20Selection%20in%20RDD%20Kettlewell%20N%20Siminski%20P_1.pdf{p_end}
{p 4 8 2}Kolesár, Michal, and Christoph Rothe. 2018. Inference in Regression Discontinuity Designs with a Discrete Running Variable. The American Economic Review 108 (8): 2277–2304. https://doi.org/10.1257/aer.20160945.{p_end}
{p 4 8 2}MacKinnon, J.G., M.O. Nielsen, and M.D. Webb. 2023. Fast and Reliable Jackknife and Bootstrap Methods for Cluster-Robust Inference. Queen's Economics Department Working Paper No. 1485.{p_end}
{p 4 8 2}Pei, Zhuan, David S. Lee, David Card, and Andrea Weber. 2022. Local Polynomial Order in Regression Discontinuity Designs. Journal of Business & Economic Statistics: A Publication of the American Statistical Association 40 (3): 1259–67. https://doi.org/10.1080/07350015.2021.1920961.{p_end}
{p 4 8 2}Roodman, D., J. MacKinnon, M. Nielsen, and M. Webb. 2019. Fast and wild: bootstrap inference in Stata using boottest. {it:Stata Journal} 19(1): 4-60.{p_end}
