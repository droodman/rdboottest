cd D:\OneDrive\Documents\Macros\rdboottest

set rmsg on
set more off

cap program drop sim
program define sim, rclass
  syntax, μ(string asis) [n(integer 1000) c(real .9) ρ(real 0) ζ(real .04) g(integer 1000) bc(integer 1) jk(integer 1) B1(integer 500) B2(integer 999)]

  drop _all
  set obs `n'
  mat C = 1, `ρ' \ `ρ', 1
  drawnorm double ut uy, corr(C)
  replace uy = .1295 * uy
  gen double X = 2 * rbeta(2,4) - 1
  gen byte T = ut <= invnormal(cond(X<0, .5-`c'/2, .5+`c'/2))
  gen double Y = `ζ'*T + uy + `μ'
  
  if `g' < `n' {
    gen int clustid = floor((_n-1)/`g') + 1
    local vceopt vce(cluster clustid)
    local clustidopt st_data(.,"clustid")
  }
  else local clustidopt J(0,1,0)

  noi rdboottest Y X, fuzzy(T) bwselect(cerrd) all `vceopt' jk  // occassionally rdrobust crashes...
  if _rc exit

  return scalar ζhatCL = _b[Conventional]
  return scalar ζhatRBC = _b[Robust]
  return scalar ILCL  = (e(ci_r_cl) - e(ci_l_cl)) / 2
  return scalar ILRBC = (e(ci_r_rb) - e(ci_l_rb)) / 2
  return scalar hCEO = e(h_l)
  return scalar bCEO = e(b_l)
  return scalar ζhatWBS = e(tau_bc_wb)
  return scalar ILWBS = e(ci_r_rb_wb) - e(ci_l_rb_wb)
  
  foreach estimator in CL RBC WBS {
    return scalar EC`estimator' = `ζ' > return(ζhat`estimator') - return(IL`estimator')/2 & `ζ' < return(ζhat`estimator') + return(IL`estimator')/2
  }
  ereturn clear
end

set seed 204375
parallel init 6  // use 6 Stata copies

cap program drop sims
program define sims
  syntax anything(name=TableName), [ρs(string) gs(string) reps(integer 5000)]

  cap postclose `TableName'
  postfile `TableName' double ζ ρ G DGP ζhatCL ζhatRBC ζhatWBS ECCL ECRBC ECWBS ILCL ILRBC ILWBS bCEO hCEO SDCL SDRBC SDWBS using sim\\`TableName', every(1) replace
  foreach g in `gs' {
    foreach ρ in `ρs' {
      forvalues DGP=1/3 {
        local ζ: word `DGP' of .04 -3.45 .04
        local μ: word `DGP' of "cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X)" ///
                               "cond(X<0, (2.3+(3.28+(1.45+(.23+.03*X)*X)*X)*X)*X, (18.49+(-54.81+(74.3+(-45.0+9.83*X)*X)*X)*X)*X)" ///
                               "cond(X<0, (1.27+(3.59+(14.147+(23.694+10.995*X)*X)*X)*X)*X, (.84+(-.3+(2.397+(-.901+3.56*X)*X)*X)*X)*X)"
        parallel sim, exp(ζhatCL=r(ζhatCL) ILCL=r(ILCL) ECCL=r(ECCL) ζhatRBC=r(ζhatRBC) ILRBC=r(ILRBC) ECRBC=r(ECRBC) ζhatWBS=r(ζhatWBS) ILWBS=r(ILWBS) ECWBS=r(ECWBS) bCEO=r(bCEO) hCEO=r(hCEO)) proc(2) reps(`reps') nodots : ///
          sim, g(`g') ζ(`ζ') ρ(`ρ') μ(`μ')
        collapse ζhat* EC* IL* *CEO (sd) SDCL=ζhatCL SDRBC=ζhatRBC SDWBS=ζhatWBS
        post `TableName' (`ζ') (`ρ') (`g') (`DGP') (ζhatCL) (ζhatRBC) (ζhatWBS) (ECCL) (ECRBC) (ECWBS) (ILCL) (ILRBC) (ILWBS) (bCEO) (hCEO) (SDCL) (SDRBC) (SDWBS)
      }
    }
  }
  postclose `TableName'
end

sims HBTable1, reps(1000) ρs(0 -.9 .9) gs(1000)
sims HBTable2, reps(1000) ρs(0) gs(5 10 25)

use sim\HBTable1, clear
qui foreach estimator in CL RBC WBS {
  gen bias`estimator' = ζhat`estimator' - ζ
  gen RMSE`estimator' = sqrt(bias`estimator'^2 + SD`estimator'^2)
  drop ζhat`estimator'
}
qui reshape long bias SD RMSE EC IL, i(ρ DGP *CEO) j(estimator) string
list ρ DGP estimator bias SD RMSE EC IL *CEO, sep(0)

use sim\HBTable2, clear
qui foreach estimator in CL RBC WBS {
  gen bias`estimator' = ζhat`estimator' - ζ
  gen RMSE`estimator' = sqrt(bias`estimator'^2 + SD`estimator'^2)
  drop ζhat`estimator'
}
qui reshape long bias SD RMSE EC IL, i(G DGP *CEO) j(estimator) string
list G DGP estimator bias SD RMSE EC IL *CEO, sep(0)

