cd D:\OneDrive\Documents\Macros\rdboottest

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

  cap noi rdboottest Y X, fuzzy(T) bwselect(cerrd) all `vceopt' // occassionally crashes...
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

set seed 1231
sim, g(10) ζ(.04) ρ(.9) jk(1) bc(0) μ(cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X))

parallel init 6

cap program drop sims
program define sims
  syntax anything(name=TableName), [ρs(string) gs(string) reps(integer 100)]

  cap postclose `TableName'
  postfile `TableName' double ζ ρ G DGP ζhatCL ζhatRBC ζhatWBS ECCL ECRBC ECWBS ILCL ILRBC ILWBS bCEO hCEO SDCL SDRBC SDWBS using `TableName', every(1) replace
  foreach g in `gs' {
    foreach ρ in `ρs' {
      forvalues DGP=1/3 {
        local ζ: word `DGP' of .04 -3.45 .04
        local μ: word `DGP' of "cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X)" ///
                               "cond(X<0, (2.3+(3.28+(1.45+(.23+.03*X)*X)*X)*X)*X, (18.49+(-54.81+(74.3+(-45.0+9.83*X)*X)*X)*X)*X)" ///
                               "cond(X<0, (1.27+(3.59+(14.147+(23.694+10.995*X)*X)*X)*X)*X, (.84+(-.3+(2.397+(-.901+3.56*X)*X)*X)*X)*X)"
        parallel sim, reps(`reps') nodots proc(2) exp(ζhatCL=r(ζhatCL) ILCL=r(ILCL) ECCL=r(ECCL) ζhatRBC=r(ζhatRBC) ILRBC=r(ILRBC) ECRBC=r(ECRBC) ζhatWBS=r(ζhatWBS) ILWBS=r(ILWBS) ECWBS=r(ECWBS) bCEO=r(bCEO) hCEO=r(hCEO)): ///
          sim, g(`g') ζ(`ζ') ρ(`ρ') μ(`μ')
        collapse ζhat* EC* IL* *CEO (sd) SDCL=ζhatCL SDRBC=ζhatRBC SDWBS=ζhatWBS
        post `TableName' (`ζ') (`ρ') (`g') (`DGP') (ζhatCL) (ζhatRBC) (ζhatWBS) (ECCL) (ECRBC) (ECWBS) (ILCL) (ILRBC) (ILWBS) (bCEO) (hCEO) (SDCL) (SDRBC) (SDWBS)
      }
    }
  }
  postclose `TableName'
  use `TableName', clear
  foreach estimator in CL RBC WBS {
    gen bias`estimator' = ζhat`estimator' - ζ
    gen RMSE`estimator' = sqrt(bias`estimator'^2 + SD`estimator'^2)
    drop ζhat`estimator'
  }
end

sims HBTable1, ρs(0 -.9 .9) gs(1000)
reshape long bias SD RMSE EC IL, i(ρ DGP *CEO) j(estimator) string
list ρ DGP estimator bias SD RMSE EC IL *CEO, sep(0)

sims HBTable2, ρs(0) gs(5 10 25)
reshape long bias SD RMSE EC IL, i(G DGP *CEO) j(estimator) string
list G DGP estimator bias SD RMSE EC IL *CEO, sep(0)

