cd D:\OneDrive\Documents\Macros\rdboottest

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

class WBSRDD {
  real scalar fuzzy, N, WWd, WWr, coarse, N_G, B₁, B₂, jk
  real colvector W, ZWd, ZWr, Wd, WrKr, halfWrKr, X, Kd, Kr, retval, MZdWrKr, WdKd, M
  real matrix info, Zr, Zd, invZZd, invZZr, ZdinvZZd, Z, ZdKd
real colvector clustid
  real colvector bc(), vs()
  void Prep()
  real matrix fold()
}

real matrix WBSRDD::fold(matrix X) return(uppertriangle(X) + lowertriangle(X,0)')  // fold matrix diagonally; returns same values as a quad form, but runs faster because of all the 0's


// one-time stuff that only depends on exog vars and cluster and kernel definitions
void WBSRDD::Prep(real scalar B₁, real scalar B₂, real colvector clustid, real colvector X, real matrix Z, real colvector Kd, real colvector Kr, real scalar fuzzy, real scalar jk) {
  real matrix ZZd, ZZr

  this.B₁ = B₁
  this.B₂ = B₂
  this.fuzzy = fuzzy
  this.jk = jk
this.clustid = clustid
  this.X = X
  this.Z = Z
  this.Kd = Kd
  this.Kr = Kr

  retval = J(B₂+1, 1, 0)
  N = rows(X)

  coarse = rows(clustid)
  if (coarse) {
    info = panelsetup(clustid,1)
    N_G = rows(info)
  } else {
    info = (1::N), (1::N)  // info = J(0,0,0)
    N_G = N
  }
  
  W = X:>=0  // RHS var of interest: ITT

  Zr = Z, J(N,1,1), X, W:*X  // expand Z to all vars to be partialled out; linear replication stage 
  Zd = X:^2; Zd = Zr, Zd, W:*Zd  // expand Z to all vars to be partialled out; quadratic DGP stage
  ZdKd = Zd :* Kd
  invZZd = invsym(ZZd = cross(ZdKd, Zd))
  invZZr = invsym(ZZr = cross(Zr, Kr, Zr))
  ZWd = cross(ZdKd, W); ZWr = cross(Zr, Kr, W)
  WWd = cross(W, Kd, W); WWr = cross(W , Kr, W)    
  Wd = W - Zd * invZZd * ZWd
  WdKd = Wd :* Kd
  WrKr = (W - Zr * invZZr * ZWr) :* Kr
  halfWrKr = WrKr/2

  ZdinvZZd = Zd * invZZd
  MZdWrKr = WrKr - ZdKd * cross(ZdinvZZd, WrKr)

  if (!fuzzy)
    WWr = WWr - ZWr ' invZZr * ZWr
  WWd   = WWd - ZWd ' invZZd * ZWd  // only cross-products FWL'd
  
  // jackknifing prep
  if (rows(clustid)) {
//     u1ddot[2].M = J(parent->Nobs, 1, 0)
//     Xg = smatrix(parent->Nstar)
//     if (parent->granularjk)
//       invMg = Xg
//     else
//       XXg = XinvHg = smatrix(parent->Nstar)
//     negR1AR1 = - *invH
//     for (g=parent->Nstar; g; g--) {
//       S = parent->NClustVar? parent->infoBootData[g,1] \ parent->infoBootData[g,2] : g\g
//       Xg[g].M = *pXS(*parent->pX1,S)
//       if (parent->granularjk) {  // for many small clusters, faster to compute jackknife errors vaia hc3-like formula
//         invMg[g].M = Xg[g].M * negR1AR1 * Xg[g].M'
//         _diag(invMg[g].M, diagonal(invMg[g].M) :+ 1)
//         invMg[g].M  = invsym(invMg[g].M)
//       } else {
//         XXg[g].M = cross(Xg[g].M, Xg[g].M)
//         XinvHg[g].M = Xg[g].M * (rows(R1perp)? R1perp * invsym(R1perp ' (H - XXg[g].M) *  R1perp) *  R1perp' : invsym(H - XXg[g].M))
//       }
//     }
  } else {
    M = sqrt(Kd) :* ((Zd, Wd))
    M = 1 :/ (1 :- rowsum(M * fold(invsym((ZZd, ZWd \ ZWd', WWd))) :* M))  // standard hc3 multipliers
  }
}

// bias-correction step. Logically similar to the variance simulation step (next), but diverges with optimization
real colvector WBSRDD::bc(real colvector Y, | real colvector T) {
  real scalar ζhat; real colvector Yd, Td, uddoty, uddott, ζddoty, ζddott; real rowvector ζhatb, Wrustby, Wystb, Wrustbt, Wtstb; real matrix v, ustby, ustbt

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  ζddoty = cross(WdKd,Yd) / WWd
  uddoty = Yd - Wd * ζddoty; if (jk) uddoty = uddoty :* M
  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    ζddott = cross(Wd,Kd,Td) / WWd
    uddott = Td - Wd * ζddott; if (jk) uddott = uddott :* M
  }

  v = (runiform(N_G, B₁ + 1) :< .5) :- .5  // Rademacher weights for *all* replications; halved for speed of generation
  v[,1] = J(N_G,1,.5)
  
  if (rows(clustid)) {
  } else {
    ustby = uddoty :* v
    Wrustby = cross(MZdWrKr,ustby) // Wr'u^{*b}_y computed instead of u^{*b}_y alone, for speed
    Wystb = (cross(halfWrKr,Y) - Wrustby[,1]) :+ Wrustby  // Wr ' (*un*-FWL'd endog vars); "half" compensates for Rademacher-halving trick
    if (fuzzy) {
      ustbt = uddott :* v
      Wrustbt = cross(MZdWrKr,ustbt)
      Wtstb = (cross(halfWrKr,T) - Wrustbt[,1]) :+ Wrustbt
    }
  }

  ζhatb = Wystb :/ (fuzzy? Wtstb : WWr/2)  // replication regressions; "/2" compensates for Rademacher-halving trick
  ζhat = ζhatb[1]  // original-sample linear estimate
  return(2 * ζhat - (rowsum(ζhatb)-ζhat)/B₁)  // bias-corrected estimate
}

// variance simulation step
real colvector WBSRDD::vs(real colvector Y, | real colvector T) {
  real scalar ζst, b; real colvector Yd, Td, uddoty, uddott, ζddoty, ζddott; real matrix v, ustby, ustbt, ystb, tstb

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  ζddoty = cross(WdKd,Yd) / WWd
  uddoty = Yd - Wd * ζddoty
  
  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    ζddott = cross(Wd,Kd,Td) / WWd
    uddott = Td - Wd * ζddott
  }

  v = 2 * (runiform(N_G, B₂ + 1) :< .5) :- 1  // Rademacher weights for *all* replications
  v[,1] = J(N_G,1,1)
  
  if (rows(clustid)) {
  } else {  // XXX eliminate O(N) operations?
    ustby = uddoty :* v; ustby = ustby - ZdinvZZd * cross(Zd,Kd,ustby)
    ystb = (Y - ustby[,1]) :+ ustby
    if (fuzzy) {
      ustbt = uddott :* v; ustbt = ustbt - ZdinvZZd * cross(Zd,Kd,ustbt)
      tstb = (T - ustbt[,1]) :+ ustbt
    }
  }

  ζst = cross(WrKr,Y) :/ (fuzzy? cross(WrKr,T) : WWr)  // in 1st entry, ζ*, uncorrected estimate for "original" sample at variance-simulating level
  retval[1] = fuzzy? bc(ystb[,1], tstb[,1]) : bc(ystb[,1])  // bias-corrected estimate in "original" bs sample, ζ* - Δ*ᵢ
  for (b=2; b<=B₂+1; b++)
    retval[b] = (fuzzy? bc(ystb[,b], tstb[,b]) : bc(ystb[,b])) - ζst // deviations of bias-corrected estimates from uncorrected estimate in "original" bs sample ζ̂ ᵢ - Δ**ᵢ - ζ* (algorithm 3.2, step 3)
  return(retval)
}

mata mlib create lbsbc, dir("`c(sysdir_plus)'l") replace  // compile to library because rdrobust clears Mata
mata mlib add lbsbc *(), dir("`c(sysdir_plus)'l")
mata mlib index
end

cap program drop sim
program define sim, rclass
  syntax, μ(string asis) [n(integer 1000) c(real .9) ρ(real 0) ζ(real .04) g(integer 1000) B₁(integer 500) B₂(integer 999)]

  drop _all
  set obs `n'
  mat C = 1, `ρ' \ `ρ', 1
  drawnorm double ut uy, corr(C)
  replace uy = .1295 * uy
  gen double X = 2 * rbeta(2,4) - 1
  gen byte T = ut <= invnormal(cond(X<0, .5-`c'/2, .5+`c'/2))
  gen double Y = `ζ'*T + uy + `μ'
  
  if `g' < `n' {
    gen int clustid = floor((_n-1)/`g'') + 1
    local vceopt vce(cluster clustid)
    local clustidopt st_data(.,"clustid")
  }
  else local clustidopt J(0,1,0)

  cap noi rdrobust Y X, fuzzy(T) bwselect(cerrd) all `vceopt' // occassionally crashes...
  if _rc exit

  return scalar ζhatCL = _b[Conventional]
  return scalar ζhatRBC = _b[Robust]
  return scalar ILCL = (e(ci_r_cl) - e(ci_l_cl)) / 2
  return scalar ILRBC = (e(ci_r_rb) - e(ci_l_rb)) / 2
  return scalar hCEO = e(h_l)
  return scalar bCEO = e(b_l)

  keep if abs(X) < max(e(b_l), e(h_l))
  gen double Kd = max(0, 1 - abs(X) / e(b_l))  // triangular kernels with rdrobust bw's
  gen double Kr = max(0, 1 - abs(X) / e(h_l))

  if `g' < `n' {
    egen long _clustid = group(clustid)  // recalculate ordinal cluster index in case whole clusters fall outside kernels
    drop clustid
    rename _clustid = clustid
    sort clustid
  }

  mata M = WBSRDD()
  mata M.Prep(`b₁', `b₂', `clustidopt', st_data(.,"X"), J(`=_N',0,0), st_data(.,"Kd"), st_data(.,"Kr"), fuzzy=1, jk=1)
  mata dist = M.vs(st_data(.,"Y"), st_data(.,"T"))
  mata distCDR = abs(dist[|2\.|]); _sort(distCDR,1)
  mata st_numscalar("return(ζhatWBS)", dist[1])
  mata st_numscalar("return(ILWBS)", 2 * distCDR[round(.95 * `b₂')])
  
  foreach estimator in CL RBC WBS {
    return scalar EC`estimator' = `ζ' > return(ζhat`estimator') - return(IL`estimator')/2 & `ζ' < return(ζhat`estimator') + return(IL`estimator')/2
  }
  ereturn clear
end

set seed 1231
// sim, ζ(.04) ρ(.9) μ(cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X))
// ret list

cap postclose HBTable1
postfile HBTable1 double ζ ρ DGP ζhatCL ζhatRBC ζhatWBS ECCL ECRBC ECWBS ILCL ILRBC ILWBS bCEO hCEO SDCL SDRBC SDWBS using HBTable1, every(1) replace
parallel init 6
foreach ρ in 0 .9 -.9 {
  forvalues DGP=1/3 {
    local ζ: word `DGP' of .04 -3.45 .04
    local μ: word `DGP' of "cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X)" ///
                           "cond(X<0, (2.3+(3.28+(1.45+(.23+.03*X)*X)*X)*X)*X, (18.49+(-54.81+(74.3+(-45.0+9.83*X)*X)*X)*X)*X)" ///
                           "cond(X<0, (1.27+(3.59+(14.147+(23.694+10.995*X)*X)*X)*X)*X, (.84+(-.3+(2.397+(-.901+3.56*X)*X)*X)*X)*X)"
    parallel sim, reps(100) nodots proc(2) exp(ζhatCL=r(ζhatCL) ILCL=r(ILCL) ECCL=r(ECCL) ζhatRBC=r(ζhatRBC) ILRBC=r(ILRBC) ECRBC=r(ECRBC) ζhatWBS=r(ζhatWBS) ILWBS=r(ILWBS) ECWBS=r(ECWBS) bCEO=r(bCEO) hCEO=r(hCEO)): ///
      sim, ζ(`ζ') ρ(`ρ') μ(`μ')
    collapse ζhat* EC* IL* *CEO (sd) SDCL=ζhatCL SDRBC=ζhatRBC SDWBS=ζhatWBS
    post HBTable1 (`ζ') (`ρ') (`DGP') (ζhatCL) (ζhatRBC) (ζhatWBS) (ECCL) (ECRBC) (ECWBS) (ILCL) (ILRBC) (ILWBS) (bCEO) (hCEO) (SDCL) (SDRBC) (SDWBS)
  }
}
postclose HBTable1
use HBTable1, clear
foreach estimator in CL RBC WBS {
  gen bias`estimator' = ζhat`estimator' - ζ
  gen RMSE`estimator' = sqrt(bias`estimator'^2 + SD`estimator'^2)
  drop ζhat`estimator'
}
reshape long bias SD RMSE EC IL, i(ρ DGP *CEO) j(estimator) string
list ρ DGP estimator bias SD RMSE EC IL *CEO, sep(0)
