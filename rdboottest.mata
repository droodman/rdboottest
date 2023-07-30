mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

// real matrix _panelsum(X, W, info)
//   return (rows(info) & rows(info) < rows(X)? panelsum(X,W,info) : X :* W)

class WBSRDD {
  real scalar fuzzy, N, WWd, WWr, coarse, N_G, B₁, B₂
  real colvector W, ZWd, ZWr, Wd, WrKr, clustid, X, Kd, Kr, retval
  real matrix info, Zr, Zd, invZZd, invZZr, ZdinvZZd, Z

  real colvector wbsrdd()
  void Prep()
}

// one-time stuff that only depends on exog vars and cluster and kernel definitions
void WBSRDD::Prep(real scalar B₁, real scalar B₂, real colvector clustid, real colvector X, real matrix Z, real colvector Kd, real colvector Kr, real scalar fuzzy) {
  this.B₁ = B₁
  this.B₂ = B₂
  this.fuzzy = fuzzy
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
  invZZd = invsym(cross(Zd, Kd, Zd))
  invZZr = invsym(cross(Zr, Kr, Zr))
  ZWd = cross(Zd, Kd, W); ZWr = cross(Zr, Kr, W)
  WWd = cross(W , Kd, W); WWr = cross(W , Kr, W)    
  Wd = W - Zd * invZZd * ZWd
  WrKr = (W - Zr * invZZr * ZWr) :* Kr

  ZdinvZZd = Zd * invZZd

  if (!fuzzy)
    WWr = WWr - ZWr ' invZZr * ZWr
  WWd   = WWd - ZWd ' invZZd * ZWd  // only cross-products FWL'd
}

real colvector WBSRDD::wbsrdd(real scalar bc, real colvector Y, | real colvector T) {
  real scalar ζhat, ζst, b; real colvector Yd, Td, uddoty, uddott, ζddoty, ζddott; real rowvector ζhatb; real matrix v, ustby, ustbt, ystb, tstb

  Yd = Y - ZdinvZZd * cross(Zd, Kd, Y)  // FWL
  if (fuzzy)
    Td = T - ZdinvZZd * cross(Zd, Kd, T)

  ζddoty = cross(Wd,Kd,Yd) / WWd
  uddoty = Yd - Wd * ζddoty
  if (fuzzy) {
    ζddott = cross(Wd,Kd,Td) / WWd
    uddott = Td - Wd * ζddott
  }

  v = 2 * (runiform(N_G, (bc ? B₁ : B₂) + 1) :< .5) :- 1  // Rademacher weights for *all* replications  XXX can halve these for speed?
  v[,1] = J(N_G,1,1)

  if (N_G < N) {
  } else {  // XXX eliminate O(N) operations?
    ustby = uddoty :* v; ustby = ustby - ZdinvZZd * cross(Zd,Kd,ustby)
    ystb = (Y - ustby[,1]) :+ ustby
    if (fuzzy) {
      ustbt = uddott :* v; ustbt = ustbt - ZdinvZZd * cross(Zd,Kd,ustbt)
      tstb = (T - ustbt[,1]) :+ ustbt
    }
  }
  
  if (bc) { // bias-correcting step
    ζhatb = cross(WrKr,ystb) :/ (fuzzy? cross(WrKr,tstb) : WWr)  // replication regressions
    ζhat = ζhatb[1]  // original-sample linear estimate
    return(2 * ζhat - (rowsum(ζhatb)-ζhat)/B₁)  // bias-corrected estimate
  } else {  // variance simulation step
    ζst = cross(WrKr,Y) :/ (fuzzy? cross(WrKr,T) : WWr)  // in 1st entry, ζ*, uncorrected estimate for "original" sample at variance-simulating level
    retval[1] = fuzzy? wbsrdd(1, ystb[,1], tstb[,1]) : wbsrdd(1, ystb[,1])  // bias-corrected estimate in "original" bs sample, ζ* - Δ*ᵢ
    for (b=2; b<=B₂+1; b++)
      retval[b] = (fuzzy? wbsrdd(1, ystb[,b], tstb[,b]) : wbsrdd(1, ystb[,b])) - ζst // deviations of bias-corrected estimates from uncorrected estimate in "original" bs sample ζ̂ ᵢ - Δ**ᵢ - ζ* (algorithm 3.2, step 3)
    return(retval)
  }
}

mata mlib create lbsbc, dir("`c(sysdir_plus)'l") replace  // compile to library because rdrobust clears Mata
mata mlib add lbsbc *(), dir("`c(sysdir_plus)'l")
mata mlib index
end

set seed 1231

cap program drop sim
program define sim, rclass
  syntax, yfn(string asis) [n(integer 1000) c(real .9) ρ(real 0) ζ(real .04) B₁(integer 500) B₂(integer 999)]

  drop _all
  set obs `n'
  mat C = 1, ρ \ ρ, 1
  drawnorm double ut uy, corr(C)
  replace uy = .1295 * uy
  gen double X = 2 * rbeta(2,4) - 1
  gen byte T = ut <= invnormal(cond(X<0, .5-c/2, .5+c/2))

  gen double Y = ζ*T + rnormal() + `yfn'

  rdrobust Y X, fuzzy(T) bwselect(cerrd) all
  return scalar CL_ζhat = _b[Conventional]
  return scalar RBC_ζhat = _b[Robust]
  return scalar CL_IL = (e(ci_r_cl) - e(ci_l_cl)) / 2
  return scalar RBC_IL = (e(ci_r_rb) - e(ci_l_rb)) / 2

  keep if abs(X) < max(e(b_l), e(h_l))
  gen double Kd = max(0, 1 - abs(X) / e(b_l))  // triangular kernels with rdrobust bw's
  gen double Kr = max(0, 1 - abs(X) / e(h_l))

  mata M = WBSRDD()
  mata M.Prep(`b₁', `b₂', J(0,1,0), st_data(.,"X"), J(`=_N',0,0), st_data(.,"Kd"), st_data(.,"Kr"), 1)
  mata dist = M.wbsrdd(0, st_data(.,"Y"), st_data(.,"T"))
  mata distCDR = abs(dist[|2\.|]); _sort(distCDR,1)
  mata st_numscalar("return(WBS_ζhat)", dist[1])
  mata st_numscalar("return(WBS_IL)", 2 * distCDR[round(.95 * `b₂')])
  
  foreach estimator in CL RBC WBS {
    return scalar `estimator'_EC = `ζ' > return(`estimator'_ζhat) - return(`estimator'_IL) / 2 & `ζ' < return(`estimator'_ζhat) + return(`estimator'_IL) / 2
  }
end

sim, ζ(.04) ρ(.9) yfn(cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X))

mata M.wbsrdd(0, st_data(.,"Y"), st_data(.,"T"))

// parallel init 6
// simulate CL_ζhat=r(CL_ζhat) CL_IL=r(CL_IL) CL_EC=r(CL_EC) RBC_ζhat=r(RBC_ζhat) RBC_IL=r(RBC_IL) RBC_EC=r(RBC_EC) WBS_ζhat=r(WBS_ζhat) WBS_IL=r(WBS_IL) WBS_EC=r(WBS_EC), nodots reps(3): sim, ζ(.04) ρ(.9) yfn(cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X))
// // scalar ζ = .04
// // scalar ζ = -3.45
// // scalar ζ = .04
//
// // gen double Y = ζ*T + cond(X<0, (2.3+(3.28+(1.45+(.23+.03*X)*X)*X)*X)*X, (18.49+(-54.81+(74.3+(-45.0+9.83*X)*X)*X)*X)*X) + rnormal()
// // gen double Y = ζ*T + cond(X<0, (1.27+(3.59+(14.147+(23.694+10.995*X)*X)*X)*X)*X, (.84+(-.3+(2.397+(-.901+3.56*X)*X)*X)*X)*X) + rnormal()
//
