cd D:\OneDrive\Documents\Macros\rdboottest

mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

struct smatrix {
  real matrix M
}

class WBSRDD {
  real scalar fuzzy, N, WWd, WWr, coarse, N_G, B₁, B₂, jk, granularjk, m, dirty, bc, ζst, ζstbc
  real colvector W, ZWd, ZWr, Wd, WrKr, halfWrKr, X, Kd, Kr, MZdWrKr, WdKd, clustid, dist
  real matrix info, Zr, Zd, invZZd, invZZr, ZdinvZZd, Z, ZdKd
  struct smatrix colvector invMg, Xg, XXg, XinvHg, uddoty, uddott

  void Prep(), vs()
  real colvector getdist()
  real rowvector getci()
  real scalar getp()
  private real colvector bc()
  private void new(), jk()
  private real matrix fold()
}

void WBSRDD::new()
  jk = bc = 1

real matrix WBSRDD::fold(matrix X) return(uppertriangle(X) + lowertriangle(X,0)')  // fold matrix diagonally; returns same values as a quad form, but runs faster because of all the 0's

real colvector triangularkernel  (real scalar bw, real colvector X) return(1:-abs(X)/bw)
real colvector epanechnikovkernel(real scalar bw, real colvector X) return((1 :- (X/bw):^2 / 5) * .75 / sqrt(5))
real colvector uniformkernel     (real scalar bw, real colvector X) return(J(rows(X), 1, 1/bw))
  

// one-time stuff that only depends on exog vars and cluster and kernel definitions
void WBSRDD::Prep(real scalar B₁, real scalar B₂, real colvector clustid, real colvector X, real matrix Z, real colvector wt, real scalar h_l, real scalar h_r, real scalar b_l, real scalar b_r, string scalar kernel, 
                  real scalar fuzzy, real scalar bc, real scalar jk) {
  real matrix ZZd, H, invH, neginvH, Z_W; real scalar g, kZ; pointer(real colvector function) kernelfn

  this.B₁ = B₁
  this.B₂ = B₂
  this.fuzzy = fuzzy
  this.jk = jk
  this.bc = bc
  this.clustid = clustid
  this.X = X
  this.Z = Z

  W = X:>=0  // RHS var of interest: ITT

  kernelfn = kernel=="triangular"? &triangularkernel() : (kernel=="epanechnikov"? &epanechnikovkernel() : &uniformkernel())
  this.Kd = (*kernelfn)(b_l, X) :* (X:<0) :* (X:>-b_l) + (*kernelfn)(b_r, X) :* (X:>=0) :* (X:<b_r)
  this.Kr = (*kernelfn)(h_l, X) :* (X:<0) :* (X:>-h_l) + (*kernelfn)(h_r, X) :* (X:>=0) :* (X:<h_r)

  if (rows(wt)) {
    this.Kd = this.Kd :* wt
    this.Kr = this.Kr :* wt
  }

  N = rows(X)

  uddoty = uddott = smatrix(1 + jk)  // for jackknife, need jk'd residuals but also non-jk'd residuals for original test stat

  coarse = rows(clustid)
  if (coarse) {
    info = panelsetup(clustid,1)
    N_G = rows(info)
  } else {
    info = (1::N), (1::N)  // info = J(0,0,0)
    N_G = N
  }
  
  Zr = Z, J(N,1,1), X, W:*X  // expand Z to all vars to be partialled out; linear replication stage 
  Zd = X:^2; Zd = Zr, Zd, W:*Zd  // expand Z to all vars to be partialled out; quadratic DGP stage
  ZdKd = Zd :* Kd
  invZZd = invsym(ZZd = cross(ZdKd, Zd))
  invZZr = invsym(      cross(Zr, Kr, Zr))
  ZWd = cross(ZdKd, W); ZWr = cross(Zr, Kr, W)
  WWd = cross(W, Kd, W); WWr = cross(W, Kr, W)

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
  if (jk) {
    Z_W = sqrt(Kd) :* ((Zd, Wd))  // all RHS vars in DGP regressions, sqrt weights folded in
    kZ = cols(Zd)
    
    uddoty[2].M = J(N,1,0)  // will hold jk'd residuals
    if (fuzzy) uddott[2].M = J(N,1,0)

    invH = invsym(H = (ZZd, ZWd \ ZWd', WWd))
    
    if (rows(clustid)) {
      m = sqrt((N_G - 1) / N_G)
      granularjk = (1+kZ)^3 + N_G * (N/N_G*(1+kZ)^2 + (N/N_G)^2*(1+kZ) + (N/N_G)^2 + (N/N_G)^3) < N_G * ((1+kZ)^2*N/N_G + (1+kZ)^3 + 2*(1+kZ)*(1+kZ + N/N_G))
      Xg = smatrix(N_G)
      if (granularjk)
        invMg = Xg
      else
        XXg = XinvHg = Xg
      neginvH = -invH
      for (g=N_G; g; g--) {  // compute jackknife errors via hc3-like block-diagonal matrix; slower for few clusters than direct jackknifing, but can be folded into once-computed objects
        Xg[g].M = Z_W[|info[g,1],. \ info[g,2],.|]
        if (granularjk) {
          invMg[g].M = Xg[g].M * neginvH * Xg[g].M'
          _diag(invMg[g].M, diagonal(invMg[g].M) :+ 1)
          invMg[g].M  = invsym(invMg[g].M)
        } else {
          XXg[g].M = cross(Xg[g].M, Xg[g].M)
          XinvHg[g].M = Xg[g].M * invsym(H - XXg[g].M)
        }
      }
    } else {
      m = sqrt((N - 1) / N)
      (invMg = smatrix()).M = 1 :/ (1 :- rowsum(Z_W * fold(invH) :* Z_W))  // standard hc3 multipliers
//       MZdWrKr = MZdWrKr :* invMg.M  // trick to pre-fold-in jackknifing in bc step
    }
  }
}

// jk-transform 1st entry in a 2-vector of residuals into the 2nd
void WBSRDD::jk(struct smatrix rowvector uddot) {
  real colvector uddotg, S; real scalar g

  if (jk)
    if (rows(clustid)) {
      if (granularjk)
        for (g=N_G; g; g--) {
          S = info[g,1] \ info[g,2]
          uddotg = uddot.M[|S|]
          uddot[2].M[|S|] = m * cross(invMg[g].M, uddotg)
        }
      else
        for (g=N_G; g; g--) {
          S = info[g,1] \ info[g,2]
          uddotg = uddot.M[|S|]
          uddot[2].M[|S|] = m * (uddotg + XinvHg[g].M * cross(Xg[g].M, uddotg))
        }
    } else
      uddot[2].M = m * invMg.M :* uddot.M
}

// bias-correction step. Logically similar to the variance simulation step (next), but diverges with optimization
real colvector WBSRDD::bc(real colvector Y, | real colvector T) {
  real scalar ζhat; real colvector Yd, Td; real rowvector ζhatb, Wrustby, Wystb, Wrustbt, Wtstb; real matrix v

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  uddoty.M = Yd - Wd * (cross(WdKd,Yd) / WWd) // jackknife adjustment already folded into MZdWrKr
  jk(uddoty)
  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    uddott.M = Td - Wd * cross(Wd,Kd,Td) / WWd
    jk(uddott)
  }

  v = (runiform(N_G, B₁ + 1) :< .5) :- .5  // Rademacher weights for *all* replications; halved for speed of generation
  v[,1] = J(N_G,1,.5)

  if (rows(clustid)) {  // boottest trick
    Wrustby = cross(panelsum(MZdWrKr, uddoty[1+jk].M, info), v)  // Wr'u^{*b}_y computed instead of u^{*b}_y alone, for speed
    if (jk) Wrustby[1] = cross(MZdWrKr, uddoty.M) / 2  // fix: original-sample ust is non-jk'd uddot; "/ 2" for consistency with Rademacher-halving
    Wystb = (cross(halfWrKr,Y) - Wrustby[1]) :+ Wrustby  // Wr ' (*un*-FWL'd endog vars); "half" compensates for Rademacher-halving trick
    if (fuzzy) {
      Wrustbt = cross(panelsum(MZdWrKr, uddott.M, info), v)
      if (jk) Wrustbt[1] = cross(MZdWrKr, uddott.M) / 2
      Wtstb = (cross(halfWrKr,T) - Wrustbt[1]) :+ Wrustbt
    }
  } else {
    Wrustby = cross(MZdWrKr, uddoty[1+jk].M, v) // Wr'u^{*b}_y computed instead of u^{*b}_y alone, for speed
    if (jk) Wrustby[1] = cross(MZdWrKr, uddoty.M) / 2
    Wystb = (cross(halfWrKr,Y) - Wrustby[1]) :+ Wrustby  // Wr ' (un-FWL'd endog vars); "half" compensates for Rademacher-halving trick
    if (fuzzy) {
      Wrustbt = cross(MZdWrKr, uddott[1+jk].M, v)
      if (jk) Wrustbt[1] = cross(MZdWrKr, uddott.M) / 2
      Wtstb = (cross(halfWrKr,T) - Wrustbt[1]) :+ Wrustbt
    }
  }

  ζhatb = Wystb :/ (fuzzy? Wtstb : WWr/2)  // replication regressions; "/2" compensates for Rademacher-halving trick
  ζhat = ζhatb[1]  // original-sample linear estimate
  return(2 * ζhat - (rowsum(ζhatb)-ζhat)/B₁)  // bias-corrected estimate
}

// variance simulation step
void WBSRDD::vs(real colvector Y, | real colvector T) {
  real scalar b; real colvector Yd, Td; real matrix v, _v, ustby, ustbt, ystb, tstb

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  uddoty.M = Yd - Wd * (cross(WdKd,Yd) / WWd)
  jk(uddoty)
  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    uddott.M = Td - Wd * (cross(WdKd,Td) / WWd)
    jk(uddott)
  }

  v = 2 * (runiform(N_G, B₂ + 1) :< .5) :- 1  // Rademacher weights for *all* replications
  v[,1] = J(N_G,1,1)

  if (rows(clustid)) {  // partial use of boottest trick
    _v = v[clustid,]
    ustby = uddoty[1+jk].M :* _v - ZdinvZZd * cross(panelsum(uddoty[1+jk].M :* ZdKd, info), v)
    if (jk) ustby[,1] = uddoty.M  // fix: original-sample ust is non-jk'd uddot
    ystb = (Y - uddoty.M) :+ ustby
    if (fuzzy) {
      ustbt = uddott[1+jk].M :* _v - ZdinvZZd * cross(panelsum(uddott[1+jk].M :* ZdKd, info), v)
      if (jk) ustbt[,1] = uddott.M
      tstb = (T - uddott.M) :+ ustbt
    }
  } else {
    ustby = uddoty[1+jk].M :* v; ustby = ustby - ZdinvZZd * cross(ZdKd,ustby)
    if (jk) ustby[,1] = uddoty.M  // fix: original-sample ust is non-jk'd uddot
    ystb = (Y - uddoty.M) :+ ustby
    if (fuzzy) {
      ustbt = uddott.M :* v; ustbt = ustbt - ZdinvZZd * cross(ZdKd,ustbt)
      if (jk) ustbt[,1] = uddott.M 
      tstb = (T - uddott.M) :+ ustbt
    }
  }

  ζst = cross(WrKr,Y) :/ (fuzzy? cross(WrKr,T) : WWr)  // ζ*, uncorrected estimate for "original" sample at variance-simulating level
"ζst"
ζst
  if (bc) {
    dist = J(B₂, 1, 0)
    ζstbc = fuzzy? bc(ystb[,1], tstb[,1]) : bc(ystb[,1])  // bias-corrected estimate in "original" bs sample, ζ* - Δ*ᵢ
    for (b=2; b<=B₂+1; b++)
      dist[b-1] = (fuzzy? bc(ystb[,b], tstb[,b]) : bc(ystb[,b])) - ζst // deviations of bias-corrected estimates from uncorrected estimate in "original" bs sample ζ̂ ᵢ - Δ**ᵢ - ζ* (algorithm 3.2, step 3)
  } else {
    ζstbc = ζst
    dist = (cross(ystb, WrKr) :/ (fuzzy? cross(tstb, WrKr) : WWr))[|2\.|] :- ζst
  }
  _sort(dist,1)
}

real colvector WBSRDD::getdist() return(dist)

real scalar WBSRDD::getp(| string scalar ptype) {
  if (ptype=="symmetric" | ptype=="")
    return(colsum(-abs(ζstbc) :> -abs(ζstbc :+ dist)) / B₂)
  if (ptype=="equaltail")
    return(2 * min((colsum(dist :< 0) , colsum(dist :> 0))) / B₂)
  if (ptype=="lower")
    return(colsum( dist :< 0) / B₂)
  if (ptype=="upper")
    return(colsum(dist :> 0) / B₂)
  _error(198, `"p type must be "symmetric", "equaltail", "lower", or "upper"."')
}

real rowvector WBSRDD::getci(real scalar level, | string scalar citype) {
  real scalar halfwidth
  if (citype=="symmetric" | citype=="") {
    halfwidth = abs(dist)[round(level/100 * B₂)]
    return((ζstbc-halfwidth, ζstbc+halfwidth))
  }
  if (citype=="equaltail")
    return((ζstbc-dist[round((1-level/100)/2 * B₂)], ζstbc+dist[round((1-(1-level/100)/2) * B₂)]))
  _error(198, `"CI type must be "symmetric" or "equaltail"."')
}

mata mlib create lbsbc, dir("`c(sysdir_plus)'l") replace  // compile to library because rdrobust clears Mata
mata mlib add lbsbc *(), dir("`c(sysdir_plus)'l")
mata mlib index
end

cap program drop rdboottest
program define rdboottest, eclass
  local cmdline `0'

  syntax varlist [if] [in], [c(real 0) scalepar(real 1) fuzzy(varname) weights(string) seed(string) nojk nobc BCREPs(integer 1000) REPs(integer 999) Level(real `c(level)') PType(string) *]
  
  local jk = "`nojk'"==""
  local bc = "`nobc'"==""

  marksample touse
  markout `touse' `weights' `fuzzy'
  
  rdrobust `varlist' `if' `in', c(`c') scalepar(`scalepar') fuzzy(`fuzzy') weights(`weights') `options'
  ereturn local fuzzy `fuzzy'

  if `:word count `ptype'' > 1 {
		di as err "The {cmd:wp:type} option must be {cmdab:sym:metric}, {cmdab:eq:qualtail}, {cmd:lower}, or {cmd:upper}."
		exit 198
	}

	if "`ptype'"'=="" local ptype symmetric
	else {
		local 0, `ptype'
		syntax, [SYMmetric EQualtail LOWer UPper]
		local ptype `symmetric'`equaltail'`lower'`upper'
	}

  preserve
  qui keep if `touse' & `e(runningvar)' > -min(e(h_l), e(b_l)) & `e(runningvar)' < min(e(h_r), e(b_r))
  
  if `c' qui replace `e(runningvar)' = `e(runningvar)' - `c'

  if "`e(clustvar)'" != "" {
    tempname clustid
    egen long `clustid' = group(`e(clustvar)')
    sort `clustid', stable
    local clustidopt st_data(.,"`clustid'")
  }
  else local clustidopt J(0,1,0)
  
  local covsopt = cond("`e(covs)'"=="", "J(`=_N',0,0)", "`st_data(.,"`e(covs)'")'")
  
  local wtopt = cond(`"`weights'"'=="", "J(0,1,0)", "`st_data(.,"`weights'")'")

  if `"`seed'"'!="" set seed `seed'
  mata M = WBSRDD()
  mata M.Prep(`bcreps', `reps', `clustidopt', st_data(.,"`e(runningvar)'"), `covsopt', `wtopt', st_numscalar("e(h_l)"), st_numscalar("e(h_r)"), st_numscalar("e(b_l)"), st_numscalar("e(b_r)"), "triangular", "`fuzzy'"!="", `bc', `jk')
  mata "`fuzzy'"=="" ? M.vs(st_data(.,"`e(depvar)'")) : M.vs(st_data(.,"`e(depvar)'"), st_data(.,"`fuzzy'"))
  restore

// mata st_numscalar("e(tau_cl_wb)", M.ζst   * `scalepar')
  mata st_numscalar("e(tau_bc_wb)", M.ζstbc * `scalepar')
  mata st_numscalar("e(p_wb)", M.getp("`ptype'"))
  mata CI = M.getci(`level', "`ptype'")
  mata st_numscalar("e(ci_l_rb_wb)", CI[1] * `scalepar' + `c')
  mata st_numscalar("e(ci_r_rb_wb)", CI[2] * `scalepar' + `c')
  mata st_matrix("e(dist_wb)", M.getdist() * `scalepar' :+ `c')

  di as inp "    Wild bootstrap {c |}" _col(22) as res %7.0g cond(`bc', e(tau_bc_wb), e(tau_bc_wb)) _col(52) %4.3f e(p_wb) _col(60) %8.0g e(ci_l_rb_wb) _col(73) %8.0g e(ci_r_rb_wb)
  di as inp "{hline 19}{c BT}{hline 60}" _n

  ereturn repost, esample(`touse')
  ereturn local cmdline `cmdline'
  ereturn local cmd rdboottest
end

cap program drop sim
program define sim, rclass
  syntax, μ(string asis) [n(integer 1000) c(real .9) ρ(real 0) ζ(real .04) g(integer 1000) bc(integer 1) jk(integer 1) B₁(integer 500) B₂(integer 999)]

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

  cap noi rdrobust Y X, fuzzy(T) bwselect(cerrd) all `vceopt' // occassionally crashes...
  if _rc exit

  return scalar ζhatCL = _b[Conventional]
  return scalar ζhatRBC = _b[Robust]
  return scalar ILCL  = (e(ci_r_cl) - e(ci_l_cl)) / 2
  return scalar ILRBC = (e(ci_r_rb) - e(ci_l_rb)) / 2
  return scalar hCEO = e(h_l)
  return scalar bCEO = e(b_l)

  keep if abs(X) < max(e(b_l), e(h_l))

  if `g' < `n' {
    egen long _clustid = group(clustid)  // recalculate ordinal cluster index in case whole clusters fall outside kernels
    drop clustid
    rename _clustid clustid
    sort clustid
  }

  mata M = WBSRDD()
  mata M.Prep(`b₁', `b₂', `clustidopt', st_data(.,"X"), J(`=_N',0,0), J(0,1,0), st_numscalar("e(b_l)"), st_numscalar("e(b_r)"), st_numscalar("e(h_l)"), st_numscalar("e(h_r)"), "triangular", fuzzy=1, `bc', `jk')
  mata M.vs(st_data(.,"Y"), st_data(.,"T"))
  mata st_numscalar("return(ζhatWBS)", M.ζstbc)
  mata CI = M.getci(95); st_numscalar("return(ILWBS)", CI[2]-CI[1])
  
  foreach estimator in CL RBC WBS {
    return scalar EC`estimator' = `ζ' > return(ζhat`estimator') - return(IL`estimator')/2 & `ζ' < return(ζhat`estimator') + return(IL`estimator')/2
  }
  ereturn clear
end

set seed 1231
local n 1000
local ρ .9
local μ cond(X<0, (1.27+(7.18+(20.21+(21.54+7.33*X)*X)*X)*X)*X, (.84+(-3+(7.99+(-9.01+3.56*X)*X)*X)*X)*X)
local ζ .04
local c .9
drop _all
set obs `n'
mat C = 1, `ρ' \ `ρ', 1
drawnorm double ut uy, corr(C)
replace uy = .1295 * uy
gen double X = 2 * rbeta(2,4) - 1
gen byte T = ut <= invnormal(cond(X<0, .5-`c'/2, .5+`c'/2))
gen double Y = `ζ'*T + uy + `μ'
gen W = X >=0 
sort X

rdboottest Y X, all 
// qui {
// scalar h_r = e(h_r)
// gen wt = max(0,1-abs(X/h_r))
// ivregress 2sls Y W##c.X [aw=wt] if abs(X) < h_r
// mvreg Y W = W#c.X [aw=wt] if abs(X) < h_r
// predict Ye if e(sample), resid
// predict We if e(sample), resid eq(#2)
// reg Ye We [aw=wt], nocons
// gen numer = Ye * We * wt if e(sample)
// gen denom = We * We * wt if e(sample)
// }
// sum numer
// di r(sum)
// sum denom
// di r(sum)
--

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

