mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

struct smatrix {
  real matrix M
}

class WBSRDD {
  real scalar fuzzy, N, WWd, WWr, N_G, B1, B2, jk, granularjk, m, dirty, bc, zetast, zetastbc, p, q
  real colvector W, ZWd, ZWr, Wd, WrKr, halfWrKr, MZdWrKr, MZdWrKrjk, WdKd, clustid, dist
  real matrix info, Zr, Zd, invZZd, invZZr, ZdinvZZd, ZdKd, vbc
  struct smatrix colvector invMg, Xg, XXg, XinvHg, uddoty, uddott

  void Prep(), vs()
  real colvector getdist()
  real rowvector getci()
  real scalar getp()
  private real colvector bc()
  private void jk()
  private real matrix fold()
}

real matrix WBSRDD::fold(matrix X) return(uppertriangle(X) + lowertriangle(X,0)')  // fold matrix diagonally; returns same values as a quad form, but runs faster because of all the 0's

real colvector triangularkernel  (real scalar bw, real colvector X) return(1:-abs(X)/bw)
real colvector epanechnikovkernel(real scalar bw, real colvector X) return((1 :- (X/bw):^2 / 5) * .75 / sqrt(5))
real colvector uniformkernel     (real scalar bw, real colvector X) return(J(rows(X), 1, 1/bw))
  

// one-time stuff that only depends on exog vars and cluster and kernel definitions
void WBSRDD::Prep(real scalar p, real scalar q, real scalar B1, real scalar B2, real colvector clustid, real colvector X, real matrix Z, real colvector wt, real scalar h_l, real scalar h_r, real scalar b_l, real scalar b_r, string scalar kernel, 
                  real scalar fuzzy, real scalar bc, real scalar jk) {
  real colvector tmp, Kd, Kr; real matrix ZZd, H, invH, neginvH, Z_W, Xp; real scalar g, kZ, invm; pointer(real colvector function) kernelfn

  this.B1 = B1
  this.B2 = B2
  this.fuzzy = fuzzy
  this.jk = jk
  this.bc = bc
  this.clustid = clustid

  W = X:>=0  // RHS var of interest: ITT

  kernelfn = kernel=="triangular"? &triangularkernel() : (kernel=="epanechnikov"? &epanechnikovkernel() : &uniformkernel())
  Kd = (*kernelfn)(b_l, X) :* (X:<0) :* (X:>-b_l) + (*kernelfn)(b_r, X) :* (X:>=0) :* (X:<b_r)
  Kr = (*kernelfn)(h_l, X) :* (X:<0) :* (X:>-h_l) + (*kernelfn)(h_r, X) :* (X:>=0) :* (X:<h_r)

  if (rows(wt)) {
    Kd = Kd :* wt
    Kr = Kr :* wt
  }

  N = rows(X)

  uddoty = uddott = smatrix(1 + jk)  // for jackknife, need jk'd residuals but also non-jk'd residuals for original test stat

  if (rows(clustid)) {
    info = panelsetup(clustid,1)
    N_G = rows(info)
  } else {
    info = (1::N), (1::N)  // info = J(0,0,0)
    N_G = N
    vbc = (runiform(N, B1 + 1) :< .5) :- .5  // Rademacher weights for all bias-correcting replications; halved for speed of generation
    vbc[,1] = J(N,1,.5)
  }
  
  tmp = X
  if (p) {  // expand Z to all vars to be partialled out; normally-linear (p=1) replication stage
    Xp = X
    for (g=p-1;g;g--) {
      tmp = tmp :* X
      Xp = Xp, tmp  // polynomials in running var
    }
    Zr = Z, J(N,1,1), Xp, W:*Xp 
  } else
    Zr = Z, J(N,1,1)
  Xp = J(N,0,0)
  for (g=q;g>p;g--) {  // same for DGP stage with higher-order poly in running var
    tmp = tmp :* X
    Xp = Xp, tmp
  }
  Zd = Zr, Xp, W:*Xp  // expand Z to all vars to be partialled out; normally-linear (p=1) replication stage 
  
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

  // jk prep
  if (jk) {
    Z_W = sqrt(Kd) :* ((Zd, Wd))  // all RHS vars in DGP regressions, sqrt weights folded in
    kZ = cols(Zd)
    
    uddoty[2].M = J(N,1,0)  // will hold jk'd residuals
    if (fuzzy) uddott[2].M = J(N,1,0)

    invH = invsym(H = (ZZd, ZWd \ ZWd', WWd))
    
    if (rows(clustid)) {
      granularjk = (1+kZ)^3 + N_G * (N/N_G*(1+kZ)^2 + (N/N_G)^2*(1+kZ) + (N/N_G)^2 + (N/N_G)^3) < N_G * ((1+kZ)^2*N/N_G + (1+kZ)^3 + 2*(1+kZ)*(1+kZ + N/N_G))
      Xg = smatrix(N_G)
      if (granularjk)
        invMg = Xg
      else
        XXg = XinvHg = Xg
      if (granularjk) {  // when clusters are many/small, compute jackknife errors via hc3-like block-diagonal matrix
        invm = sqrt(N_G / (N_G - 1))
        neginvH = -invm * invH
        for (g=N_G; g; g--) {  // compute jackknife errors via hc3-like block-diagonal matrix; slower for few clusters than direct jackknifing, but can be folded into once-computed objects
          Xg[g].M = Z_W[|info[g,1],. \ info[g,2],.|]
          invMg[g].M = Xg[g].M * neginvH * Xg[g].M'
          _diag(invMg[g].M, diagonal(invMg[g].M) :+ invm)
          invMg[g].M  = invsym(invMg[g].M)  // jk small-sample factor m backed in
        }
      } else {  // when clusters are few, prep for faster, direct jack-knifing, i.e., running OLS regressions without each cluster in turn
        m = sqrt((N_G - 1) / N_G)
        neginvH = -invH
        for (g=N_G; g; g--) {
          Xg[g].M = Z_W[|info[g,1],. \ info[g,2],.|]
          XXg[g].M = cross(Xg[g].M, Xg[g].M)  // jk small-sample factor m not baked in
          XinvHg[g].M = Xg[g].M * invsym(H - XXg[g].M)  // jk small-sample factor m not baked in
        }
      }
    } else {
      (invMg = smatrix()).M = sqrt((N - 1) / N) :/ (1 :- rowsum(Z_W * fold(invH) :* Z_W))  // standard hc3 multipliers
      MZdWrKrjk = MZdWrKr :* invMg.M  // trick: when no clustering, jk hat matrix is diagonal; multiply it one time against unchanging factor
    }
  }
}

// jk-transform 1st entry in a 2-vector of residuals into the 2nd
// vs = 1 means variance-simulation stage & no trick to pre-compute jk-ing, so do it here
void WBSRDD::jk(struct smatrix rowvector uddot, real scalar vs) {
  real colvector uddotg, S; real scalar g

  if (jk)
    if (rows(clustid)) {
      if (granularjk)
        for (g=N_G; g; g--) {
          S = info[g,1] \ info[g,2]
          uddotg = uddot.M[|S|]
          uddot[2].M[|S|] = cross(invMg[g].M, uddotg)
        }
      else
        for (g=N_G; g; g--) {
          S = info[g,1] \ info[g,2]
          uddotg = uddot.M[|S|]
          uddot[2].M[|S|] = m * (uddotg + XinvHg[g].M * cross(Xg[g].M, uddotg))
        }
    } else if (vs)
      uddot[2].M = invMg.M :* uddot.M
}

// bias-correction step. Logically similar to the variance simulation step (next), but diverges with optimization
real colvector WBSRDD::bc(real colvector Y, | real colvector T) {
  real scalar zetahat; real colvector Yd, Td; real rowvector zetahatb, Wrustby, Wystb, Wrustbt, Wtstb; real matrix v

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  uddoty.M = Yd - Wd * (cross(WdKd,Yd) / WWd)
  jk(uddoty, 0)

  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    uddott.M = Td - Wd * cross(WdKd,Td) / WWd
    jk(uddott, 0)
  }

  if (rows(clustid)) {  // boottest trick
    v = (runiform(N_G, B1 + 1) :< .5) :- .5  // Rademacher weights for *all* replications; halved for speed of generation
    v[,1] = J(N_G,1,.5)

    Wrustby = cross(panelsum(MZdWrKr, uddoty[1+jk].M, info), v)  // Wr'u^{*b}_y computed instead of u^{*b}_y alone, for speed

    if (jk) Wrustby[1] = cross(MZdWrKr, uddoty.M) / 2  // fix: original-sample ust is non-jk'd uddot; "/ 2" for consistency with Rademacher-halving
    Wystb = (cross(halfWrKr,Y) - Wrustby[1]) :+ Wrustby  // Wr ' (*un*-FWL'd endog vars); "half" compensates for Rademacher-halving trick
    if (fuzzy) {
      Wrustbt = cross(panelsum(MZdWrKr, uddott.M, info), v)
      if (jk) Wrustbt[1] = cross(MZdWrKr, uddott.M) / 2
      Wtstb = (cross(halfWrKr,T) - Wrustbt[1]) :+ Wrustbt
    }
  } else {
    _jumble(vbc)  // for speed, jumble the old Rademacher weights instead of regenerating from scratch  XXX jumble the smaller matrices instead
    Wrustby = cross(MZdWrKrjk, uddoty.M, vbc)  // Wr'u^{*b}_y computed instead of u^{*b}_y alone, for speed; jk-ing already folded into MZdWrKrjk
    if (jk) Wrustby[1] = cross(MZdWrKr, uddoty.M) / 2
    Wystb = (cross(halfWrKr,Y) - Wrustby[1]) :+ Wrustby  // Wr ' (un-FWL'd endog vars); "half" compensates for Rademacher-halving trick

    if (fuzzy) {
      Wrustbt = cross(MZdWrKrjk, uddott.M, vbc)
      if (jk) Wrustbt[1] = cross(MZdWrKr, uddott.M) / 2
      Wtstb = (cross(halfWrKr,T) - Wrustbt[1]) :+ Wrustbt
    }
  }
  zetahatb = Wystb :/ (fuzzy? Wtstb : WWr/2)  // replication regressions; "/2" compensates for Rademacher-halving trick
  zetahat = zetahatb[1]  // original-sample linear estimate
  return(2 * zetahat - (rowsum(zetahatb)-zetahat)/B1)  // bias-corrected estimate
}

// variance simulation step
void WBSRDD::vs(real colvector Y, | real colvector T) {
  real scalar b; real colvector Yd, Td; real matrix v, _v, ustby, ustbt, ystb, tstb

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  uddoty.M = Yd - Wd * (cross(WdKd,Yd) / WWd)
  jk(uddoty, 1)
  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    uddott.M = Td - Wd * (cross(WdKd,Td) / WWd)
    jk(uddott, 1)
  }

  v = 2 * (runiform(N_G, B2 + 1) :< .5) :- 1  // Rademacher weights for *all* replications
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
    if (jk) ustby[,1] = uddoty.M  // fix: original-sample u* is non-jk'd uddot
    ystb = (Y - uddoty.M) :+ ustby
    if (fuzzy) {
      ustbt = uddott.M :* v; ustbt = ustbt - ZdinvZZd * cross(ZdKd,ustbt)
      if (jk) ustbt[,1] = uddott.M 
      tstb = (T - uddott.M) :+ ustbt
    }
  }

  zetast = cross(WrKr,Y) :/ (fuzzy? cross(WrKr,T) : WWr)  // zeta*, uncorrected estimate for "original" sample at variance-simulating level
  if (bc) {
    dist = J(B2, 1, 0)
    zetastbc = fuzzy? bc(ystb[,1], tstb[,1]) : bc(ystb[,1])  // bias-corrected estimate in "original" bs sample, zeta* - Δ*ᵢ
    for (b=2; b<=B2+1; b++) {
//       printf("."); if (!mod(b-1,50)) printf("\n")
//       displayflush()
      dist[b-1] = (fuzzy? bc(ystb[,b], tstb[,b]) : bc(ystb[,b])) - zetast // deviations of bias-corrected estimates from uncorrected estimate in "original" bs sample zetâ ᵢ - Δ**ᵢ - zeta* (algorithm 3.2, step 3)
    }
//     printf("\n")
  } else {
    zetastbc = zetast
    dist = (cross(ystb, WrKr) :/ (fuzzy? cross(tstb, WrKr) : WWr))[|2\.|] :- zetast
  }
  _sort(dist,1)
}

real colvector WBSRDD::getdist() return(dist)

real scalar WBSRDD::getp(| string scalar ptype) {
  if (ptype=="symmetric" | ptype=="")
    return(colsum(-abs(zetastbc) :> -abs(zetastbc :+ dist)) / B2)
  if (ptype=="equaltail")
    return(2 * min((colsum(dist :< 0) , colsum(dist :> 0))) / B2)
  if (ptype=="lower")
    return(colsum(dist :< 0) / B2)
  if (ptype=="upper")
    return(colsum(dist :> 0) / B2)
  _error(198, `"p type must be "symmetric", "equaltail", "lower", or "upper"."')
}

real rowvector WBSRDD::getci(real scalar level, | string scalar citype) {
  real scalar halfwidth
  if (citype=="symmetric" | citype=="") {
    halfwidth = abs(dist)[round(level/100 * B2)]
    return((zetastbc-halfwidth, zetastbc+halfwidth))
  }
  if (citype=="equaltail")
    return((zetastbc-dist[round((1-level/100)/2 * B2)], zetastbc+dist[round((1-(1-level/100)/2) * B2)]))
  _error(198, `"CI type must be "symmetric" or "equaltail"."')
}

mata mlib create lrdboottest, dir("`c(sysdir_plus)'l") replace
mata mlib add lrdboottest *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
//
// use "C:\Users\drood\Downloads\tmp.dta", clear
// rdboottest Y X, fuzzy(T) bwselect(cerrd) nojk reps(50) bcreps(50) seed(1231)
// matlist e(dist_wb)