mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

struct smatrix {
  real matrix M
}

class WBSRDD {
  real scalar fuzzy, N, WWd, WWr, N_G, B1, B2, jk, granularjk, m, dirty, bc, zetast, zetastbc, p, q, v_sd, auxwttype, hasclust
  real colvector W, ZWd, ZWr, Wd, WrKr, halfWrKr, MZdWrKr, WdKd, clustid, dist
  pointer(real colvector) scalar pMZdWrKrjk
  real matrix info, Zr, Zd, invZZd, invZZr, ZdinvZZd, ZdKd, vbc, vvs
  struct smatrix colvector invMg, Xg, XXg, XinvHg, uddoty, uddott

  void Prep(), vs()
  real colvector getdist()
  real rowvector getci()
  real scalar getp()
  private real colvector bc()
  private void new(), jk()
  private real matrix fold(), MakeWildWeights(), count_binary()
}

void WBSRDD::new() {
  v_sd = 1
}

real matrix WBSRDD::fold(matrix X) return(uppertriangle(X) + lowertriangle(X,0)')  // fold matrix diagonally; returns same values as a quad form, but runs faster because of all the 0's

real colvector triangularkernel  (real scalar bw, real colvector X) return(1:-abs(X)/bw)
real colvector epanechnikovkernel(real scalar bw, real colvector X) return((1 :- (X/bw):^2 / 5) * .75 / sqrt(5))
real colvector uniformkernel     (real scalar bw, real colvector X) return(J(rows(X), 1, 1/bw))

// Return matrix that counts from 0 to 2^N-1 in binary, one column for each number, one row for each binary digit
// except use provided lo and hi values for 0 and 1
real matrix WBSRDD::count_binary(real scalar N, real scalar lo, real scalar hi) {
  real matrix tmp
  if (N<=1) return (lo , hi)
  tmp = count_binary(N-1, lo, hi)
  return (J(1, cols(tmp), lo), J(1, cols(tmp), hi) \ tmp, tmp)
}


// draw wild weight matrix
real matrix WBSRDD::MakeWildWeights(real scalar B, real scalar scaletrickok) {
  real matrix v

  if ((B & auxwttype==0 & N_G*ln(2) < ln(B)+1e-6))
    v = J(N_G,1,1), count_binary(N_G, -m, m)  // complete Rademacher set
  else if (auxwttype==3)
    v = rnormal(N_G, B+1, 0, m)  // normal weights
  else if (auxwttype==4)
    v = m * (rgamma(N_G, B+1, 4, .5) :- 2)  // Gamma weights
  else if (auxwttype==2) {
    v = sqrt(ceil(runiform(N_G, B+1) * 3)) :* ((runiform(N_G, B+1):>=.5):-.5)  // Webb weights, divided by sqrt(2)
    v_sd = 1.6a09e667f3bcdX-001 /*sqrt(.5)*/
  } else if (auxwttype) {
    v = ( rdiscrete(N_G, B+1, 1.727c9716ffb76X-001 \ 1.1b06d1d200914X-002 /*.5+sqrt(.05) \ .5-sqrt(.05)*/) :- 1.5 ) :+ 1.c9f25c5bfedd9X-003 /*.5/sqrt(5)*/  // Mammen weights, divided by sqrt(5)
    v_sd = 1.c9f25c5bfedd9X-002 /*sqrt(.2)*/
  } else {
    v = (runiform(N_G, B+1) :>= .5) :- .5   // Rademacher weights, divided by 2
    v_sd = .5
  }

  if (scaletrickok)
    v[,1] = J(N_G, 1, v_sd)  // keep original residuals in 1 entry to compute base model stat
  else {
    if (v_sd != 1) v = v / v_sd
    v[,1] = J(N_G, 1, 1)
  }
  return(v)
}

// one-time stuff that only depends on exog vars and cluster and kernel definitions
void WBSRDD::Prep(real scalar p, real scalar q, real scalar B1, real scalar B2, string scalar auxwttype, real colvector clustid, real colvector X, real matrix Z, real colvector wt, real scalar h_l, real scalar h_r, real scalar b_l, real scalar b_r, string scalar kernel, 
                  real scalar fuzzy, real scalar bc, real scalar jk) {
  real colvector tmp, Kd, Kr; real matrix ZZd, H, invH, neginvH, Z_W, Xp; real scalar g, kZ, invm; pointer(real colvector function) kernelfn

  this.B1 = B1
  this.B2 = B2
  this.fuzzy = fuzzy
  this.jk = jk
  this.bc = bc
  this.clustid = clustid

  W = X:>=0  // RHS var of interest: ITT

  kernelfn = kernel=="Triangular"? &triangularkernel() : (kernel=="Epanechnikov"? &epanechnikovkernel() : &uniformkernel())
  Kd = (*kernelfn)(b_l, X) :* (X:<0) :* (X:>-b_l) + (*kernelfn)(b_r, X) :* (X:>=0) :* (X:<b_r)
  Kr = (*kernelfn)(h_l, X) :* (X:<0) :* (X:>-h_l) + (*kernelfn)(h_r, X) :* (X:>=0) :* (X:<h_r)

  if (rows(wt)) {
    Kd = Kd :* wt
    Kr = Kr :* wt
  }

  N = rows(X)

  uddoty = uddott = smatrix(1 + jk)  // for jackknife, need jk'd residuals but also non-jk'd residuals for original test stat

  if (hasclust = rows(clustid)) {
    info = panelsetup(clustid,1)
    N_G = rows(info)
  } else {
    info = (1::N), (1::N)  // info = J(0,0,0)
    N_G = N
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

  ZdinvZZd = Zd * invZZd
  MZdWrKr = WrKr - ZdKd * cross(ZdinvZZd, WrKr)

  if (!fuzzy)
    WWr = WWr - ZWr ' invZZr * ZWr
  WWd   = WWd - ZWd ' invZZd * ZWd  // only cross-products FWL'd

  auxwttype = strlower(auxwttype)
  if (.==(this.auxwttype = auxwttype=="rademacher" ? 0 : (auxwttype=="mammen" ? 1 : (auxwttype=="webb" ? 2 : (auxwttype=="normal" ? 3 : (auxwttype=="gamma" ? 4 : .))))))
    _error(198, `"Wild type must be "Rademacher", "Mammen", "Webb", "Normal", or "Gamma"."')

  vvs = MakeWildWeights(B2, 0)  // one-time: make aux weights for variance-simulation level and set v_sd for v-scaling trick at bc level
  if (!hasclust) vbc = MakeWildWeights(B1, 1)  // if not clustering, generate bc weights once and jumble them indirectly by instead jumbling what they multiply against
  halfWrKr = v_sd != 1? WrKr * v_sd : WrKr

  // jk prep
  if (jk) {
    Z_W = sqrt(Kd) :* ((Zd, Wd))  // all RHS vars in DGP regressions, sqrt weights folded in
    kZ = cols(Zd)
    
    uddoty[2].M = J(N,1,0)  // will hold jk'd residuals
    if (fuzzy) uddott[2].M = J(N,1,0)

    invH = invsym(H = (ZZd, ZWd \ ZWd', WWd))
    
    if (hasclust) {
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
      pMZdWrKrjk = &(MZdWrKr :* invMg.M)  // trick: when no clustering, jk hat matrix is diagonal; multiply it one time against unchanging factor
    }
  } else
    pMZdWrKrjk = &MZdWrKr
}

// jk-transform 1st entry in a 2-vector of residuals into the 2nd
// vs = 1 means variance-simulation stage & no trick to pre-compute jk-ing, so do it here
void WBSRDD::jk(struct smatrix rowvector uddot, real scalar vs) {
  real colvector uddotg, S; real scalar g

  if (jk)
    if (hasclust) {
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
real colvector WBSRDD::bc(real scalar b, real colvector Y, | real colvector T) {
  real scalar zetahat; real colvector Yd, Td, u; real rowvector zetahatb, Wrustby, Wystb, Wrustbt, Wtstb; pointer(real matrix) scalar pMZdWrKrjku, pMZdWrKru

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  uddoty.M = Yd - Wd * (cross(WdKd,Yd) / WWd)
  jk(uddoty, 0)

  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    uddott.M = Td - Wd * cross(WdKd,Td) / WWd
    jk(uddott, 0)
  }

  if (hasclust) {  // boottest trick
    vbc = MakeWildWeights(B1, 1)  // Rademacher weights for *all* replications; halved for speed of generation

    Wrustby = cross(panelsum(MZdWrKr, uddoty[1+jk].M, info), vbc)  // Wr'u^{*b}_y computed instead of u^{*b}_y alone, for speed

    if (jk) Wrustby[1] = cross(MZdWrKr, uddoty.M) * v_sd  // fix: original-sample ust is non-jk'd uddot; "* v_sd" for consistency with weight-scaling trick
    Wystb = (cross(halfWrKr,Y) - Wrustby[1]) :+ Wrustby  // Wr ' (*un*-FWL'd endog vars); "half" compensates for Rademacher-halving trick
    if (fuzzy) {
      Wrustbt = cross(panelsum(MZdWrKr, uddott.M, info), vbc)
      if (jk) Wrustbt[1] = cross(MZdWrKr, uddott.M) * v_sd
      Wtstb = (cross(halfWrKr,T) - Wrustbt[1]) :+ Wrustbt
    }
  } else {
    if (b == 1) {
      pMZdWrKrjku = pMZdWrKrjk
      pMZdWrKru = &MZdWrKr
    } else {
      u = unorder(N)  // for additional speed, instead jumble the matrices v is multiplied against
      _collate(uddoty.M, u)
      pMZdWrKrjku = &(*pMZdWrKrjk)[u,]
      if (jk) pMZdWrKru = &MZdWrKr[u,]
    }

    Wrustby = cross(*pMZdWrKrjku, uddoty.M, vbc)  // Wr'u^{*b}_y computed instead of u^{*b}_y alone, for speed; jk-ing already folded into MZdWrKrjk
    if (jk) Wrustby[1] = cross(*pMZdWrKru, uddoty.M) * v_sd
    Wystb = (cross(halfWrKr,Y) - Wrustby[1]) :+ Wrustby  // Wr ' (un-FWL'd endog vars); "half" compensates for Rademacher-halving trick

    if (fuzzy) {
      Wrustbt = cross(*pMZdWrKrjku, uddott.M, vbc)
      if (jk) Wrustbt[1] = cross(*pMZdWrKru, uddott.M) * v_sd
      Wtstb = (cross(halfWrKr,T) - Wrustbt[1]) :+ Wrustbt
    }
  }
  zetahatb = Wystb :/ (fuzzy? Wtstb : WWr * v_sd)  // replication regressions; "* v_sd" compensates for Rademacher-halving trick
  zetahat = zetahatb[1]  // original-sample linear estimate
  return(2 * zetahat - (rowsum(zetahatb)-zetahat)/B1)  // bias-corrected estimate
}

// variance simulation step
void WBSRDD::vs(real colvector Y, | real colvector T) {
  real scalar b; real colvector Yd, Td; real matrix _v, ustby, ustbt, ystb, tstb

  Yd = Y - ZdinvZZd * cross(ZdKd, Y)  // FWL
  uddoty.M = Yd - Wd * (cross(WdKd,Yd) / WWd)
  jk(uddoty, 1)
  if (fuzzy) {
    Td = T - ZdinvZZd * cross(ZdKd, T)
    uddott.M = Td - Wd * (cross(WdKd,Td) / WWd)
    jk(uddott, 1)
  }

  if (hasclust) {  // partial use of boottest trick
    _v = vvs[clustid,]
    ustby = uddoty[1+jk].M :* _v - ZdinvZZd * cross(panelsum(uddoty[1+jk].M :* ZdKd, info), vvs)
    if (jk) ustby[,1] = uddoty.M  // fix: original-sample ust is non-jk'd uddot
    ystb = (Y - uddoty.M) :+ ustby
    if (fuzzy) {
      ustbt = uddott[1+jk].M :* _v - ZdinvZZd * cross(panelsum(uddott[1+jk].M :* ZdKd, info), vvs)
      if (jk) ustbt[,1] = uddott.M
      tstb = (T - uddott.M) :+ ustbt
    }
  } else {
    ustby = uddoty[1+jk].M :* vvs; ustby = ustby - ZdinvZZd * cross(ZdKd,ustby)
    if (jk) ustby[,1] = uddoty.M  // fix: original-sample u* is non-jk'd uddot
    ystb = (Y - uddoty.M) :+ ustby
    if (fuzzy) {
      ustbt = uddott.M :* vvs; ustbt = ustbt - ZdinvZZd * cross(ZdKd,ustbt)
      if (jk) ustbt[,1] = uddott.M 
      tstb = (T - uddott.M) :+ ustbt
    }
  }

  zetast = cross(WrKr,Y) :/ (fuzzy? cross(WrKr,T) : WWr)  // zeta*, uncorrected estimate for "original" sample at variance-simulating level
  if (bc) {
    dist = J(B2, 1, 0)
    zetastbc = fuzzy? bc(1, ystb[,1], tstb[,1]) : bc(1, ystb[,1])  // bias-corrected estimate in "original" bs sample, zeta* - Δ*ᵢ
    for (b=2; b<=B2+1; b++) {
//       printf("."); if (!mod(b-1,50)) printf("\n")
//       displayflush()
      dist[b-1] = (fuzzy? bc(b, ystb[,b], tstb[,b]) : bc(b, ystb[,b])) - zetast // deviations of bias-corrected estimates from uncorrected estimate in "original" bs sample zetâ ᵢ - Δ**ᵢ - zeta* (algorithm 3.2, step 3)
    }
//     printf("\n")
  } else {
    zetastbc = zetast
    dist = (cross(ystb, WrKr) :/ (fuzzy? cross(tstb, WrKr) : WWr))[|2\.|] :- zetast
  }
  _sort(dist,1)
}

real colvector WBSRDD::getdist() return(dist :+ zetastbc)

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
