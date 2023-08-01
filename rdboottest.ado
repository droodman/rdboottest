cap program drop rdboottest
program define rdboottest, eclass
  version 13

  local cmdline `0'

  syntax varlist [if] [in], [c(real 0) scalepar(real 1) fuzzy(varname) weights(string) seed(string) nojk nobc BCREPs(integer 500) REPs(integer 999) Level(real `c(level)') PType(string) *]
  
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

  local jk = "`jk'"==""
  local bc = "`bc'"==""

  marksample touse
  markout `touse' `weights' `fuzzy'
  
  rdrobust `varlist' `if' `in', c(`c') scalepar(`scalepar') fuzzy(`fuzzy') weights(`weights') `options'
  ereturn local fuzzy `fuzzy'

  preserve
  qui keep if `touse' & `e(runningvar)' > 0`c'-max(e(h_l), e(b_l)) & `e(runningvar)' < 0`c'+max(e(h_r), e(b_r))
  
  if `c' {
    tempvar runningvar
    qui gen double `runningvar' = `e(runningvar)' - `c'
  }
  else local runningvar `e(runningvar)'

  if "`e(clustvar)'" != "" {
    tempname clustid
    egen long `clustid' = group(`e(clustvar)')
    sort `clustid', stable
    local clustidopt st_data(.,"`clustid'")
  }
  else local clustidopt J(0,1,0)
  
  local covsopt = cond("`e(covs)'"=="", "J(`=_N',0,0)", "`st_data(.,"`e(covs)'")'")
  
  local wtopt = cond(`"`weights'"'=="", "J(0,1,0)", "`st_data(.,"`weights'")'")

  if `"`seed'"'!="" {
    set seed `seed'
    ereturn local seed `seed'
  }
  else ereturn local seed `c(seed)'

  mata _rdboottestM = WBSRDD()
  mata _rdboottestM.Prep(`e(p)', `e(q)', `bcreps', `reps', `clustidopt', st_data(.,"`runningvar'"), `covsopt', `wtopt', st_numscalar("e(h_l)"), st_numscalar("e(h_r)"), st_numscalar("e(b_l)"), st_numscalar("e(b_r)"), "triangular", "`fuzzy'"!="", `bc', `jk')
  mata "`fuzzy'"=="" ? _rdboottestM.vs(st_data(.,"`e(depvar)'")) : _rdboottestM.vs(st_data(.,"`e(depvar)'"), st_data(.,"`fuzzy'"))
  restore

  mata st_numscalar("e(tau_bc_wb)", _rdboottestM.zetastbc * `scalepar')
  mata st_numscalar("e(p_wb)", _rdboottestM.getp("`ptype'"))
  mata CI = _rdboottestM.getci(`level', "`ptype'")
  mata st_numscalar("e(ci_l_rb_wb)", CI[1+(`scalepar'<0)] * `scalepar')
  mata st_numscalar("e(ci_r_rb_wb)", CI[2-(`scalepar'<0)] * `scalepar')
  mata st_matrix("e(dist_wb)", _rdboottestM.getdist() * `scalepar')

  di _n as inp "{hline 19}{c TT}{hline 60}"
  di as inp "    Wild bootstrap {c |}" _col(22) as res %7.0g cond(`bc', e(tau_bc_wb), e(tau_bc_wb)) _col(52) %4.3f e(p_wb) _col(60) %8.0g e(ci_l_rb_wb) _col(73) %8.0g e(ci_r_rb_wb)
  di as inp "{hline 19}{c BT}{hline 60}"
  if `bc' di "Bias-corrected. " _c
  di "Bootstrap method of He and Bartalotti (2020)" _n

  ereturn repost, esample(`touse')
  ereturn local cmdline `cmdline'
  ereturn local cmd rdboottest
end
