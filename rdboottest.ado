cap program drop rdboottest
program define rdboottest, eclass
  version 13

  local cmdline `0'

  syntax varlist [if] [in], [c(real 0) scalepar(real 1) Level(real `c(level)') fuzzy(varname) weights(string) covs(string) ///
                             seed(string) JACKknife jk nobc REPs(integer 999) BCREPs(integer 500) WEIGHTtype(string) PType(string) *]
  
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

  if `:word count `weighttype'' > 1 {
		di as err "The {cmd:weight:type} option must be {cmdab:rad:emacher}, {cmdab:mam:men}, {cmdab:nor:mal}, {cmdab:web:b}, or {cmdab:gam:ma}."
		exit 198
	}
	if "`weighttype'"'=="" local weighttype rademacher
	else {
		local 0, `weighttype'
		syntax, [RADemacher MAMmen NORmal WEBb GAMma]
		local weighttype `rademacher'`mammen'`normal'`webb'`gamma'
	}


  local jk = "`jk'`jackknife'"!=""
  local bc = "`bc'"==""

  marksample touse
  markout `touse' `weights' `fuzzy' `covs'
  
  rdrobust `varlist' `if' `in', c(`c') scalepar(`scalepar') fuzzy(`fuzzy') weights(`weights') covs(`covs') level(`level') `options'
  ereturn local fuzzy `fuzzy'

  preserve
  qui keep if `touse' & `e(runningvar)' > `c'-max(e(h_l), e(b_l)) & `e(runningvar)' < `c'+max(e(h_r), e(b_r))
  
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
  
  local covsopt = cond("`covs'"   =="", "J(`=_N',0,0)", "`st_data(.,"`covs'")'")
  local wtopt   = cond("`weights'"=="", "J(0,1,0)"    , "`st_data(.,"`weights'")'")

  if `"`seed'"'!="" {
    set seed `seed'
    ereturn local seed `seed'
  }
  else ereturn local seed `c(seed)'

  mata _rdboottestM = WBSRDD()
  mata _rdboottestM.Prep(`e(p)', `e(q)', `bcreps', `reps', "`weighttype'", `clustidopt', st_data(.,"`runningvar'"), `covsopt', `wtopt', st_numscalar("e(h_l)"), st_numscalar("e(h_r)"), st_numscalar("e(b_l)"), st_numscalar("e(b_r)"), "`e(kernel)'", "`fuzzy'"!="", `bc', `jk')
  mata "`fuzzy'"=="" ? _rdboottestM.vs(st_data(.,"`e(depvar)'")) : _rdboottestM.vs(st_data(.,"`e(depvar)'"), st_data(.,"`fuzzy'"))
  restore

  if `bc' mata st_numscalar("e(tau_bc_wb)", _rdboottestM.zetastbc * `scalepar')
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

  ereturn local weighttype_wb `weighttype'
  ereturn local ptype_wb `ptype'
  ereturn scalar reps = `reps'
  ereturn scalar bcreps = `bcreps'
  ereturn local jk_wb `jk'
  ereturn local bc_wb `jk'
  ereturn repost, esample(`touse')
  ereturn local cmdline `cmdline'
  ereturn local cmd rdboottest
end
