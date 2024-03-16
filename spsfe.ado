*! version 0.2, 15Mar2024
* by Kerry Du, kerrydu@xmu.edu.cn

* todo  genweff // mata:genweff [I-rhoW]^(-1)
capture program drop spsfe
program define spsfe, eclass sortpreserve
version 16

cap which lxtsfsp.mlib 
if _rc {
    display as red "xtsfsp package is required"
    display as red "Install xtsfsp via: "
	di "  net install xtsfsp,from(https://raw.githubusercontent.com/kerrydu/spxtsfa/main/ado)"
    exit
}

qui mata mata mlib index

spsfe_vparse `0'
local version `r(version)'
if("`version'"!=""){
	dis "The installed version of sdsfe is `version'"
	checkupdate spsfe
	exit
}

if replay() {
	if (`"`e(cmd)'"' != "spsfe") error 301
		Replay `0'
    }
else	Estimate `0'

end

///////////////////////


program Estimate, eclass sortpreserve

syntax varlist,  [Id(varname) Time(varname) ///
				  EXOGVars(varlist) ENDVars(varlist) iv(varlist) LEAVEout(varlist) ///
				  WXMat(string) wxvars(varlist) WYMat(string) WMat(string) NORMalize(string) ///
				  vhet(varlist) uhet(string) /// 
                  NOCONstant COST TRUNCate ///
				  GENWVARS mldisplay(string) NOLOG level(real 95) mex(varlist) ///
				  DELmissing MLPLOT NOGraph ///
				  INItial(name)  MLMODELopt(string)   ///
				  MLSEarch(string) MLMAXopt(string) DELVE CONSTraints(string)  ///
				  lndetfull lndetmc(numlist >0 min=2 max=2) /* undocumented options */ ] 

							  
local cmdline spsfe `0'

if `"`exogvars'"'!=""{
	if `"`iv'"'!=""{
		di as red "Warning: iv() is ignored as exogvars is specified."
	}
	if `"`leaveout'"'!=""{
		di as red "Warning: leaveout() is ignored as exogvars is specified."
	}
	local iv `exogvars'
}

if("`te'"!="") confirm new var `te'
if("`genwvars'"!=""){
	foreach v in `wxvars'{
		confirm new var Wx_`v'
	}
	if ("`wmat'"!="" | "`wymat'"!=""){
		local yvar: word 1 of `varlist'
		confirm new var Wy_`yvar'
	}

}

if("`mex'"!="") & ("`wmat'"=="") & ("`wymat'"==""){
	di as error " mex() should be combined with wymat() or wmat()"
	error 198
}

if("`mex'"!="" &("`wmat'"!=""|"`wymat'"!="") ){
	local missvar: list mex - varlist
	if (`"`missvar'"'!=""){
		di as error " Variables speicifed in mex() are not included in the frontier function"
		error 198
	}
}


if("`wmat'"!="" & "`wymat'"!=""){
	di as error "wmat() can not be combined with wymat()"
	error 198
}
if("`wmat'"!="" & "`wxmat'"!=""){
	di as error "wmat() can not be combined with wxmat()"
	error 198
}

if "`vhet'"!=""{
	nocparse `vhet'
	local vhetvars `r(varlist)'
}
if "`uhet'"!=""{
	nocparse `uhet'
	local uhetvars  `r(varlist)'
}

if (`"`endvars'"'!="") checkendviv `varlist'  `uhetvars', endvars(`endvars') iv(`iv') exogvars(`exogvars')

gettoken yvar xvars: varlist 

local endvars: list uniq endvars 
local checkendv: list yvar & endvars
if "`checkendv'"!=""{
	di as error "dep. variable should not be in endvars()."
	erorr 198
}

local endx: list  endvars - xvars
local endu: list endx - uhetvars  

global ivlist
local iv: list uniq iv 
//local iv: list iv - endvars
local otheriv: list xvars - endvars 
local iv `iv' `otheriv'
local otheriv: list muvars - endvars 
local iv `iv' `otheriv'
local otheriv: list uhetvars - endvars 
local iv `iv' `otheriv'
local otheriv: list vhetvars - endvars 
local iv `iv' `otheriv'

local iv: list uniq iv 
local iv: list iv - leaveout

global ivlist `iv'

if("`wmat'"!=""){
	spsfe00 `0'
	exit
}

if("`wmat'"=="" & "`wymat'"==""){
	spsfe01 `0'
	exit
}

if ("`wxvars'"!="" & "`wxmat'"==""){
	di as error "wxvars() should be combined with wxmat()"
	error 198
}
if ("`wxvars'"=="" & "`wxmat'"!=""){
	di as error "wxmat() should be combined with wxvars()"
	error 198
}


///////////////////////////////////////////////////////////////////

if ("`nolog'"!="") local nolog qui

global end1 
global end2 
global tranparametrs

preserve
marksample touse 

foreach v in `endvars'{
	local etaterm `etaterm' /eta_`v'
}

foreach v in `endvars'{
	local surterm `surterm' (`v': `v'=`iv')
}



markout `touse' `uhetvars' `vhetvars'  `iv' `endvars'


local diopts level(`level') `mldisplay'
mlopts std, `mlmodelopt' `mlmaxopt' `constraints'
local cns constraints(`constraints')

//_fv_check_depvar `yvar'
if ("`initial'"!="" & "`delve'"!=""){
	di "Warning: initial(`initial') overrides delve"
}
if ("`initial'"!="" & "`mlsearch'"!=""){
	di "Warning: initial(`initial') overrides mlsearch(`mlsearch')"
}
if ("`delve'"!="" & "`mlsearch'"!=""){
	di "Warning: delve overrides mlsearch(`mlsearch')"
}


	parsespmat0 `wymat' 
	parsespmat1 `wymat' `r(ldot)' aname(w_ina)
	local nw = r(nw)	




// examine the data and spmatrix
if "`id'"==""{
	tempvar id 
	qui gen int `id'=_n
}
if "`time'"==""{
	tempvar time 
	qui gen int `time'=1
}

    qui keep `varlist' `wxvars' `id' `time' `uhetvars' `touse' `muvars' `vhetvars' `endvars' `iv'
    tempvar order0
    qui gen int `order0' =_n
// sort data	
	qui issorted `time' `id'	
//tempvar time2

	qui distinct2 `time'
	local T = r(ndistinct)	

//	mata: marksuse = st_data(.,"`touse'")
    mata: _order_0   = st_data(.,"`order0'","`touse'") // record the row# in the original data	

	if ( `nw'!=1 & `nw'!=`T') {
		di as error "Spatial weight matrixs in wymat() are specified as time-varying, but # of spmatrix != # of periods"
		exit 198
	}    
	tempvar time2
	qui egen `time2' = group(`time')
	//global paneltvar `time2'
	//mata mata describe

	qui count if `touse'==0
	local nummissing = r(N)
	if(`nummissing'>0){
		mata: marksuse = st_data(.,"`touse'")
	}

	local Nobs = _N

	//if(`nw'==1 & (`nummissing'>0 | mod(`Nobs',`T')!=0 )){
	//	repeatspw w_ina `T'
	//}


	checkspmat w_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')

    scalar rmin = max(-0.9999,r(min_w_ina))
	scalar rmax = min(0.9999,r(max_w_ina))
	global rmin = rmin
	global rmax = rmax

    qui keep if `touse'
    mata: _pan_tvar =st_data( .,"`time2'")

	if ("`lndetfull'"!=""){
		local bp bp 
		mata: _rho_lndet_ = panlndetfull(w_ina,$rmin,$rmax,`T')
	}
	if ("`lndetmc'"!=""){
		local bp bp 
		tokenize `lndetmc'
		mata: _rho_lndet_ = panlndetmc(`1',`2',w_ina,$rmin,$rmax,`T')
	}
**************

   * generating Wx
	if(`"`wxvars'"'!=""){
	  	
		  parsespmat0 `wxmat' 
		  parsespmat1 `wxmat' `r(ldot)' aname(wx_ina)
		  local nw = r(nw) 
			if ( `nw'!=1 & `nw'!=`T') {
				di as error "Spatial weight matrixs in wxmat() are specified as time-varying, but # of spmatrix != # of periods"
				exit 198
			} 

		  checkspmat wx_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')
		  qui genwvars `wxvars', aname(wx_ina) tvar(`time2') pref(Wx_)
		  local wxvars2  `r(wxnames)'
		  mata: _order_wx = st_data(.,"`wxvars2'","`touse'")

	}
 qui genwvars `yvar', aname(w_ina) tvar(`time2') pref(__W_y_)   
 mata: __wy__ = st_data(.,"__W_y_`yvar'","`touse'")	
************

    //di "`cost'"
	if ("`cost'"!=""){
		mata: _cost = -1
	} 
	else{
		mata: _cost = 1
	}	

************

if "`truncate'"==""{
	local dist h 
	local mu 
}
else{
	local dist t 
	local mu /mu 
}

if `"`endvars'"'!=""{
	local epost _end	
	local endivs `","`endvars'","`iv'" "'
	local nedx:  word count `endvars'
	mata: nedx = `nedx'
	global end1 "Endogeneous variables: `endvar'"
	global end2 "Instrumental variables: `iv'"
}

	if("`initial'"=="" & "`delve'"!="") { 
		//qui sfpanel `yvar' `xvars' `wxvars2',`noconstant' m(bc95) usigma(`uhet') vsigma(`vhet') emean(`wmuvars2' `mu') iterate(50) `cns'
		ml model lf sfscal`dist'() (`yvar'=  `wxvars2' `xvars', `noconstant') (`vhet') (`uhet') `mu',  `cns'
		qui ml max, iterate(50)
	    mat b0 =e(b)
		qui corr `yvar' __W_y_`yvar'
		mat b0 = b0, `=r(rho)'		
		
		foreach v in `endvars'{
		      mat b0=b0,0.5
		 }		
		
		foreach v in `endvars'{
			qui reg `v' `iv'
			mat b0=b0,e(b)
		}		
		
	}


	
//local modeltype = cond("`wxvars'"=="","y-SAR","yx-SAR")
local title Spatial frontier model(SPSF_scaling)

local lnsigv lnsigv
if (`"`endvars'"'!="") local lnsigv lnsigw
//mata mata describe
	ml model d0 spsfscal`dist'`epost'() ///
            (frontier:`yvar' =  `wxvars2' `xvars',`noconstant') ///
            (`lnsigv': `vhet')  (lnsigu: `uhet') `mu' ///
	        (Wy:) `etaterm' `surterm', nopreserve `cns' `mlmodelopt' title(`title')
	
	
	if("`initial'"=="" & "`delve'"!="") { 
		ml init b0,copy
	}
	if ("`initial'"=="" & "`delve'"=="") `nolog' ml search, `mlsearch'
	if ("`initial'"!="") ml init `initial', copy
	if ("`mlplot'"!=""){
		if "`nograph'"!="" set graphics off
		`nolog' ml plot Wy:_cons
		if "`nograph'"!="" set graphics on
	}

   local mlmaxopt `mlmaxopt' noout difficult
   local mlmaxopt: list uniq mlmaxopt   
   `nolog' ml max, `mlmaxopt' 

  foreach v in `wxvars'{
	local sendoutvar `sendoutvar' Wx_`v'
  }

foreach v in `wuvars'{
	local sendoutvar `sendoutvar' Wu_`v'
  }




   ereturn local cmd spsfe
   ereturn local spatialwvars wy_`yvar' `sendoutvar'
   ereturn local endvars `endvars'
   ereturn local cmdbase ml
   ereturn local cmdline `cmdline'
   ereturn local depvar `yvar'
   ereturn local function = cond("`cost'"!="","cost","production")
   ereturn local distribution=cond("`truncate'"=="","half normal","truncated normal") 
   
   global tranparametrs
   Replay , `diopts'

   //global tranparametrs diparm(Wy, label("rho") prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@))) d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2)))  
   //Replay , `diopts'

     di " "
     di in gre "Note: rho = 1/(rmin+rmin*exp(Wy:_cons))+exp(Wy:_cons)/(rmax+rmax*exp(Wy:_cons)); "
     di in gre "      where rmin and rmax are the minimum and maximum eigenvalues of sp matrix"
	 di " "
     di in gre "   ---convert Wy:_cons to the original form---  "
     di " "
          di in smcl in gr abbrev("variable",12) _col(14) "{c |}" /*
               */ _col(21) "Coef." _col(29) "Std. Err." _col(44) "t" /*
               */ _col(49) "P>|t|" _col(59) "[95% Conf. Interval]"
          di in smcl in gr "{hline 13}{c +}{hline 64}"
		  _diparm Wy, label("rho") prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@))) d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2))

  
   local rho = r(est)
   local rhose = r(se)
   ereturn scalar rho = `rho'
   ereturn scalar rho_se = `rhose' 
   ereturn local predict = "spsfe_p"  
  
  tempvar uuhat
  qui spsfepost `uuhat', yvar(`yvar') d(`distribution') `cost' endvars(`endvars') rho(`rho')
  qui genweff `uuhat', aname(w_ina) tvar(`time2') 

/*
   if(`"`te'"'!=""){
		tempname bml
		mat `bml' = e(b)
		mata: _b_ml = st_matrix("`bml'")	
	    local nx: word count `xvars' `wxvars2'
		local nz: word count `uhet'
		if("`noconstant'"=="") local noconstant constant
		if "`vhet'"=="" local vhet ,
		if "`uhet'"=="" local uhet ,
		//important `endivs' = , "`endvars'","`iv'"
		mata:_te_order=spsfscal`dist'`epost'_te(_b_ml,"`yvar'","`wxvars2' `xvars'","`vhet'","`uhet'","`noconstant'" `endivs')
   }	
*/

/////// total marginal effect, direct and indirect marginal effects

 	    mata: totalemat = J(0,2,.)
		mata: diremat = J(0,2,.)
		mata: indiremat = J(0,2,.)
 if("`mex'"!=""){
		tempname bml V0
		mat `bml' = e(b)
		mat `V0'  = e(V)
		mata: _b_ml = st_matrix("`bml'")
		mata:  rhocon = _b_ml[length(_b_ml)]
		mata: V0 = st_matrix("`V0'")
		local varnames: colnames `bml'
		mata: varnames=tokens(st_local("varnames"))
	foreach v in `mex'{
		mata: k = select(1..length(varnames),varnames:==`"`v'"')
		mata: k1 = k[1]
		mata: k = select(1..length(varnames),varnames:==`"Wx_`v'"')
		mata: st_numscalar("kk",length(k))
		if(kk>0){
			mata: k2 = k[1]	
			mata: k1 = k1, k2
		}
		mata: bx =  _b_ml[k1]
		mata: rhocon = _b_ml[length(_b_ml)]
		mata: k1 = k1, length(_b_ml)
		if ("`wxmat'"==""){
		   mata: toteff = meff_sdsfbc(rhocon, bx, V0[k1,k1],  `T',  w_ina,  ///
                        dire=.,  indire=.)
		}
		else{
		   mata: toteff = meff_sdsfbc(rhocon, bx, V0[k1,k1],  `T',  w_ina,  ///
                        dire=.,  indire=.,wx_ina)			
		}
	
		mata: totalemat = totalemat \ toteff 
		mata: diremat   = diremat \ dire 
		mata: indiremat = indiremat \ indire 
	}
 }
 
 if("`mex'"!="" ){
	display _n(2) in gr "Marginal effects are reported as follows."
	di "Note: The standard errors are estimated with the Delta method."
    mata: totalemat =totalemat, ((totalemat[.,1]):/ totalemat[.,2]),2*normal(-abs((totalemat[.,1]):/totalemat[.,2]))
	mata: st_matrix("totaleff",totalemat)
	foreach v in `mex'{
		local rnames `rnames' `"`v'"'
	}
    mat rownames totaleff = `rnames'
	mat colnames totaleff = "Coeff" "se" "z" "P"	

	
		local r = rowsof(totaleff)-1
		local rf "--"
		forvalues i=1/`r' {
			local rf "`rf'&"
		}
		local rf "`rf'-"
		local cf "&  %10s | %12.4f & %12.4f &  %12.4f  & %12.4f &"
		dis _n in gr "Total marginal effects:"
		matlist totaleff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
	
		mata: diremat =diremat, ((diremat[.,1]):/diremat[.,2]),2*(normal(-abs((diremat[.,1]):/diremat[.,2])))
		mata: st_matrix("directeff",diremat)
		mat rownames directeff = `rnames'
		mat colnames directeff = "Coeff" "se" "z" "P"	

		dis _n in gr "Direct marginal effect:"
		matlist directeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")

		mata: indiremat =indiremat, (indiremat[.,1]):/indiremat[.,2],2*(normal(-abs((indiremat[.,1]):/indiremat[.,2])))
		mata: st_matrix("indirecteff",indiremat)		
		mat rownames indirecteff = `rnames'
		mat colnames indirecteff = "Coeff" "se" "z"	"P"


		dis _n in gr "Indirect marginal effect:"
		matlist indirecteff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
		

		ereturn matrix totalmargins = totaleff
		ereturn matrix directmargins = directeff
		ereturn matrix indirectmargins = indirecteff

 }

  restore
  
	if ("`genwvars'"!=""){
		qui gen double Wy_`yvar'=.
		mata: getdatafmata(__wy__,_order_0,"Wy_`yvar'")
		label var Wy_`yvar'  "Wy*`yvar'"
		qui gen double __u_sp__ = .
		mata: getdatafmata(__u_sp__,_order_0,"__u_sp__")
		qui gen double __dir_u_sp__ = .
		label var __dir_u_sp__ "direct inefficiency component"
		mata: getdatafmata(__dir_u_sp__,_order_0,"__dir_u_sp__")
		qui gen double __indir_u_sp__ = __u_sp__ - __dir_u_sp__
		label var _indir_u_sp__ "indirect inefficiency component"
	}

  	if(`"`wxvars'"'!=""&"`genwvars'"!=""){
      foreach v in `wxvars'{
        qui gen double Wx_`v' = .
        label var Wx_`v' `"Wx*`v'"'
        local wxall `wxall' Wx_`v'
      }
	  mata: getdatafmata(_order_wx,_order_0,"`wxall'")
      cap mata mata drop  _order_wx

	}

 	
   	if(`nummissing'>0){
		di "Missing values found"
		di "The regression sample recorded by variable __e_sample__"
		cap drop __e_sample__
		qui cap gen byte __e_sample__ = 0
		label var __e_sample__ "e(sample)"
		mata: getdatafmata(J(length(_order_0),1,1),_order_0,"__e_sample__")
		//cap mata mata drop  _touse		
	}	 	
   cap mata mata drop _order_0
end


/////////////////////////

cap program drop spsfe01
program spsfe01, eclass sortpreserve

syntax varlist,  [Id(varname) Time(varname) ///
				  EXOGVars(varlist) ENDVars(varlist) iv(varlist) LEAVEout(varlist) ///
				  WXMat(string) wxvars(varlist) NORMalize(string) ///
				  vhet(varlist) uhet(string) /// 
                  NOCONstant COST TRUNCate ///
				  GENWVARS mldisplay(string) NOLOG level(real 95) mex(varlist) ///
				  DELmissing MLPLOT NOGraph ///
				  INItial(name)  MLMODELopt(string)   ///
				  MLSEarch(string) MLMAXopt(string) DELVE ///
				  CONSTraints(string)  ///
				  lndetfull lndetmc(numlist >0 min=2 max=2) /* undocumented options */ ] 				

if ("`wxvars'"!="" & "`wxmat'"==""){
	di as error "wxvars() should be combined with wxmat()"
	error 198
}
if ("`wxvars'"=="" & "`wxmat'"!=""){
	di as error "wxmat() should be combined with wxvars()"
	error 198
}


///////////////////////////////////////////////////////////////////
local cmdline spsfe `0'
if ("`nolog'"!="") local nolog qui
global end1 
global end2 
global tranparametrs

preserve
marksample touse 

if "`vhet'"!=""{
	nocparse `vhet'
	local vhetvars `r(varlist)'
}
if "`uhet'"!=""{
	nocparse `uhet'
	local uhetvars  `r(varlist)'
}

local iv $ivlist

gettoken yvar xvars: varlist 

markout `touse' `uhetvars' `vhetvars' `muvars' `iv' `endvars'

local diopts level(`level') `mldisplay'
mlopts std, `mlmodelopt' `mlmaxopt' `constraints'
local cns constraints(`constraints')

//_fv_check_depvar `yvar'
if ("`initial'"!="" & "`delve'"!=""){
	di "Warning: initial(`initial') overrides delve"
}
if ("`initial'"!="" & "`mlsearch'"!=""){
	di "Warning: initial(`initial') overrides mlsearch(`mlsearch')"
}
if ("`delve'"!="" & "`mlsearch'"!=""){
	di "Warning: delve overrides mlsearch(`mlsearch')"
}

if "`id'"==""{
	tempvar id 
	qui gen int `id'=_n
}
if ("`time'"==""){
	tempvar time
	qui gen int `time'=1
}
    qui keep `varlist' `wxvars' `id' `time' `uhetvars' `touse' `muvars' `vhetvars' `iv' `endvars'
    tempvar order0
    qui gen int `order0' =_n
    // sort data	
	qui issorted `time' `id'	
	//tempvar time2
	qui distinct2 `time'
	local T = r(ndistinct)	
//mata mata describe
//	mata: marksuse = st_data(.,"`touse'")
    mata: _order_0   = st_data(.,"`order0'","`touse'") // record the row# in the original data	
   
	tempvar time2
	qui egen `time2' = group(`time')
	//global paneltvar `time2'
	//mata mata describe

	qui count if `touse'==0
	local nummissing = r(N)
	if(`nummissing'>0){
		mata: marksuse = st_data(.,"`touse'")
	}

	local Nobs = _N

    qui keep if `touse'
    mata: _pan_tvar =st_data( .,"`time2'")


**************

   * generating Wx
	if(`"`wxvars'"'!=""){
	  	
		  parsespmat0 `wxmat' 
		  parsespmat1 `wxmat' `r(ldot)' aname(wx_ina)
		  local nw = r(nw) 
			if ( `nw'!=1 & `nw'!=`T') {
				di as error "Spatial weight matrixs in wxmat() are specified as time-varying, but # of spmatrix != # of periods"
				exit 198
			} 

		  checkspmat wx_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')
		  qui genwvars `wxvars', aname(wx_ina) tvar(`time2') pref(Wx_)
		  local wxvars2  `r(wxnames)'
		  mata: _order_wx = st_data(.,"`wxvars2'","`touse'")

	}
   	
************


if "`truncate'"==""{
	local dist h 
	local mu 
}
else{
	local dist t 
	local mu /mu 
}

	foreach v in `endvars'{
		local etaterm `etaterm' /eta_`v'
	}

	foreach v in `endvars'{
		local surterm `surterm' (`v': `v'=`iv')
	}
	if `"`endvars'"'!=""{
		local epost _end	
		local endivs `","`endvars'","`iv'" "'
		local nedx:  word count `endvars'
		mata: nedx = `nedx'
		global end1 "Endogeneous variables: `endvar'"
		global end2 "Instrumental variables: `iv'"		
	}	

	if ("`cost'"!=""){
		mata: _cost = -1
	} 
	else{
		mata: _cost = 1
	}	


	if("`initial'"=="" & "`delve'"!="") { 

		qui frontier `yvar' `wxvars2' `xvars',`noconstant' vhet(`vhet') uhet(`uhet') iterate(50) `cns'  `cost'
	    mat b0 =e(b)
		if `"`dist'"'=="t" mat b0=b0,0 
		
		foreach v in `endvars'{
		      mat b0=b0,0.5
		 }		
		
		foreach v in `endvars'{
			qui reg `v' `iv'
			mat b0=b0,e(b)
		}		
		
		
	}


//local modeltype = cond("`wxvars'"=="","y-SAR","yx-SAR")
local title Stoc. frontier model(SF_scaling)

local lnsigv lnsigv
if (`"`endvars'"'!="") local lnsigv lnsigw
local lfd0 = cond("`epost'"=="","lf","d0")
//mata mata describe
	ml model `lfd0' sfscal`dist'`epost'() (`yvar'= `wxvars2' `xvars', `noconstant') ///
	 (`lnsigv': `vhet') (lnsigu: `uhet') `mu' `etaterm' `surterm', nopreserve `cns' `mlmodelopt' title(`title')

	if("`initial'"=="" & "`delve'"!="") { 
		ml init b0,copy
	}	
	if ("`initial'"=="" ) `nolog' ml search, `mlsearch'
	if ("`initial'"!="") ml init `initial', copy

   local mlmaxopt `mlmaxopt' noout difficult
   local mlmaxopt: list uniq mlmaxopt   
   `nolog' ml max, `mlmaxopt' 

   ereturn local cmd spsfe
   ereturn local endvars `endvars'
   ereturn local cmdbase ml
   ereturn local cmdline `cmdline'
   ereturn local depvar `yvar'
   ereturn local function = cond("`cost'"!="","cost","production")
   ereturn local distribution=cond("`truncate'"=="","half normal","truncated normal")  
   ereturn local predict = "spsfe_p"  
   Replay , `diopts'	

  restore
  
  	if(`"`wxvars'"'!=""&"`genwvars'"!=""){
      foreach v in `wxvars'{
        qui gen double Wx_`v' = .
        label var Wx_`v' `"Wx*`v'"'
        local wxall `wxall' Wx_`v'
      }
	  mata: getdatafmata(_order_wx,_order_0,"`wxall'")
      cap mata mata drop  _order_wx

	} 
 	
   	if(`nummissing'>0){
		di "Missing values found"
		di "The regression sample recorded by variable __e_sample__"
		cap drop __e_sample__
		qui cap gen byte __e_sample__ = 0
		label var __e_sample__ "e(sample)"
		mata: getdatafmata(J(length(_order_0),1,1),_order_0,"__e_sample__")
		//cap mata mata drop  _touse		
	}	 	
   cap mata mata drop _order_0
end



//////////////////////////
cap program drop spsfe00
program spsfe00, eclass sortpreserve

syntax varlist, WMat(string)  [INItial(name) NOCONstant NORMalize(string) ///
                Id(varname)  Time(varname) te(name) GENWVARS mldisplay(string)  ///
                DELmissing MLPLOT NOGraph MLMODELopt(string) level(real 95) COST  ///
				MLSEarch(string) MLMAXopt(string) DELVE CONSTraints(string)  ///
				lndetfull lndetmc(numlist >0 min=2 max=2) NOLOG  TRUNCate ///
                wxvars(varlist)  LEAVEout(varlist) vhet(string) uhet(string) /// 
				ENDVars(varlist) iv(varlist)  mex(varlist) EXOGVars(varlist)] 

local cmdline spsfe `0'
if ("`nolog'"!="") local nolog qui
global end1 
global end2 
global tranparametrs
local iv $ivlist
preserve
marksample touse 

if "`vhet'"!=""{
	nocparse `vhet'
	local vhetvars `r(varlist)'
}
if "`uhet'"!=""{
	nocparse `uhet'
	local uhetvars  `r(varlist)'
}


gettoken yvar xvars: varlist 

markout `touse' `uhetvars' `vhetvars' `muvars' `iv' `endvars'
local diopts level(`level') `mldisplay'
mlopts std, `mlmodelopt' `mlmaxopt' `constraints'
local cns constraints(`constraints') 
//_fv_check_depvar `yvar'
if ("`initial'"!="" & "`delve'"!=""){
	di "Warning: initial(`initial') overrides delve"
}
if ("`initial'"!="" & "`mlsearch'"!=""){
	di "Warning: initial(`initial') overrides mlsearch(`mlsearch')"
}
if ("`delve'"!="" & "`mlsearch'"!=""){
	di "Warning: delve overrides mlsearch(`mlsearch')"
}



parsespmat0 `wmat' 
parsespmat1 `wmat' `r(ldot)' aname(w_ina)
local nw = r(nw)

// 检查权重矩阵与数据是否匹配
if "`id'"==""{
	tempvar id 
	qui gen int `id'=_n
}
if ("`time'"==""){
	tempvar time 
	qui gen int `time'=1
}
    qui keep `varlist' `wxvars' `id' `time' `uhetvars' `touse' `muvars' `vhetvars' `iv' `endvars'
    tempvar order0
    qui gen int `order0' =_n
// sort data	
	qui issorted `time' `id'	
	//tempvar time2

	qui distinct2 `time'
	local T = r(ndistinct)	
//mata mata describe
//	mata: marksuse = st_data(.,"`touse'")
    mata: _order_0   = st_data(.,"`order0'","`touse'") // record the row# in the original data	

	if (`nw'!=1 & `nw'!=`T') {
		di as error "Spatial weight matrixs in wmat() are specified as time-varying, but # of spmatrix != # of periods"
		exit 198
	}    
	tempvar time2
	qui egen `time2' = group(`time')
	//global paneltvar `time2'
	//mata mata describe

	qui count if `touse'==0
	local nummissing = r(N)
	if(`nummissing'>0){
		mata: marksuse = st_data(.,"`touse'")
	}

	local Nobs = _N

	//if(`nw'==1 & (`nummissing'>0 | mod(`Nobs',`T')!=0 )){
	//	repeatspw w_ina `T'
	//}


	checkspmat w_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')

    scalar rmin = max(-0.9999,r(min_w_ina))
	scalar rmax = min(0.9999,r(max_w_ina))
	global rmin = rmin
	global rmax = rmax

    qui keep if `touse'
    mata: _pan_tvar =st_data( .,"`time2'")

	if ("`lndetfull'"!=""){
		local bp bp 
		mata: _rho_lndet_ = panlndetfull(w_ina,$rmin,$rmax,`T')
	}
	if ("`lndetmc'"!=""){
		local bp bp 
		tokenize `lndetmc'
		mata: _rho_lndet_ = panlndetmc(`1',`2',w_ina,$rmin,$rmax,`T')
	}
**************

   * generating Wx
	if(`"`wxvars'"'!=""){
	  if ("`wxmat'"!=""){	  	
		  parsespmat0 `wxmat' 
		  parsespmat1 `wxmat' `r(ldot)' aname(wx_ina)

		  checkspmat wx_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')
		  qui genwvars `wxvars', aname(wx_ina) tvar(`time2') pref(Wx_)
		  local wxvars2  `r(wxnames)'
		  mata: _order_wx = st_data(.,"`wxvars2'","`touse'")
	  }
	  else{
		  genwvars `wxvars', aname(w_ina) tvar(`time2') pref(Wx_)
		  local wxvars2  `r(wxnames)'
		  mata: _order_wx = st_data(.,"`wxvars2'","`touse'")	  	
	  }

	}

qui genwvars `yvar', aname(w_ina) tvar(`time2') pref(__W_y_)
mata: __wy__ = st_data(.,"__W_y_`yvar'","`touse'")	
*****************
	
	
if "`truncate'"==""{
	local dist h 
	local mu 
}
else{
	local dist t 
	local mu /mu 
}

	foreach v in `endvars'{
		local etaterm `etaterm' /eta_`v'
	}

	foreach v in `endvars'{
		local surterm `surterm' (`v': `v'=`iv')
	}
	if `"`endvars'"'!=""{
		local epost _end	
		local endivs `","`endvars'","`iv'" "'
		local nedx:  word count `endvars'
		mata: nedx = `nedx'
		global end1 "Endogeneous variables: `endvar'"
		global end2 "Instrumental variables: `iv'"		
	}	
    	
************

	if ("`cost'"!=""){
		mata: _cost = -1
	} 
	else{
		mata: _cost = 1
	}
	
	
	
	if("`initial'"=="" & "`delve'"!="") { 
		ml model lf sfscal`dist'() (`yvar'=  `wxvars2' `xvars', `noconstant') (`vhet') (`uhet') `mu',  `cns'
		qui ml search
		qui ml max, iterate(50) difficult
	    mat b0 =e(b)
		qui corr `yvar' __W_y_`yvar'
		mat b0 = b0, `=r(rho)'
		
		foreach v in `endvars'{
		      mat b0=b0,0.5
		 }		
		
		foreach v in `endvars'{
			qui reg `v' `iv'
			mat b0=b0,e(b)
		}		
		
	}


	
//local modeltype = cond("`wxvars'"=="","y-SAR","yx-SAR")
local title Spatial frontier model(SPSF_scaling)
local lnsigv lnsigv
if (`"`endvars'"'!="") local lnsigv lnsigw
//mata mata describe
	ml model d0 spsfscal`dist'`epost'() ///
            (frontier:`yvar' =  `wxvars2' `xvars',`noconstant') ///
            (`lnsigv': `vhet')  (lnsigu: `uhet') `mu' ///
	        (Wy:) `etaterm' `surterm', nopreserve `cns' `mlmodelopt' title(`title')
	
	
	if("`initial'"=="" & "`delve'"!="") { 
		ml init b0,copy
	}
	if ("`initial'"=="" & "`delve'"=="") `nolog' ml search, `mlsearch'
	if ("`initial'"!="") ml init `initial', copy
	if ("`mlplot'"!=""){
		if "`nograph'"!="" set graphics off
		`nolog' ml plot Wy:_cons
		if "`nograph'"!="" set graphics on
	}

   local mlmaxopt `mlmaxopt' noout difficult
   local mlmaxopt: list uniq mlmaxopt   
   `nolog' ml max, `mlmaxopt' 

   ereturn local spatialwvars Wy_`yvar' `sendoutvar'
   ereturn local cmd spsfe
   ereturn local cmdbase ml
   ereturn local cmdline `cmdline'
   ereturn local endvars `endvars'
   ereturn local depvar `yvar'
   ereturn local function = cond("`cost'"!="","cost","production")
   ereturn local distribution=cond("`truncate'"=="","half normal","truncated normal") 

   global tranparametrs
   Replay , `diopts'

     di " "
     di in gre "Note: rho = 1/(rmin+rmin*exp(Wy:_cons))+exp(Wy:_cons)/(rmax+rmax*exp(Wy:_cons)); "
     di in gre "      where rmin and rmax are the minimum and maximum eigenvalues of sp matrix"
	 di " "
     di in gre "   ---convert Wy:_cons to the original form---  "
     di " "
          di in smcl in gr abbrev("variable",12) _col(14) "{c |}" /*
               */ _col(21) "Coef." _col(29) "Std. Err." _col(44) "t" /*
               */ _col(49) "P>|t|" _col(59) "[95% Conf. Interval]"
          di in smcl in gr "{hline 13}{c +}{hline 64}"
		  _diparm Wy, label("rho") prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@))) d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2))

   local rho = r(est)
   local rhose = r(se)
   ereturn scalar rho = `rho'
   ereturn scalar rho_se = `rhose' 
   ereturn local predict = "spsfe_p"


  tempvar uuhat
  qui spsfepost `uuhat', yvar(`yvar') d(`distribution') `cost' endvars(`endvars') rho(`rho')
  qui genweff `uuhat', aname(w_ina) tvar(`time2') 
 /////// total marginal effect, direct and indirect marginal effects

 	    mata: totalemat = J(0,2,.)
		mata: diremat = J(0,2,.)
		mata: indiremat = J(0,2,.)
 if("`mex'"!=""){
		tempname bml V0
		mat `bml' = e(b)
		mat `V0'  = e(V)
		mata: _b_ml = st_matrix("`bml'")
		mata:  rhocon = _b_ml[length(_b_ml)]
		//mata: _b_ml[length(_b_ml)]= $rmin/(1+exp(rhocon))+$rmax*exp(rhocon)/(1+exp(rhocon))
		mata: V0 = st_matrix("`V0'")
		//mata: III = I(length(_b_ml))
		//mata: III[length(_b_ml),length(_b_ml)] = exp(rhocon)*(($rmax-$rmin)/(1+exp(rhocon))^2)
		//mata: V = III
		local varnames: colnames `bml'
		mata: varnames=tokens(st_local("varnames"))
	foreach v in `mex'{
		mata: k = select(1..length(varnames),varnames:==`"`v'"')
		mata: k1 = k[1]
		mata: k = select(1..length(varnames),varnames:==`"Wx_`v'"')
		mata: st_numscalar("kk",length(k))
		if(kk>0){
			mata: k2 = k[1]	
			mata: k1 = k1, k2
		}
		mata: bx =  _b_ml[k1]
		mata: rhocon = _b_ml[length(_b_ml)]
		mata: k1 = k1, length(_b_ml)
		if ("`wxmat'"==""){
		   mata: toteff = meff_sdsfbc(rhocon, bx, V0[k1,k1],  `T',  w_ina,  ///
                        dire=.,  indire=.)
		}
		else{
		   mata: toteff = meff_sdsfbc(rhocon, bx, V0[k1,k1],  `T',  w_ina,  ///
                        dire=.,  indire=.,wx_ina)			
		}
	
		mata: totalemat = totalemat \ toteff 
		mata: diremat   = diremat \ dire 
		mata: indiremat = indiremat \ indire 
	}
 }
 


 if("`mex'"!="" ){
	display _n(2) in gr "Marginal effects are reported as follows."
	di "Note: The standard errors are estimated with the Delta method."
    mata: totalemat =totalemat, ((totalemat[.,1]):/ totalemat[.,2]), 2*normal(-abs((totalemat[.,1]):/ totalemat[.,2]))
	mata: st_matrix("totaleff",totalemat)
	foreach v in `mex'{
		local rnames `rnames' `"`v'"'
	}
    mat rownames totaleff = `rnames'
	mat colnames totaleff = "Coeff" "se" "z" "P"	

	
		local r = rowsof(totaleff)-1
		local rf "--"
		forvalues i=1/`r' {
			local rf "`rf'&"
		}
		local rf "`rf'-"
		local cf "&  %10s | %12.4f & %12.4f &  %12.4f  & %12.4f &"
		dis _n in gr "Total marginal effects:"
		matlist totaleff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
	
		mata: diremat =diremat, ((diremat[.,1]):/diremat[.,2]),2*normal(-abs((diremat[.,1]):/diremat[.,2]))
		mata: st_matrix("directeff",diremat)
		mat rownames directeff = `rnames'
		mat colnames directeff = "Coeff" "se" "z" "P"	

		dis _n in gr "Direct marginal effect:"
		matlist directeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")

		mata: indiremat =indiremat, (indiremat[.,1]):/indiremat[.,2],2*normal(-abs((indiremat[.,1]):/indiremat[.,2]))
		mata: st_matrix("indirecteff",indiremat)		
		mat rownames indirecteff = `rnames'
		mat colnames indirecteff = "Coeff" "se" "z"	"P"


		dis _n in gr "Indirect marginal effect:"
		matlist indirecteff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
		

		ereturn matrix totalmargins = totaleff
		ereturn matrix directmargins = directeff
		ereturn matrix indirectmargins = indirecteff

 }

  restore


	if ("`genwvars'"!=""){
		qui gen double Wy_`yvar'=.
		mata: getdatafmata(__wy__,_order_0,"Wy_`yvar'")
		label var Wy_`yvar'  "Wy*`yvar'"
		qui gen double __u_sp__ = .
		mata: getdatafmata(__u_sp__,_order_0,"__u_sp__")
		qui gen double __dir_u_sp__ = .
		label var __dir_u_sp__ "direct inefficiency component"
		mata: getdatafmata(__dir_u_sp__,_order_0,"__dir_u_sp__")
		qui gen double __indir_u_sp__ = __u_sp__ - __dir_u_sp__
		label var _indir_u_sp__ "indirect inefficiency component"
	}
  
  	if(`"`wxvars'"'!=""&"`genwvars'"!=""){
      foreach v in `wxvars'{
        qui gen double Wx_`v' = .
        label var Wx_`v' `"Wx*`v'"'
        local wxall `wxall' Wx_`v'
      }
	  mata: getdatafmata(_order_wx,_order_0,"`wxall'")
      cap mata mata drop  _order_wx

	}  
 	
   	if(`nummissing'>0){
		di "Missing values found"
		di "The regression sample recorded by variable __e_sample__"
		cap drop __e_sample__
		qui cap gen byte __e_sample__ = 0
		label var __e_sample__ "e(sample)"
		mata: getdatafmata(J(length(_order_0),1,1),_order_0,"__e_sample__")
		//cap mata mata drop  _touse		
	}	 	
   cap mata mata drop _order_0
end


///////////////////////subprograms//////////////////
cap program drop Replay
program Replay
	syntax [, Level(cilevel) * ]
	ml display , level(`level')	`options'  $tranparametrs
	tablenote       
end

cap program drop tablenote 
program define tablenote 
version 16 

if `"$tranparametrs"'!=""{
	   di "Note: rho = 1/(rmin+rmin*exp(Wy:_cons))+exp(Wy:_cons)/(rmax+rmax*exp(Wy:_cons)),"  
	   di "      where rmin and rmax are the minimum and maximum eigenvalues of sp matrix"
}
if "`$end1'"!=""{
	di "$end1"
	di "$end2"
}

end

/////////////////////
//////utility comands and function for spxtsfa////

cap program drop genwvars
program define genwvars,rclass

version 16

syntax varlist, aname(name) [tvar(varname) pref(string)]

if `"`tvar'"'==""{
	tempvar tvar 
	qui gen  byte `tavr'=1
}

if `"`pref'"'==""{
	local pref W_
}

mata: genwvars("`varlist'",`aname',"`tvar'","`pref'")
return local wxnames `wxnames'

end


cap program drop genweff
program define genweff

version 16

syntax varlist, aname(name) rho(numlist) [tvar(varname)  ]

if `"`tvar'"'==""{
	tempvar tvar 
	qui gen  byte `tavr'=1
}

if `"`pref'"'==""{
	local pref W_
}

mata: genweff("`varlist'",`aname',`rho',"`tvar'",__dir_u_sp__=.,__u_sp__=.)

end



capture program drop checkspmat
program define checkspmat, rclass

syntax namelist(name=wnames), time(varname) touse(varname) [DELMissing NORMalize(string)]

//preserve

qui count if `touse'==0
local n0 = r(N)

if `n0'>0 & "`delmissing'"==""{
	di as red "missing values found. use delmissing to remove the units from the spmatrix"
	error 198
}


if `n0'>0 & "`delmissing'"!=""{
	di  "missing values found. The corresponding units are deleted from the spmatrix" _n
}

if "`normalize'"==""{
	local ntype=0
}
else if "`normalize'"=="row"{
	local ntype=1
}
else if "`normalize'"=="col"{
	local ntype=2
}
else if "`normalize'"=="spectral"{
	local ntype=3
}
else if "`normalize'"=="minmax"{
	local ntype=4
}
else{
	di as error "errors in normalize(), one of {row,col,spectral,minmax} should be specified. "
	error 198
}


//mata mata describe
foreach w in `wnames'{
    mata: _checkspmat("`time' `touse'",`w',`ntype')
    return scalar rmin_`w' = r(rmin)
    return scalar rmax_`w' = r(rmax)
}


end


///////////////////////////////////

capture program drop parsespmat0
program define parsespmat0,rclass
syntax namelist(name=wnames),[MATA ARRAY] 
if "`mata'"=="" & "`array'"==""{
    return local ldot=","
}

end


capture program drop parsespmat1
program define parsespmat1, rclass
syntax namelist(name=wnames),aname(name) [MATA ARRAY] 
local nw: word count `wnames'

local i=1
if "`mata'"=="" & "`array'"==""{
	mata: `aname' = asarray_create("real")
	foreach w in `wnames'{
		tempname w`i' w`i'_id
		spmatrix matafromsp `w`i'' `w`i'_id' = `w'
		local matanames `matanames' `w`i''
		mata: asarray(`aname',`i',`w`i'')
		local i=`i'+1
	}
	cap mata mata drop `matanames'
	
}
else if "`mata'"!=""{
	mata: `aname' = asarray_create("real")
	local matanames `wnames'
	local i=1
	foreach w in `matanames'{
		mata: asarray(`aname',`i',`w')
		local i=`i'+1
	}

}
else{
	mata: _temparray = asarray_create("real")
	mata: keys = asarray_keys(`wnames')
    mata: keys = sort(keys,1) // sort w in time order
	mata: st_local("keytypes",eltype(keys))
	if ("`keytypes'"!="real"){
		di as error "keys in array `wnames' is not real"
		exit 198
	}
    mata: st_numscalar("r(nw)",length(keys))
	local nw = r(nw)
    forv j=1/`nw'{
		mata: asarray(_temparray,`j',asarray(`wnames',keys[`j']))
	}
	mata: `aname' = asarray_create("real")
	forv j=1/`nw'{
		mata: asarray(`aname',`j',asarray(_temparray,`j'))
	}
	cap mata mata drop _temparray

 }

return scalar nw=`nw'

end 


//////////////////////////////////////

capture program drop issorted
program define issorted
	syntax	varlist 
	
	local sorted : dis "`:sortedby'"
	if "`sorted'" != "`varlist'" {
	    noi disp `"sort data by `varlist'"'
		noi disp "make sure that each spmatrix is the same order" _n
	    sort `varlist'
	}

end

cap program drop nocparse
program define nocparse,rclass 
version 10
syntax [varlist], [NOConstant] 

return local varlist `varlist'


end

/////

cap program drop spsfe_vparse
program define spsfe_vparse,rclass

syntax [varlist], [VERSION Gitee *]

if("`version'"!=""){
	qui findfile spsfe.ado
	mata: vfile = cat("`r(fn)'")
	mata: vfile = select(vfile,vfile:!="")
	mata: vfile = usubinstr(vfile,char(9)," ",.)
	mata: vfile = select(vfile,!ustrregexm(vfile,"^( )+$"))
	mata: st_local("versionuse",vfile[1])
	local versionuse = ustrregexrf("`versionuse'","^[\D]+","")
	gettoken vers versionuse:versionuse, p(", ")
	local versionuse `vers'	
    return local version `vers'
}
if "`gitee'"!="" global gitee gitee

end

cap program drop checkupdate
program define checkupdate
local url1 https://raw.githubusercontent.com/kerrydu/sdsfe/main/ado
local url2 https://gitee.com/kerrydu/sdsfe/raw/main/ado
if `"$gitee"' != ""{
	cap mata: vfile = cat(`"`url2'/`0'.ado"')
	if _rc exit
}
else{
	cap mata: vfile = cat(`"`url1'/`0'.ado"')
	if _rc{
		cap mata: vfile = cat(`"`url2'/`0'.ado"')
	}
	if _rc exit
}


mata: vfile = select(vfile,vfile:!="")
mata: vfile = usubinstr(vfile,char(9)," ",.)
mata: vfile = select(vfile,!ustrregexm(vfile,"^( )+$"))
mata: st_local("versiongit",vfile[1])
local versiongit = ustrregexrf("`versiongit'","^[\D]+","")
gettoken vers versiongit:versiongit, p(", ")
local versiongit `vers'

qui findfile `0'.ado
mata: vfile = cat("`r(fn)'")
mata: vfile = select(vfile,vfile:!="")
mata: vfile = usubinstr(vfile,char(9)," ",.)
mata: vfile = select(vfile,!ustrregexm(vfile,"^( )+$"))
mata: st_local("versionuse",vfile[1])
local versionuse = ustrregexrf("`versionuse'","^[\D]+","")
gettoken vers versionuse:versionuse, p(", ")
local versionuse `vers'

if("`versionuse'"!="`versiongit'"){
	di "New version available, `versionuse' =>`versiongit'"
	di "It can be updated by:"
	di "  net install `0',from(`url1') replace force"
	di "or,"
	di "  net install `0',from(`url2') replace force"
}


end



cap program drop repeatspw 

program define  repeatspw

args w t 

mata: _repeatspw(`w',`t')

end


cap mata mata drop _repeatspw()
mata:

void function _repeatspw(transmorphic matrix w,real scalar t)
{
	wkeys = asarray_keys(w)
	wvalue = asarray(w,wkeys[1])

	for(i=2;i<=t;i++){
		asarray(w,wkeys[1]+i,wvalue)
	}
}

end


cap program drop checkendviv
program define checkendviv
version 16

syntax varlist, [endvars(varlist) iv(varlist) EXOGVars(varlist)]

if `"`iv'"'=="" & `"`exogvars'"'==""{
	di as error "instrumental variables must be specified in iv() or exogvars()"
	error 198
}

if `"`exogvars'"'!=""{
	local comvar: list exogvars & endvars 
	if `"`comvar'"'!=""{
		di as error "{`comvar'} specified in exogvars() must not be included in endvars()"
		error 198
	}	
	exit
}

local iv: list uniq iv 
local endvars: list uniq endvars 
local comvar: list iv & varlist 

if `"`comvar'"'!=""{
	di as error "{`comvar'} specified as instrumental variables must not be included in the explanatory variables"
	error 198
}

local ne: word count `endvars'
local niv: word count `iv'

if `ne'>`niv'{
	di as error "# of iv < # of endogenous variables"
	error 198	
	
}


end


//////////////////////
capture program drop spsfepost
program define spsfepost
version 16
syntax newvarname, yvar(varname) Distribution(string) COST  rho(numlist) [endvars(varlist)]
local cost  =cond("`cost'" == "cost", -1,1)
if strpos("`distribution'","half")>0 local mu =0
else local mu = _b[/mu]

tempvar yhat 
qui _predict double `yhat'  , xb

tempvar  lnsigmav lnsigmau sigma2 mustar sigma2s omega corterm
qui _predict double `lnsigmav', xb eq(#2)
qui _predict double `lnsigmau', xb eq(#3)

qui replace `yhat' =`yhat' + `rho'*Wy_`yvar'

qui gen double `corterm' = 0
foreach v in `endovars'{
	tempvar r`v'
	qui _predict double `r`v''  , xb eq(`v')
	qui replace `r`v'' = (`v' - `r`v'' )
	qui replace `corterm' = `corterm' + `r`v''*_b["/eta_`v'"]
}


		if "`endovars'"==""{
			qui gen double `omega' = `yvar'-`yhat' 
		}
		else{
			qui gen double `omega' = `yvar'-`yhat' - `corterm'*exp(`lnsigmav')
		}
			
		qui replace `omega' = `cost'*`omega'	
		
		qui gen double `sigma2' = exp(2*`lnsigmav') + exp(2*`lnsigmau')
		qui gen double `mustar' = (exp(2*`lnsigmav') *`mu' - exp(2*`lnsigmau')*`omega')/`sigma2'
		qui gen double `sigma2s' = exp(2*`lnsigmav')* exp(2*`lnsigmau')/`sigma2'
		
		qui gen double `varlist' = `mustar' + sqrt(`sigma2s')*normalden(`mustar'/sqrt(`sigma2s'))/normal(`mustar'/sqrt(`sigma2s'))
	

end


