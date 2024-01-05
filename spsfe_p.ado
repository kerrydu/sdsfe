
capture program drop spsfe_p
program define spsfe_p
version 16

syntax newvarname, [xb mu Residuals  u bc jlms omega]
local cost "`e(function)'"
local cost  =cond("`cost'" == "cost", -1,1)
local distribution = e(distribution)

local te `bc' `jlms'
local nopts: word count `xb'  `residuals'  `u' `te' `omega'
    if `nopts' >1 {
        display "{err}only one statistic may be specified"
        exit 498
    }
	if `nopts'==0 local xb xb
	
	local yvar = e(depvar)
	if strpos("`distribution'","half")>0 local mu =0
	else local mu = _b[/mu]
	
    if "`xb'" != "" {
        _predict `typlist' `varlist'  , xb
		local flagy = e(user)
		local flagy = (strpos("`flagy'","sdsf"))
		if `flagy'>0{
			local rho = e(rho)
			qui replace `varlist' =`varlist' + `rho'*Wy_`yvar'
		}
    }
	else{
		tempvar yhat 
        qui _predict double `yhat'  , xb
		local flagy = e(user)
		local flagy = (strpos("`flagy'","spsf"))
		if `flagy'>0{
			local rho = e(rho)
			qui replace `yhat'=`yhat' + `rho'*Wy_`yvar'
		}			
	}
	
	if "`residuals'"!=""{
		qui gen `typlist' `varlist' = `yvar'-`yhat'
		//qui gen `typlist' `varlist' = `rho'*Wy_`yvar'
		label var `varlist' "prediction of residuals"
	}
	
	if "`omega'"!="" | "`u'"!="" | "`te'"!=""{
		local endovars  "`e(endvars)'"
		if "`endovars'"=="" & "`omega'"!=""{
			 display "{err}omega can only be specified in models with endogenous variables"
             exit 498
		}
		tempvar corterm
		qui gen double `corterm' = 0
		foreach v in `endovars'{
			tempvar r`v'
			qui _predict double `r`v''  , xb eq(`v')
			qui replace `r`v'' = (`v' - `r`v'' )
			qui replace `corterm' = `corterm' + `r`v''*_b["/eta_`v'"]
		}
		
	}
	
	if "`u'"!=""| "`te'"!="" | "`omega'"!=""{
		
		tempvar  lnsigmav lnsigmau sigma2 mustar sigma2s
		qui _predict double `lnsigmav', xb eq(#2)
		qui _predict double `lnsigmau', xb eq(#3)
		if "`omega'"!=""{
			qui gen double `varlist' = `yvar'-`yhat' - `corterm'*exp(`lnsigmav')
			exit 
		}
		tempvar omega
		local endovars  "`e(endvars)'"
		if "`endovars'"==""{
			qui gen double `omega' = `yvar'-`yhat' 
		}
		else{
			qui gen double `omega' = `yvar'-`yhat' - `corterm'*exp(`lnsigmav')
		}
		////
		qui replace `omega' = `cost'*`omega'
		
		qui gen double `sigma2' = exp(2*`lnsigmav') + exp(2*`lnsigmau')
		//qui gen double `lambda' = exp(2*`lnsigmau')/ `sigma'
		qui gen double `mustar' = (exp(2*`lnsigmav') *`mu' - exp(2*`lnsigmau')*`omega')/`sigma2'
		qui gen double `sigma2s' = exp(2*`lnsigmav')* exp(2*`lnsigmau')/`sigma2'
		
		if "`u'"!=""{
			qui gen `typlist' `varlist' = `mustar' + sqrt(`sigma2s')*normalden(`mustar'/sqrt(`sigma2s'))/normal(`mustar'/sqrt(`sigma2s'))
			label var `varlist' "Prediction of inefficiency: u"
		}
		if "`jlms'"!=""{
			qui gen `typlist' `varlist' = exp(-`mustar' - sqrt(`sigma2s')*normalden(`mustar'/sqrt(`sigma2s'))/normal(`mustar'/sqrt(`sigma2s')))
			label var `varlist' "Prediction of efficiency: JLMS estimator"
		}
		
		if "`bc'"!=""{ 
			
			qui gen `typlist' `varlist'  = exp(-`mustar'+0.5*`sigma2s')*normal(`mustar'/sqrt(`sigma2s')-sqrt(`sigma2s'))/normal(`mustar'/sqrt(`sigma2s'))
			label var `varlist' "Prediction of efficiency: BC estimator"
		}
		
		
	}
	


end


