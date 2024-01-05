{smcl}
{cmd:help spsfe}
{hline}

{title:Title}

{p2colset 5 13 15 2}{...}
{p2col :{hi:spsfe} {hline 2}}Spatial stochastic frontier models with endogeneity in the style of {help spsfe##Kutlu2020:{bind:Kutlu et al. (2020)}} {p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Estimation syntax

{p 8 17 2}
{cmd:spsfe} {depvar} [{indepvars}] [{cmd:,} {it:options}]

{pstd}
Version syntax

{p 8 17 2}
{cmd:spsfe}{cmd:,} {opt ver:sion}

{pstd}
Replay syntax

{p 8 17 2}
{cmd:spsfe} [{cmd:,} {cmdab:l:evel(}{help level##remarks:{it:#}}{cmd:)}]


{synoptset 31 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Data structure}

{phang}
{cmd: id({it:varname})} specifies cross-sectional id variable. 

{phang}
{cmd: time({it:varname})} specifies time variable. It must be specified for panel data. If not, the data is assumed to be cross-sectional. 

{syntab :Frontier}
{synopt :{opt nocons:tant}}suppress constant term{p_end}
{synopt :{opt cost}}fit cost frontier model; default is {cmd:production}{p_end}
{synopt :{cmd:wxvars({it:varlist})}}spatially lagged independent variables 
in the frontier function{p_end}

{syntab : Ancillary equations}
{synopt :{cmd:uhet(}{it:varlist} [{cmd:,} {opt nocons:tant}]{cmd:)}}
specify explanatory
variables for the inefficiency variance function; 
use {opt noconstant}
to suppress constant term{p_end}

{synopt :{cmd:vhet({it:varlist})}}
specify explanatory
variables for the idiosyncratic error variance function.{p_end}

{syntab :Spatial setting}
{synopt :{cmdab:wm:at(}{it:W1 [W2...WT][,mata array]}{cmd:)}}specify spatial weight matrix{p_end}

{synopt :{cmdab:wym:at(}{it:W1 [W2...WT][,mata array]}{cmd:)}}specify spatial weight matrix for 
spatial lag term{p_end}

{synopt :{cmdab:wxm:at(}{it:W1 [W2...WT][,mata array]}{cmd:)}}specify spatial weight matrix for 
spatial Durbin terms in the frontier{p_end}


{syntab :Regression}
{synopt :{cmdab:init:ial(}{it:{help matrix:matname}}{cmd:)}}specify initial values matrix{p_end}
{synopt :{cmdab:endv:ars(}{it:varlist}{cmd:)}}specify endogeneous variables{p_end}
{synopt :{cmd:iv(}{it:varlist}{cmd:)}}specify instrumental variables{p_end}
{synopt :{cmdab:leaveo:ut(}{it:varlist}{cmd:)}}specify included exogenous variables to be left out{p_end}
{synopt :{cmd:mlsearch(}{it:{help ml##model_options:search_options}}{cmd:)}}specify options for searching initial values{p_end}
{synopt :{opt delve}}delve into maximization problem to find initial values{p_end}
{synopt :{opt mlplot}}use ml plot to find better initial values{p_end}
{synopt :{cmd:mlmodel(}{it:{help ml##model_options:model_options}}{cmd:)}}control {cmd:ml model} options{p_end}
{synopt :{cmd:mlmax({it:{help ml##ml_maximize_options:maximize_options}})}}control {cmd:ml maximize} options{p_end}
{synopt :{opt delmissing}}delete the units with missing observations from spmatrix{p_end}

{syntab :Reporting}
{synopt :{cmd:nolog}}omit the display of the criterion function iteration log{p_end}
{synopt :{cmd:mex({it:varlist})}}reports marginal effects of variables in the frontier{p_end}
{synopt :{cmdab:mldis:play(}{it:{help ml##display_options:display_options}}{cmd:)}}control {cmd:ml display} options; seldom used{p_end}
{synopt :{cmd:te(}{it:{help newvar:effvar}}{cmd:)}}create efficiency variables{p_end}
{synopt :{opt genwxvars}}generate the spatial Durbin terms in the frontier function{p_end}

{syntab :Other}
{synopt :{cmdab:constraints(}{it:{help estimation options##constraints():constraints}}{cmd:)}}apply specified linear constraints{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}


{marker description}{...}
{title:Description}

{pstd}
{opt spsfe} fits spatial stochastic production or cost frontier models
following the methodology of 
{help spsfe##Kutlu2020:{bind:Kutlu et al. (2020)}}. See 
{help spsfe##Kutlu2020:{bind:Kutlu et al. (2020)}} for a detailed
explanation of the methodology and empirical analyses.


{marker options}{...}
{title:Options for the estimation syntax}

{dlgtab:Data structure}

{phang}
{cmd: id({it:varname})} specifies cross-sectional id variable. 

{phang}
{cmd: time({it:varname})} specifies time variable. It must be specified for panel data. If not, the data is assumed to be cross-sectional. 

{dlgtab:Frontier}

{phang}
{cmd:noconstant} suppresses the constant term (intercept) in the frontier function.

{phang}
{cmd:cost} specifies that the model to be fit is a cost frontier model.  The
default is {cmd:production}.

{phang}
{cmd: wxvars({it:varlist})} specifies spatially lagged independent variables in the frontier function.


{dlgtab:Ancillary equations}

{phang}{cmd:uhet(}{it:varlist}[{cmd:,} {cmd:noconstant}]{cmd:)}
specifies that the technical inefficiency component is heteroskedastic,
with the variance expressed as a function of the covariates 
defined in {it:varlist}. Specifying {cmd:noconstant} suppresses 
the constant in this function.

{phang}{cmd:vhet(}{it:varlist}{cmd:)}
specifies that the idiosyncratic error component is heteroskedastic,
with the variance expressed as a function of the covariates defined in
{it:varlist}. 


{dlgtab:Spatial weight matrix}

{phang}
{opt wmat(W1 [W2 ... WT][,mata array])} specifies that the spatial weight matrices. If specified, 
all spatial terms are assumed with the same weight matrices. 

{phang}
By default, the weight matrices are Sp object. mata declares weight matrices are mata matrix. 
If one weight matrix is specified, it assumes time-constant weight matrix. 
For time-varying cases, T weight matrices should be specified in time order. 
Alternatively, using array to declares weight matrices are store in a array.  
If only one matrix is stored in the specified array, the time-constant weight matrix is assumed.  
Otherwise, the keys of the array specifies time information and 
the values store time-specific weight matrices.  


{phang}
{opt wymat(W1 [W2 ... WT][,mata array])} specifies that the spatial weight matrices 
for the Spatial lag of the independent variable. 

{phang}
{opt wxmat(W1 [W2 ... WT][,mata array])} specifies that the spatial weight matrices 
for the Spatial Durbin terms in the frontier function. 



{phang}
{opt normalize(row|col|spectral|minmax)} specifies the normalized method of spatial weight matrixs. 
By default, the command would not normalization the spatial weight matrixs. normalize(row) is row normalisation;
normalize(col) is collumn normalisation; normalize(spectral) is spectral normalisation;
normalize(minmax) is minmax normalisation.


{dlgtab:Regression}

{phang}
{cmd:initial(}{it:{help matrix:matname}}{cmd:)} specifies that {it:matname} is
the initial value matrix.

{phang}
{cmd:endvars(}{it:varlist}{cmd:)} specifies  the variables  are to be treated as endogenous. If this option is not specified, all variables are 
assumed to be exogenous.

{phang}
{cmd:iv(}{it:varlist}{cmd:)} specifies that the variables in ivarlist 
are to be used as instrumental variables to handle endogeneity.

{phang}
{cmd:leaveout(}{it:varlist}{cmd:)} specifies that the variables are to be taken
out of the default iv list.

{phang}
{cmd:mlsearch(}{it:{help ml##model_options:search_options}}{cmd:)} specifies ml search options for searching initial values.

{phang}
{cmd:delve} provides a regression-based methodology to search for 
initial values. The default is to use {helpb ml search:ml search} with default options.

{phang}
{opt mlplot} specifes using ml plot to find better initial values.

{phang}
{cmd:mlmodel({it:{help ml##mlmode:model_options}})} controls the {cmd:ml}
{cmd:model} options; it is seldom used.

{phang}
{cmd:mlmax({it:{help ml##ml_max_descript:maximize_options}})} controls the
{cmd:ml max} options; it is seldom used.

{phang}
{opt delmissing} deletes the units with missing observations from spmatrix.

{dlgtab:Reporting}

{phang}
{cmd:nolog} suppresses the display of the criterion function iteration log.

{phang}
{cmd:mex({it:varlist})} reports total, direct and indirect 
marginal effects of variables in the frontier function.

{phang}
{cmd:meu({it:varlist})} reports total, direct and indirect 
marginal effects of variables in the inefficiency function.

{phang}
{cmd:mldisplay({it:{help ml##mldisplay:display_options}})} controls the
{cmd:ml display} options; it is seldom used.

{phang}
{cmd:te(}{it:{help newvar:effvar}}{cmd:)} generates
the production or cost efficiency variable.

{phang}
{cmd:genwxvars} generates the spatial Durbin terms. 
The option automatically extends any specified variable name in wxvars() 
with Wx_ prefix.


{marker optionsversion}{...}
{title:Options for the version and replay syntax}

{phang}
{cmd:version} displays the version of {cmd:spsfe} installed on Stata and the
program author information.  This option can be used only in version syntax.

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals.  The default is {cmd:level(95)} or as set by {helpb set level}.
This option can be used in the replay syntax or in 
{cmd:mldisplay(}{it:{help ml##mldisplay:display_options}}{cmd:)}.

{marker optionsversion}{...}
{title:Other}

{phang}
{cmdab:constraints(}{it:{help estimation options##constraints():constraints}}{cmd:)}
specifies linear constraints for the estimated model.


    {marker examples}{...}
    {title:Examples}
    
        {title:SP-SF model with time-constant spatial weight matrix}
    
    {pstd}
    Setup{p_end}
    {phang2}{bf:. {stata "mata mata matuse w_ex1,replace"}}{p_end}
    {phang2}{bf:. {stata "use spsfe1.dta"}}{p_end}
    
    {pstd}
    Stochastic Stoch. production model {p_end}
    {phang2}{bf:. {stata "spsfe y x,  uhet(z) noconstant  id(id) time(t)"}}{p_end}
    
    
   
       {title:SP-SF model with endogeneous variables}

    {pstd}
    Setup{p_end}
    {phang2}{bf:. {stata "mata mata matuse w_ex1,replace"}}{p_end}
    {phang2}{bf:. {stata "use spsfe2.dta"}}{p_end}

    {pstd}
    Stochastic Stoch. production model {p_end}
    {phang2}{bf:. {stata "spsfe y x, wmat(wm,mata) wxvars(x) uhet(z) noconstant endvars(x z) iv(q1 q2) id(id) time(t)"}}{p_end}

  

{marker results}{...} 
{title:Stored results}

{pstd}
{cmd:spsfe} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(p)}}significance{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}command used for estimation{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared test{p_end}
{synopt:{cmd:e(vce)}}{it:oim}{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(which)}}{cmd:max} or {cmd:min}; whether optimizer is to perform maximization or minimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(cmdbase)}}base command{p_end}
{synopt:{cmd:e(function)}}production or cost{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker acknowledgments}{...}
{title:Acknowledgments}

{pstd}
Kerui Du acknowledges financial support from the National Natural Science Foundation of China(Grant no. 72074184).



{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
{cmd:spsfe} is not an official Stata command.  It is a third-party command
programmed as a free contribution
to the research society.  By choosing to download, install, and use the
{cmd:spsfe} package, users assume all the liability for any
{cmd:spsfe}-package-related risk.  If you encounter any problems with the
{cmd:spsfe} package, or if you have comments, suggestions, or questions, please
send an email to Kerui Du at 
{browse "mailto:kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}.


{marker citation}{...}
{title:References}

{marker Kutlu2020}{...}
{phang}Kutlu, Levent, Kien C. Tran, and Mike G. Tsionas (2020).
 “A Spatial Stochastic Frontier Model with Endogenous Frontier and Environmental Variables.”  European Journal of Operational Research,
286, 1: 389–99. https://doi.org/10.1016/j.ejor.2020.03.020.




{marker author}{...}
{title:Author}

{pstd}
Kerui Du{break}
Xiamen University{break}
School of Management{break}
China{break}
{browse "kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}{break}


{pstd}
Federica Galli{break}
University of Bologna{break}
Department of Statistical Sciences “Paolo Fortunati”{break}
Italy{break}
{browse "federica.galli14@unibo.it":federica.galli14@unibo.it}{break}


{marker see}{...}
{title:Also see}

{p 7 14 2}{manhelp frontier R}, 
{manhelp xtfrontier XT}{p_end} 
{p 7 14 2}{helpb sfpanel}, {helpb sfkk},  {helpb nwxtregress} (if installed)
