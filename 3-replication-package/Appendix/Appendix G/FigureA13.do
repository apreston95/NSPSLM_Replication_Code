clear all
set more 1
set matsize 800
clear all
set more 1
set matsize 800
set mem 400m
cd "C:\Users\344278\Documents\NSPSLM_replication_package\Appendix\Appendix G"
import excel "C:\Users\344278\Documents\NSPSLM_replication_package\Section 4/Main Datasheet.xlsx", sheet("extended_bs") firstrow case(lower)
gen qdate = q(1960q1) + _n-1
tsset qdate, q

set scheme s1color

graph set window fontface "Georgia"

gen t = _n
gen t2 = t^2
gen t3 = t^3
gen t4 = t^4



* LET H BE THE HORIZON
gen h = t - 1 


* CREATE VARIABLES FOR ESTIMATING ;

local p = 4  



global shock bs_tfp_news_extend
global trendlist t


foreach var in log_real_gdp log_tfp_util log_real_consumption log_unemployment log_vacancies  ffr  {
 quietly gen b`var' = .
  quietly gen up90`var' = .
  quietly gen lo90`var' = .
 
  format b`var' up90`var' lo90`var' %5.1f


}

foreach var in inflation  {
 quietly gen b`var' = .
  quietly gen up90`var' = .
  quietly gen lo90`var' = .


}


foreach var in log_real_gdp log_tfp_util log_real_consumption log_unemployment log_vacancies {

 forvalues i = 0/20 {
 
ivreg2 F`i'.`var' (f24.log_tfp_util = $shock ) L(1/`p').log_tfp_util L(1/`p').log_real_gdp  L(1/`p').`var'  $trendlist, robust		
				
  gen b`var'h`i' = -1*_b[f24.log_tfp_util]

  gen se`var'h`i' = -1*_se[f24.log_tfp_util]
  quietly replace b`var' = b`var'h`i' if h==`i' 
 
  
  quietly replace up90`var' = b`var'h`i' + 1.65*se`var'h`i' if h==`i'
  quietly replace lo90`var' = b`var'h`i' - 1.65*se`var'h`i' if h==`i'
  
  }
  }
  
foreach var in ffr inflation  {

 forvalues i = 0/20 {
 
ivreg2 F`i'.`var' (f24.log_tfp_util = $shock ) L(1/`p').log_tfp_util L(1/`p').log_real_gdp  L(1/`p').`var'  $trendlist, robust		
				
  gen b`var'h`i' = -1*_b[f24.log_tfp_util]

  gen se`var'h`i' = -1*_se[f24.log_tfp_util]
  quietly replace b`var' = b`var'h`i' if h==`i' 
 
  
  quietly replace up90`var' = b`var'h`i' + 1.65*se`var'h`i' if h==`i'
  quietly replace lo90`var' = b`var'h`i' - 1.65*se`var'h`i' if h==`i'
  
  }
  }


* LABELLLING VARIABLES FOR GRAPH TITLES
label variable log_real_gdp "Real GDP"
label variable ffr "Fed. Funds rate"
label variable inflation "Inflation"
label variable log_unemployment "Unemployment"
label variable log_vacancies "Vacancies"
label variable log_tfp_util "TFP"
label variable log_real_consumption "Consumption"


label variable h "Quarters"

rename h Quarters


* PRODUCING IRFs 
foreach var in    log_unemployment log_vacancies inflation ffr { 
tw (rarea up90`var' lo90`var' Quarters, bcolor(green*.3) clw(medthin medthin) title(`: variable label `var'') legend(off))  (line b`var' Quarters, lw(medthick) lc(green) ytitle("% Deviation") xtitle("Quarters") ) (scatteri 0 0 0 20, recast(line) lp(longdash) lw(medthin) lc(green)) if Quarters<=20 ,saving(`shock'_`var'_LPIV.gph,replace) 

}


foreach var in log_real_gdp log_tfp_util log_real_consumption  { 
tw (rarea up90`var' lo90`var' Quarters, bcolor(green*.3) clw(medthin medthin) title(`: variable label `var'') legend(off))  (line b`var' Quarters, ysc(r(-4.0 1.0)) ylabel(-4.0[1.0]1.0) lw(medthick) lc(green) ytitle("% Deviation") xtitle("Quarters")) (scatteri 0 0 0 20, recast(line) lp(longdash) lw(medthin) lc(green)) if Quarters<=20 ,saving(`shock'_`var'_LPIV.gph,replace) 
}

graph combine   `shock'_log_tfp_util_LPIV.gph `shock'_log_real_gdp_LPIV.gph  `shock'_log_real_consumption_LPIV.gph `shock'_log_vacancies_LPIV.gph  `shock'_log_unemployment_LPIV.gph `shock'_ffr_LPIV.gph  `shock'_inflation_LPIV.gph , rows(2) saving(Figure3B.gph, replace)

graph export FigureA13.png, as(png) width(2000) replace
