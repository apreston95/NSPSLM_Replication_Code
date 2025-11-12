clear all
set more 1
set matsize 800
set mem 400m
cd "C:\Users\344278\Documents\NSPSLM_replication_package\Appendix\Appendix G"
import excel "C:\Users\344278\Documents\NSPSLM_replication_package\Section 4/Main Datasheet.xlsx", sheet("techdat") firstrow case(lower)
gen qdate = q(1960q1) + _n-1
tsset qdate, q

set scheme s1color

graph set window fontface "Georgia"



drop if _n>276


gen t = _n
gen t2 = t^2
gen t3 = t^3
gen t4 = t^4
gen rstockp = stockp_sh/pgdp


foreach var in rgdp rstockp vacancies unemployment realndcserv cpi {
  
  gen l`var' = ln(`var')

} 

g inflation = ln(gdpdeflator/l.gdpdeflator )


* LET H BE THE HORIZON
gen h = t - 1 


* CREATE VARIABLES FOR ESTIMATING ;

local p = 4  

egen shock_std = std(bp_tfp_news_sr)
replace bp_tfp_news_sr = shock_std

global shock bp_tfp_news_sr
global trendlist t 


foreach var in lrgdp ltfp_util lrealndcserv lunemployment lvacancies  ffr  lrstockp {
 quietly gen b`var' = .
  quietly gen up90`var' = .
  quietly gen lo90`var' = .
 
  format b`var' up90`var' lo90`var' %5.1f


}

foreach var in inflation {
 quietly gen b`var' = .
  quietly gen up90`var' = .
  quietly gen lo90`var' = .
 

}



foreach var in lrgdp lunemployment lvacancies lrstockp ltfp_util inflation lrealndcserv {

 forvalues i = 0/20 {
 
reg F`i'.`var' $shock L(1/`p').ltfp_util L(1/`p').lrgdp  L(1/`p').`var'  $trendlist, robust		
				
  gen b`var'h`i' = -100*_b[$shock]

  gen se`var'h`i' = -100*_se[$shock]
  quietly replace b`var' = b`var'h`i' if h==`i' 
 
  
  quietly replace up90`var' = b`var'h`i' + 1.65*se`var'h`i' if h==`i'
  quietly replace lo90`var' = b`var'h`i' - 1.65*se`var'h`i' if h==`i'
  
  }
  }
  
foreach var in  ffr  {

 forvalues i = 0/20 {
 
reg F`i'.`var' $shock L(1/`p').ltfp_util L(1/`p').lrgdp  L(1/`p').`var'  $trendlist, robust		
				
  gen b`var'h`i' = -1*_b[$shock]

  gen se`var'h`i' = -1*_se[$shock]
  quietly replace b`var' = b`var'h`i' if h==`i' 
 
  
  quietly replace up90`var' = b`var'h`i' + 1.65*se`var'h`i' if h==`i'
  quietly replace lo90`var' = b`var'h`i' - 1.65*se`var'h`i' if h==`i'
  
  }
  }


* LABELLLING VARIABLES FOR GRAPH TITLES
label variable lrgdp "Real GDP"
label variable lrstockp "Stock prices"
label variable ffr "Fed. Funds rate"
label variable inflation "Inflation"
label variable lunemployment "Unemployment"
label variable lvacancies "Vacancies"
label variable ltfp_util "TFP"
label variable lrealndcserv "Consumption"

label variable h "Quarters"

rename h Quarters


* PRODUCING IRFs 
foreach var in  lrstockp  lunemployment lvacancies inflation  ffr  { 
tw (rarea up90`var' lo90`var' Quarters, bcolor(cranberry*.6) clw(medthin medthin) title(`: variable label `var'') legend(off))  (line b`var' Quarters, lw(medthick) lc(cranberry) ytitle("% Deviation") xtitle("Quarters") ) (scatteri 0 0 0 20, recast(line) lp(longdash) lw(medthin) lc(cranberry)) if Quarters<=20 ,saving(`shock'_`var'_LP.gph,replace) 

}



foreach var in lrgdp ltfp_util lrealndcserv  { 
tw (rarea up90`var' lo90`var' Quarters, bcolor(cranberry*.6) clw(medthin medthin) title(`: variable label `var'') legend(off))  (line b`var' Quarters, ysc(r(-1.0 0.5)) ylabel(-1.0[0.5]0.5) lw(medthick) lc(cranberry) ytitle("% Deviation") xtitle("Quarters")) (scatteri 0 0 0 20, recast(line) lp(longdash) lw(medthin) lc(cranberry)) if Quarters<=20 ,saving(`shock'_`var'_LP.gph,replace) 
}

graph combine   `shock'_ltfp_util_LP.gph `shock'_lrgdp_LP.gph  `shock'_lrealndcserv_LP.gph `shock'_lvacancies_LP.gph  `shock'_lunemployment_LP.gph `shock'_ffr_LP.gph  `shock'_inflation_LP.gph , rows(2) saving(Figure3B.gph, replace)
// graph export LP_BS_Main.pdf, replace
graph export FigureA11.png, as(png) width(2000) replace

