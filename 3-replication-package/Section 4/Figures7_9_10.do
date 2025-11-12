clear all
set more 1
set matsize 800
set mem 400m
cd "C:\Users\344278\Documents\NSPSLM_replication_package\Section 4"
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


foreach var in log_real_gdp log_tfp_util log_real_consumption log_unemployment log_vacancies  ffr pvs {
 quietly gen b`var' = .
  quietly gen up90`var' = .
  quietly gen lo90`var' = .
 
  format b`var' up90`var' lo90`var' %5.1f


}

foreach var in inflation jl_rate jf_rate {
 quietly gen b`var' = .
  quietly gen up90`var' = .
  quietly gen lo90`var' = .


}

foreach var in durables_mich unemp_index_mich {
  egen std_`var' = std(`var')	
  replace `var' = std_`var'
  quietly gen b`var' = .
  quietly gen up90`var' = .
  quietly gen lo90`var' = .	
	
	format b`var' up90`var' lo90`var' %5.1f
}



foreach var in log_real_gdp log_tfp_util log_real_consumption log_unemployment log_vacancies  pvs jl_rate jf_rate {

 forvalues i = 0/20 {
 
reg F`i'.`var' $shock L(1/`p').log_tfp_util L(1/`p').log_real_gdp  L(1/`p').`var'  $trendlist, robust		
				
  gen b`var'h`i' = -100*_b[$shock]

  gen se`var'h`i' = -100*_se[$shock]
  quietly replace b`var' = b`var'h`i' if h==`i' 
 
  
  quietly replace up90`var' = b`var'h`i' + 1.65*se`var'h`i' if h==`i'
  quietly replace lo90`var' = b`var'h`i' - 1.65*se`var'h`i' if h==`i'
  
  }
  }
  
foreach var in ffr inflation durables_mich unemp_index_mich {

 forvalues i = 0/20 {
 
reg F`i'.`var' $shock L(1/`p').log_tfp_util L(1/`p').log_real_gdp  L(1/`p').`var'  $trendlist, robust		
				
  gen b`var'h`i' = -1*_b[$shock]

  gen se`var'h`i' = -1*_se[$shock]
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
label variable pvs "PVS"
label variable durables_mich "Durable purchase sentiment"
label variable unemp_index_mich "Unemployment expectations index"
label variable jl_rate "Job loss rate"
label variable jf_rate "Job finding rate"

label variable h "Quarters"

rename h Quarters


* PRODUCING IRFs 
foreach var in    log_unemployment log_vacancies inflation ffr pvs durables_mich unemp_index_mich jl_rate jf_rate { 
tw (rarea up90`var' lo90`var' Quarters, bcolor(cranberry*.3) clw(medthin medthin) title(`: variable label `var'') legend(off))  (line b`var' Quarters, lw(medthick) lc(cranberry) ytitle("% Deviation") xtitle("Quarters") ) (scatteri 0 0 0 20, recast(line) lp(longdash) lw(medthin) lc(cranberry)) if Quarters<=20 ,saving(`shock'_`var'_LP.gph,replace) 

}


foreach var in log_real_gdp log_tfp_util log_real_consumption  { 
tw (rarea up90`var' lo90`var' Quarters, bcolor(cranberry*.3) clw(medthin medthin) title(`: variable label `var'') legend(off))  (line b`var' Quarters, ysc(r(-1.0 0.5)) ylabel(-1.0[0.5]0.5) lw(medthick) lc(cranberry) ytitle("% Deviation") xtitle("Quarters")) (scatteri 0 0 0 20, recast(line) lp(longdash) lw(medthin) lc(cranberry)) if Quarters<=20 ,saving(`shock'_`var'_LP.gph,replace) 
}

graph combine   `shock'_log_tfp_util_LP.gph `shock'_log_real_gdp_LP.gph  `shock'_log_real_consumption_LP.gph `shock'_log_vacancies_LP.gph  `shock'_log_unemployment_LP.gph `shock'_ffr_LP.gph  `shock'_inflation_LP.gph , rows(2) saving(Figure3B.gph, replace)

graph export Figure7.png, as(png) width(2000) replace

graph combine   `shock'_durables_mich_LP.gph `shock'_unemp_index_mich_LP.gph, rows(1) saving(Figure3B.gph, replace)
graph export Figure9.pdf, replace


graph combine   `shock'_pvs_LP.gph , rows(1) saving(Figure3B.gph, replace)
graph export Figure10.pdf, replace


