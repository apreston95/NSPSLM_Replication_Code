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

*reg log_real_gdp t 
*predict lhs_var, r

//////////////////////////////////////////////////////////////////////////// OUTPUT ///////////////////////////////////////////////////////////////////////


tsfilter hp lhs_var = log_real_gdp , smooth(1600)

g xx = bs_tfp_news_extend

local nlags 4
local irflag = 20

* compute necessary leads and lags
forvalues lag = 1/`nlags' {
    foreach yy of varlist xx lhs_var  {
        cap drop `yy'L`lag'
        g `yy'L`lag' = `yy'[_n-`lag']
    }
}

foreach yy of varlist xx lhs_var {
    forvalues hh = 0/`irflag' {

        cap drop `yy'F`hh'
        g `yy'F`hh' = `yy'[_n+`hh']
    }
}

forvalues hh = 1/`irflag' {
    cap drop xxL`hh'
    g xxL`hh' = xx[_n-`hh']
}

capture program drop lpirf
program lpirf, eclass
    args nlags irflag
    
    * compute IRFs
    local irflagp1 = `irflag'+1
    mat bbb = J(`irflagp1',1,0)

    forvalues h = 0/`irflag' {
        local hp1 = `h'+1
        
        qui reg lhs_varF`h' xx xxL* lhs_varL*, r
        mat bbb[`hp1',1] = _b[xx]
    }
    
    tempname bb
    matrix `bb'=bbb'
    ereturn clear
    ereturn post `bb'
    ereturn local cmd="bootstrap"

end

lpirf `nlags' `irflag'

matrix list e(b)

*bootstrap _b, reps(10) nodrop: lpirf `nlags' `irflag'

capture program drop histdeco
program histdeco, eclass
    args nlags irflag
    
    * compute IRFs
    local irflagp1 = `irflag'+1
    mat bbb = J(`irflagp1',1,0)

    forvalues h = 0/`irflag' {
        local hp1 = `h'+1
        
        qui reg lhs_varF`h' xx xxL* lhs_varL*, r
        mat bbb[`hp1',1] = _b[xx]
    }
    
    
    * compute HD
    cap drop xx_hd
    gen xx_hd = 0
    
    forvalues h = 1/`irflag' {
        qui replace xx_hd = xx_hd + bbb[`h',1]*xxL`h'

    }
		qui replace xx_hd = xx_hd + bbb[1,1]*xx 

    mkmat xx_hd if xx_hd != ., matrix(histdec)

    tempname bb
    matrix `bb'=histdec'
    ereturn clear
    ereturn post `bb'
    ereturn local cmd="bootstrap"

end

lpirf `nlags' `irflag'

matrix list e(b)

histdeco `nlags' `irflag'

summ lhs_var

local min = r(min)
g recession_ind = r(max) if recession == 1 
replace recession_ind = r(min) if recession == 0 


* Overlaid line plot: actual 'lhs_var' vs. the shock's contribution
twoway ( area recession_ind quarter if lhs_var <. & xx_hd<.,    color(gs14)  base(`min')  ) ///
    (line lhs_var quarter if lhs_var <. & xx_hd<., lcolor(blue) lwidth(medthick) ///
        lpattern(solid) ///
        ) ///
    (line xx_hd quarter if lhs_var <. & xx_hd<., lcolor(red) lwidth(medthick) ///
        lpattern(dash)) ///
    , ///
    title("Real GDP") ///
    ytitle("Levels") ///
	xlabel(1960(10)2020) ///
    xtitle("Year", margin(small)) ///
    legend(off)
	
graph save output_hd.gph, replace

graph export "HD_Output.png", as(png) replace



////////////////////////////////////////////////////////////////////////// CONSUMPTION ///////////////////////////////////////////////////////////////////////

keep quarter bs_tfp_news_extend log_tfp_util log_tfp log_real_gdp log_real_consumption log_real_stock_price log_hours log_vacancies log_unemployment ffr inflation recession qdate t

tsfilter hp lhs_var = log_real_consumption , smooth(1600)

g xx = bs_tfp_news_extend

local nlags 4
local irflag = 20

* compute necessary leads and lags
forvalues lag = 1/`nlags' {
    foreach yy of varlist xx lhs_var  {
        cap drop `yy'L`lag'
        g `yy'L`lag' = `yy'[_n-`lag']
    }
}

foreach yy of varlist xx lhs_var {
    forvalues hh = 0/`irflag' {

        cap drop `yy'F`hh'
        g `yy'F`hh' = `yy'[_n+`hh']
    }
}

forvalues hh = 1/`irflag' {
    cap drop xxL`hh'
    g xxL`hh' = xx[_n-`hh']
}



lpirf `nlags' `irflag'

matrix list e(b)

*bootstrap _b, reps(10) nodrop: lpirf `nlags' `irflag'


lpirf `nlags' `irflag'

matrix list e(b)

histdeco `nlags' `irflag'

summ lhs_var

local min = 1.1*r(min)
g recession_ind = r(max) if recession == 1 
replace recession_ind = 1.1*r(min) if recession == 0 


twoway ///
    ( area recession_ind quarter if lhs_var <. & xx_hd<., color(gs14) base(`min') ) ///
    ( line lhs_var quarter if lhs_var <. & xx_hd<., lcolor(blue) lwidth(medthick) lpattern(solid) ) ///
    ( line xx_hd quarter if lhs_var <. & xx_hd<., lcolor(red) lwidth(medthick) lpattern(dash) ) ///
    , ///
    title("Consumption") ///
    legend(off) ///
    ytitle("Levels") ///
    xlabel(1960(10)2020) ///
    xtitle("Year", margin(small)) ///
    yscale(noextend) ///
	ylabel(-0.03(0.01)0.02)
	
graph save consumption_hd.gph, replace

////////////////////////////////////////////////////////////////////////// UNEMPLOYMENT ///////////////////////////////////////////////////////////////////////

keep quarter bs_tfp_news_extend log_tfp_util log_tfp log_real_gdp log_real_consumption log_real_stock_price log_hours log_vacancies log_unemployment ffr inflation recession qdate t

tsfilter hp lhs_var = log_unemployment , smooth(1600)

g xx = bs_tfp_news_extend

local nlags 4
local irflag = 20

* compute necessary leads and lags
forvalues lag = 1/`nlags' {
    foreach yy of varlist xx lhs_var  {
        cap drop `yy'L`lag'
        g `yy'L`lag' = `yy'[_n-`lag']
    }
}

foreach yy of varlist xx lhs_var {
    forvalues hh = 0/`irflag' {

        cap drop `yy'F`hh'
        g `yy'F`hh' = `yy'[_n+`hh']
    }
}

forvalues hh = 1/`irflag' {
    cap drop xxL`hh'
    g xxL`hh' = xx[_n-`hh']
}


lpirf `nlags' `irflag'

matrix list e(b)

*bootstrap _b, reps(10) nodrop: lpirf `nlags' `irflag'


lpirf `nlags' `irflag'

matrix list e(b)

histdeco `nlags' `irflag'

summ lhs_var

local min = r(min)
g recession_ind = r(max) if recession == 1 
replace recession_ind = r(min) if recession == 0 


twoway ///
    ( area recession_ind quarter if lhs_var <. & xx_hd<., color(gs14) base(`min') ) ///
    ( line lhs_var quarter if lhs_var <. & xx_hd<., lcolor(blue) lwidth(medthick) lpattern(solid) ) ///
    ( line xx_hd quarter if lhs_var <. & xx_hd<., lcolor(red) lwidth(medthick) lpattern(dash) ) ///
    , ///
    title("Unemployment") ///
    legend(off) ///
    ytitle("Levels") ///
    xlabel(1960(10)2020) ///
    xtitle("Year", margin(small)) ///
    
	
graph save unemployment_hd.gph, replace

////////////////////////////////////////////////////////////////////////// VACANCIES ///////////////////////////////////////////////////////////////////////

keep quarter bs_tfp_news_extend log_tfp_util log_tfp log_real_gdp log_real_consumption log_real_stock_price log_hours log_vacancies log_unemployment ffr inflation recession qdate t

tsfilter hp lhs_var = log_vacancies , smooth(1600)

g xx = bs_tfp_news_extend

local nlags 4
local irflag = 20

* compute necessary leads and lags
forvalues lag = 1/`nlags' {
    foreach yy of varlist xx lhs_var  {
        cap drop `yy'L`lag'
        g `yy'L`lag' = `yy'[_n-`lag']
    }
}

foreach yy of varlist xx lhs_var {
    forvalues hh = 0/`irflag' {

        cap drop `yy'F`hh'
        g `yy'F`hh' = `yy'[_n+`hh']
    }
}

forvalues hh = 1/`irflag' {
    cap drop xxL`hh'
    g xxL`hh' = xx[_n-`hh']
}


lpirf `nlags' `irflag'

matrix list e(b)

*bootstrap _b, reps(10) nodrop: lpirf `nlags' `irflag'

lpirf `nlags' `irflag'

matrix list e(b)

histdeco `nlags' `irflag'

summ lhs_var

local min = r(min)
g recession_ind = r(max) if recession == 1 
replace recession_ind = r(min) if recession == 0 


twoway ///
    ( area recession_ind quarter if lhs_var <. & xx_hd<., color(gs14) base(`min') ) ///
    ( line lhs_var quarter if lhs_var <. & xx_hd<., lcolor(blue) lwidth(medthick) lpattern(solid) ) ///
    ( line xx_hd quarter if lhs_var <. & xx_hd<., lcolor(red) lwidth(medthick) lpattern(dash) ) ///
    , ///
    title("Vacancies") ///
    legend(off) ///
    ytitle("Levels") ///
    xlabel(1960(10)2020) ///
    xtitle("Year", margin(small)) ///
    
	
graph save vacancies_hd.gph, replace


graph combine output_hd.gph consumption_hd.gph unemployment_hd.gph vacancies_hd.gph, ///
    saving(HD_Combined.gph, replace) ///
    col(2)              ///

graph export "HD_Combined.png", as(png) replace

graph combine output_hd.gph  unemployment_hd.gph , ///
    saving(HD_Combined2.gph, replace) ///
    col(1)              ///

* 2) Export the combined graph as PNG
graph export "Figure8.png", as(png) replace



