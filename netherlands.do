cap log close
clear all

cd "C:\Users\anaaz\Documents\Master of Economics\Economics of Sex Work\Data"

log using netherlands.log, replace

/* Do-file for population estimation of sex workers in the Netherlands */

* Heterogeneous truncated Poisson regression

use netherlands_dummies

set matsize 11000
gen one=1
mkmat one

tempname NPresults

postfile `NPresults' counter ll_null ll_model df AIC BIC  ///
LRchi_null LRp_null LRdf_null LMtest					 ///
T_Poiss T_left T_right using "NPresults.dta", replace
 
#delimit ;
local covsets "	
				"one"

				"L1 L2 L3 L4"
				"AP1 AP2 AP3 AP4"
				"C1 C2 C3 C4 C5 C6 C7"
				"G1 G2 G3" 
				"P1 P2 P3 P4 P5"
				"H1 H2 H3 H4 H5"
				"N1 N2 N3 N4 N5 N6 N7" 
				"PI1 PI2 PI3"
				"LA1 LA2 LA3 LA4"
				"A1 A2 A3 A4 A5 A6"
				"NS1 NS2 NS3"
				"SI1 SI2 SI3"
				"TA1 TA2 TA3"
				"HY1 HY2 HY3 HY4"
				"R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15"
				"S1 S2 S3 S4 S5"
				
				" S1 S2 S3 S4 S5
				LA1 LA2 LA3 LA4
				N1 N2 N3 N4 N5 N6 N7
				HY1 HY2 HY3 HY4
				R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15
				H1 H2 H3 H4 H5 
				P1 P2 P3 P4 P5 "
				
				"AP1 AP2 AP3 AP4
				C1 C2 C3 C4 C5 C6 C7
				G1 G2 G3
				PI1 PI2 PI3
				A1 A2 A3 A4 A5 A6
				NS1 NS2 NS3
				SI1 SI2 SI3
				TA1 TA2 TA3 
				L1 L2 L3 L4 "	
				
				" A1 A2 A3 A4 A5 A6
				G1 G2 G3
				LA1 LA2 LA3 LA4
				P1 P2 P3 P4 P5
				N1 N2 N3 N4 N5 N6 N7
				L1 L2 L3 L4
				C1 C2 C3 C4 C5 C6 C7
				H1 H2 H3 H4 H5
				SI1 SI2 SI3
				HY1 HY2 HY3 HY4
				TA1 TA2 TA3
				PI1 PI2 PI3
				NS1 NS2 NS3
				AP1 AP2 AP3 AP4
				R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15
				S1 S2 S3 S4 S5 "	
				
				" ;
#delimit cr

local counter = 1

foreach x of local covsets {

	* Truncated Poisson regression
	tpoisson reviews `x'
	
	estimates store mo_`counter'
	
	gen result7 = e(chi2)
	gen result8 = e(p)
	gen result9 = e(df_m)
	
	estat ic
	
	matrix list r(S)
	matrix S = r(S)
	svmat S, names(result)
	
	predict LP_`counter', xb
	matrix V_`counter'=get(VCE)
	
	* Calculating lambda_i:
	gen lambda_`counter' = exp(LP_`counter')
	
	* LM test of overdispersion (poisson vs negative binomial):
	gen err_`counter' = reviews - lambda_`counter'
	egen ybar_`counter'=mean(reviews)
	mkmat ybar_`counter'
	matrix nybar_`counter'=one'*ybar_`counter'
	mkmat err_`counter', matrix(err_`counter')
	mkmat lambda_`counter', matrix(lambda_`counter')
	
	matrix lmtest_`counter'=((err_`counter''*err_`counter'-nybar_`counter')*(err_`counter''*err_`counter'-nybar_`counter'))*inv(2*lambda_`counter''*lambda_`counter')
	matrix list lmtest_`counter'
	scalar result10 = lmtest_`counter'[1,1]
	
	* Calculating p_i = 1-exp(-lambda) = P(Y>0):
	gen p_`counter' = 1-exp(-lambda_`counter')

	* Calculating pinv = 1/p:
	gen pinv_`counter' = 1/(p_`counter')

	mkmat pinv_`counter'

	* Computing the truncated poisson estimator:
	matrix TPoiss_`counter' = one'*pinv_`counter'

	* Computing (1-p_i)/p_i² - the first part of the variance:
	gen a0_`counter'=(1-p_`counter')/(p_`counter'*p_`counter') 
	mkmat a0_`counter'

	matrix var1_`counter'=one'*a0_`counter'

	* Computing a1 - the second part of the variance:
	gen a1_`counter'=exp(log(lambda_`counter')-lambda_`counter')/((1-exp(-lambda_`counter'))^2) 
	mkmat a1_`counter' 

	* Creating a (nx(k+1)) matrix "grad" with the data of the covariates:
	mkmat `x' one, matrix(grad_`counter') 
	
	* Creating a (1x(k+1)) matrix "fullgrad" with a1*x (where x is all covariates): 
	matrix fullgrad_`counter' =a1_`counter''*grad_`counter' 

	* Creating variance matrices:
	* Note: var21 is a (1x(k+1))matrix, var2, var1 and var are (1x1) matrices
	matrix var21_`counter' = fullgrad_`counter'*V_`counter' 
	matrix var2_`counter' = var21_`counter'*fullgrad_`counter'' 
	matrix var_`counter' = var1_`counter'+var2_`counter'
	scalar var_`counter' = trace(var_`counter')
	
	* Computing point estimate:
	scalar TPoiss_`counter'=trace(TPoiss_`counter') 
	
	* Computing confidence intervals using point estimate and variance:
	scalar TP_left_`counter' = TPoiss_`counter'-1.96*sqrt(var_`counter') 
	scalar TP_right_`counter' = TPoiss_`counter'+1.96*sqrt(var_`counter') 
	scalar list TPoiss_`counter' TP_left_`counter' TP_right_`counter'

	post `NPresults' (`counter') (result2) (result3) (result4) (result5)  ///
	 (result6) (result7) (result8) (result9) (result10)							   ///
	(TPoiss_`counter') (TP_left_`counter') (TP_right_`counter')
	
	drop result*
	
local counter = `counter' + 1
	}
	
	postclose `NPresults'

clear
********************************************************************************

* Zelterman estimation:

use netherlands_dummies

gen reviews_z = reviews
replace reviews_z=0 if reviews==1
replace reviews_z=1 if reviews==2
replace reviews_z=. if reviews>2

set matsize 11000

gen one=1
mkmat one

tempname NZresults

postfile `NZresults' counter ll_null ll_model df AIC BIC  ///
LRchi_null LRp_null LRdf_null							 ///
Zelter_ Z_left_ Z_right_ using "NZresults.dta", replace

#delimit ;
local covsets "	
				"one"

				"L1 L2 L3 L4"
				"AP1 AP2 AP3 AP4"
				"C1 C2 C3 C4 C5 C6 C7"
				"G1 G2 G3" 
				"P1 P2 P3 P4 P5"
				"H1 H2 H3 H4 H5"
				"N1 N2 N3 N4 N5 N6 N7" 
				"PI1 PI2 PI3"
				"LA1 LA2 LA3 LA4"
				"A1 A2 A3 A4 A5 A6"
				"NS1 NS2 NS3"
				"SI1 SI2 SI3"
				"TA1 TA2 TA3"
				"HY1 HY2 HY3 HY4"
				"R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15"
				"S1 S2 S3 S4 S5"
				
				" S1 S2 S3 S4 S5
				LA1 LA2 LA3 LA4
				L1 L2 L3 L4
				N1 N2 N3 N4 N5 N6 N7
				HY1 HY2 HY3 HY4
				R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15 "
				
				" 
				AP1 AP2 AP3 AP4
				C1 C2 C3 C4 C5 C6 C7
				G1 G2 G3
				H1 H2 H3 H4 H5
				PI1 PI2 PI3
				A1 A2 A3 A4 A5 A6
				NS1 NS2 NS3
				SI1 SI2 SI3
				TA1 TA2 TA3 
				P1 P2 P3 P4 P5 "	
				
				" A1 A2 A3 A4 A5 A6
				G1 G2 G3
				LA1 LA2 LA3 LA4
				P1 P2 P3 P4 P5
				N1 N2 N3 N4 N5 N6 N7
				L1 L2 L3 L4
				C1 C2 C3 C4 C5 C6 C7
				H1 H2 H3 H4 H5
				SI1 SI2 SI3
				HY1 HY2 HY3 HY4
				TA1 TA2 TA3
				PI1 PI2 PI3
				NS1 NS2 NS3
				AP1 AP2 AP3 AP4
				R1 R2 R3 R4 R5 R6 R7 R8 R9 R10 R11 R12 R13 R14 R15
				S1 S2 S3 S4 S5 "
				
				" ;
#delimit cr

local counter = 1

foreach x of local covsets {
	
	* Logit regression
	logit reviews_z `x'
	
	estimates store mo_`counter'
	
	gen result7 = e(chi2)
	gen result8 = e(p)
	gen result9 = e(df_m)
	
	estat ic
	
	matrix list r(S)
	matrix S = r(S)
	svmat S, names(result)
	
	predict LP_`counter', xb
	matrix V_`counter'=get(VCE)
	
	* Calculating w_i and 1/w_i:
	gen v_`counter'=-2*exp(LP_`counter') 
	gen w_`counter' =1-exp(v_`counter') 
	gen winv_`counter'=1/w_`counter'

	mkmat winv_`counter'

	* Computing the Zelterman estimator:
	matrix Zelter_`counter'= one'*winv_`counter' 
	matrix list Zelter_`counter' 

	* Computing (1-w_i)/w_i² - the first part of the variance:
	gen a0_`counter' =(1-w_`counter' )/(w_`counter'*w_`counter') 
	mkmat a0_`counter' 

	matrix var1_`counter'=one'*a0_`counter' 

	* Computing (1-w_i)*v_i/w_i² - the second part of the variance:
	gen a1_`counter'=(1-w_`counter')*v_`counter'/(w_`counter'*w_`counter') 
	mkmat a1_`counter'

	* Creating a (nx(k+1)) matrix "grad" with the data of the covariates:
	mkmat `x' one, matrix(grad_`counter') 
	
	* Creating a (1x(k+1)) matrix "fullgrad" with a1*x (where x is all covariates): 
	matrix fullgrad_`counter' =a1_`counter''*grad_`counter' 

	* Creating variance matrices:
	* Note: var21 is a (1x(k+1))matrix, var2, var1 and var are (1x1) matrices
	matrix var21_`counter' = fullgrad_`counter'*V_`counter' 
	matrix var2_`counter' = var21_`counter'*fullgrad_`counter'' 
	matrix var_`counter' = var1_`counter'+var2_`counter'
	scalar var_`counter' = trace(var_`counter')
	
	* Computing point estimate:
	scalar Zelter_`counter'=trace(Zelter_`counter') 
	
	* Computing confidence intervals using point estimate and variance:
	scalar Z_left_`counter' = Zelter_`counter'-1.96*sqrt(var_`counter') 
	scalar Z_right_`counter' = Zelter_`counter'+1.96*sqrt(var_`counter') 
	scalar list Zelter_`counter' Z_left_`counter' Z_right_`counter'
	
	post `NZresults' (`counter') (result2) (result3) (result4) (result5)  ///
	 (result6) (result7) (result8) (result9)							///
	(Zelter_`counter') (Z_left_`counter') (Z_right_`counter')
	
	drop result*
	
local counter = `counter' + 1
	}

	postclose `NZresults'
	
	
log close
