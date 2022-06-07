

capture log close

log using "biprobit_new.log", replace

set seed 1234
use "auditor_with_stratio.dta", clear

tab auditor, gen(auditor)
sum stratio
gen stratio_r = (stratio - `r(mean)')/`r(sd)'


global controls logassets_r lev_r investments_r cash_r roa_r logprice_r Intangible_r RD_r RDmissing numyears fracnonfees_r feemissing downgrade NoRate RateC retvol_r stratio_r sumlogret_r
biprobit (bankrptobs = going_concern $controls auditor2-auditor8 year3-year15 ) (going_concern =  $controls auditor2-auditor8 year3-year15)

predict prob_b1g1, p11 
predict prob_b1g0, p10 

sum prob_b1g1 prob_b1g0

saveold "auditor_with_stratio_biprob_s1234.dta", version(12) replace





set seed 25678
use "auditor_with_stratio.dta", clear

tab auditor, gen(auditor)
sum stratio
gen stratio_r = (stratio - `r(mean)')/`r(sd)'


global controls logassets_r lev_r investments_r cash_r roa_r logprice_r Intangible_r RD_r RDmissing numyears fracnonfees_r feemissing downgrade NoRate RateC retvol_r stratio_r sumlogret_r
biprobit (bankrptobs = going_concern $controls auditor2-auditor8 year3-year15 ) (going_concern =  $controls auditor2-auditor8 year3-year15)

predict prob_b1g1, p11 
predict prob_b1g0, p10 

sum prob_b1g1 prob_b1g0

saveold "auditor_with_stratio_biprob_s25678.dta", version(12) replace

log close
