/*

This code replicates the analysis of the paper.
The replication data is "data_for_replication.dta". 
Simply load this data in Stata, then:
1. Load the functions rum_mlogt_weights and run_models_weights
2. Start the analys from line 349


*/

capture program drop run_mlogit_weights
program define run_mlogit_weights, rclass
    syntax [if], KEYvars(string) Categorical(string) COVariates(string)

    marksample touse

    // Count keyvars and set up results matrix
    local k : word count `keyvars'
    local totalcols = `k' * 3 + 1  // 3 stats per keyvar + 1 for N

    tempname mRES_cat
    matrix `mRES_cat' = J(300, `totalcols', .) // room for many outcomes
    local i_cat = 1
    local rownames_cat

    // Loop over categorical outcomes
    foreach var of local categorical {
        di as text "---- Processing variable: `var' ----"

        // Define missingness indicator
        capture drop obs_`var'
        gen obs_`var' = !missing(`var')

        // Estimate attrition model (probit on being observed)
        // Use only baseline covariates (no outcomes)
        quietly probit obs_`var' `keyvars' `covariates' if `touse'

        // Predict probability of being observed
        capture drop phat_`var'
        predict phat_`var' if e(sample)

        // Stabilised and trimmed IPW
        replace phat_`var' = 0.0001 if phat_`var' < 0.0001
        replace phat_`var' = 0.9999 if phat_`var' > 0.9999
        gen ipw_`var' = 1 / phat_`var'
        quietly sum ipw_`var' if `touse'
        scalar mean_ipw_`var' = r(mean)
        gen ipwstab_`var' = ipw_`var' / mean_ipw_`var'

        // Trim weights (1st–99th percentile)
        quietly summ ipwstab_`var', detail
        scalar p1_`var' = r(p1)
        scalar p99_`var' = r(p99)
        gen ipwtrim_`var' = ipwstab_`var'
        replace ipwtrim_`var' = p99_`var' if ipwtrim_`var' > p99_`var'
        replace ipwtrim_`var' = p1_`var'  if ipwtrim_`var' < p1_`var'

        // Weighted mlogit for observed cases
        quietly mlogit `var' `keyvars' `covariates' [pweight=ipwtrim_`var'] ///
            if obs_`var' & `touse', cluster(bcsid)

        // Get outcome levels (excluding base)
        levelsof `var' if e(sample), local(levels)
        local base = e(baseoutcome)

        // Extract marginal effects for each non-base level
        foreach lvl of local levels {
            if "`lvl'" != "`base'" {
                quietly margins, dydx(`keyvars') predict(outcome(`lvl'))
                matrix aux = r(table)
                local marginvars : colnames aux

                local rowlabel = "`var'_level`lvl'"
                local rownames_cat `"`rownames_cat' `rowlabel'"'

                local col = 1
                foreach kv of local keyvars {
                    local found = 0
                    forvalues j = 1/`=colsof(aux)' {
                        local mv : word `j' of `marginvars'
                        if "`mv'" == "`kv'" {
                            matrix `mRES_cat'[`i_cat', `col']     = aux[1, `j']
                            matrix `mRES_cat'[`i_cat', `col'+1]   = aux[2, `j']
                            matrix `mRES_cat'[`i_cat', `col'+2]   = aux[4, `j']
                            local found = 1
                            continue, break
                        }
                    }
                    if `found' == 0 {
                        matrix `mRES_cat'[`i_cat', `col']     = .
                        matrix `mRES_cat'[`i_cat', `col'+1]   = .
                        matrix `mRES_cat'[`i_cat', `col'+2]   = .
                        di as error "Warning: `kv' not found in margins output for `var' level `lvl'"
                    }
                    local col = `col' + 3
                }

                matrix `mRES_cat'[`i_cat', `totalcols'] = r(N)
                local i_cat = `i_cat' + 1
            }
        }

        // Clean up temporary variables
        drop obs_`var' phat_`var' ipw_`var' ipwstab_`var' ipwtrim_`var'
    }

    // Finalise matrix
    local rowcount_cat = `i_cat' - 1
    matrix `mRES_cat' = `mRES_cat'[1..`rowcount_cat', 1..`totalcols']
    matrix rownames `mRES_cat' = `rownames_cat'

    // Set column names
    local colnames
    foreach kv of local keyvars {
        local colnames `colnames' dydx_`kv' se_`kv' pval_`kv'
    }
    local colnames `colnames' N
    matrix colnames `mRES_cat' = `colnames'

    // Return results
    return matrix categorical = `mRES_cat'
end









capture program drop run_models_weights
program define run_models_weights, rclass
    syntax [if], KEYvars(string) Binary(string) Continuous(string) COVariates(string)

    marksample touse

    // Count keyvars
    local k : word count `keyvars'
    local totalcols = `k' * 3 + 1  // 3 stats per keyvar + 1 for N

    // Create result matrix
    tempname mRES
    matrix `mRES' = J(100, `totalcols', .)
    local i = 1
    local rownames

    // =====================================================
    //                   BINARY OUTCOMES
    // =====================================================
    foreach var of local binary {
        // Step 1. Model missingness of dependent variable
		
		
       // Define missingness indicator
        capture drop obs_`var'
        gen obs_`var' = !missing(`var')

        // Estimate attrition model (probit on being observed)
        // Use only baseline covariates (no outcomes)
        quietly probit obs_`var' `keyvars' `covariates' if `touse'

        // Predict probability of being observed
        capture drop phat_`var'
        predict phat_`var' if e(sample)

        // Stabilised and trimmed IPW
        replace phat_`var' = 0.0001 if phat_`var' < 0.0001
        replace phat_`var' = 0.9999 if phat_`var' > 0.9999
        gen ipw_`var' = 1 / phat_`var'
        quietly sum ipw_`var' if `touse'
        scalar mean_ipw_`var' = r(mean)
        gen ipwstab_`var' = ipw_`var' / mean_ipw_`var'

        // Trim weights (1st–99th percentile)
        quietly sum ipwstab_`var', detail
        scalar p1_`var' = r(p1)
        scalar p99_`var' = r(p99)
        gen ipwtrim_`var' = ipwstab_`var'
        replace ipwtrim_`var' = p99_`var' if ipwtrim_`var' > p99_`var'
        replace ipwtrim_`var' = p1_`var'  if ipwtrim_`var' < p1_`var'
		

        // Step 3. Run main model
        quietly logit `var' `keyvars' `covariates' [pweight=ipwtrim_`var']  if obs_`var' & `touse', cluster(bcsid)
        margins, dydx(`keyvars')

        // Step 4. Store results
        matrix aux = r(table)
        local marginvars : colnames aux
        local col = 1

        foreach kv of local keyvars {
            local found = 0
            forvalues j = 1/`=colsof(aux)' {
                local mv : word `j' of `marginvars'
                if "`mv'" == "`kv'" {
                    matrix `mRES'[`i', `col']     = aux[1, `j']   // dydx
                    matrix `mRES'[`i', `col'+1]   = aux[2, `j']   // se
                    matrix `mRES'[`i', `col'+2]   = aux[4, `j']   // p-value
                    local found = 1
                    continue, break
                }
            }
            if !`found' {
                matrix `mRES'[`i', `col']     = .
                matrix `mRES'[`i', `col'+1]   = .
                matrix `mRES'[`i', `col'+2]   = .
                di as error "Warning: `kv' not found in margins output for `var'"
            }
            local col = `col' + 3
        }

        matrix `mRES'[`i', `totalcols'] = r(N)
        local rownames `"`rownames' `var'"'
        local i = `i' + 1

        drop observed p_obs ipw
    }

    // =====================================================
    //                 CONTINUOUS OUTCOMES
    // =====================================================
    foreach var of local continuous {
        // Step 1. Model missingness of dependent variable
        // Define missingness indicator
        capture drop obs_`var'
        gen obs_`var' = !missing(`var')

        // Estimate attrition model (probit on being observed)
        // Use only baseline covariates (no outcomes)
        quietly probit obs_`var' `keyvars' `covariates' if `touse'

        // Predict probability of being observed
        capture drop phat_`var'
        predict phat_`var' if e(sample)

        // Stabilised and trimmed IPW
        replace phat_`var' = 0.0001 if phat_`var' < 0.0001
        replace phat_`var' = 0.9999 if phat_`var' > 0.9999
        gen ipw_`var' = 1 / phat_`var'
        quietly sum ipw_`var' if `touse'
        scalar mean_ipw_`var' = r(mean)
        gen ipwstab_`var' = ipw_`var' / mean_ipw_`var'

        // Trim weights (1st–99th percentile)
        quietly sum ipwstab_`var', detail
        scalar p1_`var' = r(p1)
        scalar p99_`var' = r(p99)
        gen ipwtrim_`var' = ipwstab_`var'
        replace ipwtrim_`var' = p99_`var' if ipwtrim_`var' > p99_`var'
        replace ipwtrim_`var' = p1_`var'  if ipwtrim_`var' < p1_`var'

        // Step 3. Run main model
        quietly regress `var' `keyvars' `covariates' [pweight=ipwtrim_`var']  if obs_`var' & `touse', cluster(bcsid)
        margins, dydx(`keyvars')

        // Step 4. Store results
        matrix aux = r(table)
        local marginvars : colnames aux
        local col = 1

        foreach kv of local keyvars {
            local found = 0
            forvalues j = 1/`=colsof(aux)' {
                local mv : word `j' of `marginvars'
                if "`mv'" == "`kv'" {
                    matrix `mRES'[`i', `col']     = aux[1, `j']
                    matrix `mRES'[`i', `col'+1]   = aux[2, `j']
                    matrix `mRES'[`i', `col'+2]   = aux[4, `j']
                    local found = 1
                    continue, break
                }
            }
            if !`found' {
                matrix `mRES'[`i', `col']     = .
                matrix `mRES'[`i', `col'+1]   = .
                matrix `mRES'[`i', `col'+2]   = .
                di as error "Warning: `kv' not found in margins output for `var'"
            }
            local col = `col' + 3
        }

        matrix `mRES'[`i', `totalcols'] = r(N)
        local rownames `"`rownames' `var'"'
        local i = `i' + 1

        drop observed p_obs ipw
    }

    // =====================================================
    //                   FINALIZE RESULTS
    // =====================================================
    local rowcount = `i' - 1
    matrix `mRES' = `mRES'[1..`rowcount', 1..`totalcols']
    matrix rownames `mRES' = `rownames'

    local colnames
    foreach kv of local keyvars {
        local colnames `colnames' dydx_`kv' se_`kv' pval_`kv'
    }
    local colnames `colnames' N
    matrix colnames `mRES' = `colnames'

    matrix list `mRES'
    return matrix results = `mRES'
end



 
**# HOW TO REPLICATE THE DATA. 
/*

FOR REPLICATION LEAVE THIS UNCOMMENTED AND SIMPLY LOAD `data_for_replication.dta' before running the code below. 

qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_3.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_1.do"
qui do "/Users/user/Dropbox/Econometrics/perserverance/datapreparation_2.do"



merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-7473-stata/stata/stata13/bcs70_2012_flatfile.dta", keepusing(B9SCQ6 B9SCQ3A B9SCQ3B B9SCQ3C B9SCQ3D B9SCQ3E B9SCQ3F B9SCQ3G B9SCQ3H B9SCQ3I B9SCQ3J B9SCQ4 B9SCQ5A B9SCQ5B B9SCQ5C B9SCQ6 B9SCQ7 B9SCQ8A B9SCQ8B B9SCQ8C B9SCQ8D B9SCQ8E B9SCQ8F B9SCQ8G B9SCQ8H B9SCQ8I B9SCQ8J B9SCQ8K B9SCQ8L B9SCQ8M B9SCQ8N B9SCQ8O B9SCQ8P B9SCQ9 B9SCQ11A)

drop _merge

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2000.dta", keepusing(vote97 votewho97 interestPolitics97)
drop _merge

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2004.dta", keepusing(vote* votewho* interestPolitics*)
drop _merge

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2012.dta", keepusing(vote* votewho* interestPolitics*)
drop _merge

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2016.dta", keepusing(vote* votewho* )
drop _merge

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2021.dta", keepusing(vote* votewho* interestPolitics* b11*)
drop _merge

/*merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-8547-stata/stata/stata13/bcs_age46_main.dta", keepusing(B10VOTE01 B10VOTEWO1 B10VOTEFLAG B10VOTE02 B10VOTEWO2)
*/
*/



keep if isInPerinatalWave==1 & pc1 !=. & (perserverance_1_sd !=. | perserverance_2_sd!=.)


gen novoted2010 = B9SCQ6==16
replace novoted2010 =. if B9SCQ6<0

gen votedConservative2010 = B9SCQ6==1
gen votedLabour2010 = B9SCQ6==2
gen votedLibDem2010= B9SCQ6==3

replace votedLabour2010=. if B9SCQ6<0 | B9SCQ6>=16
replace votedConservative2010=. if B9SCQ6<0 | B9SCQ6>=16
replace votedLibDem2010=. if B9SCQ6<0 | B9SCQ6>=16


gen votedLeft2010 = (B9SCQ6==2 | B9SCQ6==4| B9SCQ6==5 | B9SCQ6==6 | B9SCQ6==11)
gen votedRight2010 = (B9SCQ6==1 | B9SCQ6==7| B9SCQ6==8 | B9SCQ6==10 )
replace votedLeft2010 = . if B9SCQ6<0 | B9SCQ6>=16
replace votedRight2010 = . if B9SCQ6<0 | B9SCQ6>=16


 
gen novoted2005 = B9SCQ7==16
replace novoted2005 =. if B9SCQ7<0

gen votedConservative2005 = B9SCQ7==1
gen votedLabour2005 = B9SCQ7==2
gen votedLibDem2005= B9SCQ7==3

replace votedLabour2005=. if B9SCQ7<0 | B9SCQ7>=16
replace votedConservative2005=. if B9SCQ7<0 | B9SCQ7>=16
replace votedLibDem2005=. if B9SCQ7<0 | B9SCQ7>=16

gen votedLeft2005 =  (B9SCQ7==2 | B9SCQ7==4| B9SCQ7==5 | B9SCQ7==6 | B9SCQ7==11)
gen votedRight2005 = (B9SCQ7==1 | B9SCQ7==7| B9SCQ7==8 | B9SCQ7==10 )
replace votedLeft2005 = . if  B9SCQ7<0 | B9SCQ7>=16
replace votedRight2005 = . if B9SCQ7<0 | B9SCQ7>=16




gen novoted1997 = vote97 == 2
replace novoted1997 =. if vote97 >2

gen votedConservative1997 = votewho97== 1
gen votedLabour1997 = votewho97== 2
gen votedLibDem1997= votewho97== 3

replace votedConservative1997 = . if votewho97<0 |votewho97 >=98
replace votedLabour1997 = . if votewho97< 0 |votewho97 >=98
replace votedLibDem1997= . if votewho97<0|votewho97 >=98

gen votedLeft1997 = (votewho97==2 | votewho97==4| votewho97==5 | votewho97==6)
gen votedRight1997 = (votewho97==1 )
replace votedLeft1997 = . if votewho97<0 | votewho97>=98
replace votedRight1997 = . if votewho97<0 | votewho97>=98

 


gen novoted2001 = vote01 ==2
replace novoted2001 =. if vote01 <0

gen votedConservative2001 = votewho01== 1
gen votedLabour2001= votewho01== 2
gen votedLibDem2001= votewho01== 3

replace votedConservative2001 = . if votewho01== .|votewho01<0
replace votedLabour2001= . if votewho01== .|votewho01<0
replace votedLibDem2001= . if votewho01== .|votewho01<0

gen votedLeft2001 = (votewho01==2 | votewho01==4| votewho01==5 | votewho01==6)
gen votedRight2001 = (votewho01==1 |votewho01== 7 )
replace votedLeft2001 = . if votewho01<0 | votewho01>=8
replace votedRight2001 = . if votewho01<0 | votewho01>=8

 




gen novoted2015 = vote15 ==2
replace novoted2015 =. if vote15 <0

gen votedConservative2015 = votewho15== 1
gen votedLabour2015= votewho15== 2
gen votedLibDem2015= votewho15== 3

replace votedConservative2015 =. if votewho15==. |votewho15<0 |votewho15>=9
replace votedLabour2015=. if votewho15==. |votewho15<0 |votewho15>=9
replace votedLibDem2015=. if votewho15==. |votewho15<0 |votewho15>=9

gen votedLeft2015 = (votewho15==2 | votewho15==4| votewho15==5 | votewho15==6)
gen votedRight2015 = (votewho15==1 |votewho15== 7 )
replace votedLeft2015 = . if votewho15<0 | votewho15>=9
replace votedRight2015 = . if votewho15<0 | votewho15>=9



gen novoted2017 = vote17 ==2
replace novoted2017 =. if vote17 <0

gen votedConservative2017 = votewho17== 1
gen votedLabour2017= votewho17== 2
gen votedLibDem2017= votewho17== 3

replace votedConservative2017 = . if votewho17<0
replace votedLabour2017= . if votewho17<0
replace votedLibDem2017= . if votewho17<0

gen votedLeft2017 = (votewho17==2 | votewho17==4| votewho17==5 | votewho17==6)
gen votedRight2017 = (votewho17==1 |votewho17== 7 )
replace votedLeft2017 = . if votewho17<0 | votewho17>=8
replace votedRight2017 = . if votewho17<0 | votewho17>=8




gen novoted2019 = vote19 ==2
replace novoted2019 =. if vote19 <0

gen votedConservative2019 = votewho19== 1
gen votedLabour2019= votewho19== 2
gen votedLibDem2019= votewho19== 3

replace votedConservative2019 = . if votewho19<0 
replace votedLabour2019= . if votewho19<0 
replace votedLibDem2019= . if votewho19<0 

gen votedLeft2019 = (votewho19==2 | votewho19==4| votewho19==5 | votewho19==6)
gen votedRight2019 = (votewho19==1 |votewho19== 7 |votewho19== 9 )
replace votedLeft2019 = . if votewho19<0 | votewho19>=10
replace votedRight2019 = . if votewho19<0 | votewho19>=10


replace interestPolitics21=. if interestPolitics21<0
replace interestPolitics12=. if interestPolitics12<0
replace interestPolitics01=. if interestPolitics01<0
replace interestPolitics97=. if interestPolitics97>=8

rename (interestPolitics97 interestPolitics01 interestPolitics12 interestPolitics21)(interestPolitics1997 interestPolitics2001 interestPolitics2012 interestPolitics2021)


replace voteRef = . if voteRef <0
replace voteRefWhat = . if voteRefWhat <0
recode voteRef (2=0)
gen votedForBrexit=voteRefWhat == 1 
replace votedForBrexit = . if voteRefWh==.

	* Reverse-code where higher value = lower trait level
gen j129_r = 2 - j129   // cannot concentrate
gen j134_r = 2 - j134   // temper
gen j135_r = 2 - j135   // teases
gen j138_r = 2 - j138   // bored
gen j142_r = 2 - j142   // interferes
gen j147_r = 2 - j147   // mood swings
gen j149_r = 2 - j149   // anxious
gen j156_r = 2 - j156   // unhappy/tearful
gen j160_r = 2 - j160   // quarrels
gen j163_r = 2 - j163   // destroys belongings
gen j169_r = 2 - j169   // bullies
gen j170_r = 2 - j170   // sullen
gen j177_r = 2 - j177   // fails to finish
gen j158_r = 2 - j158   // forgetful
gen j175_r = 2 - j175   // frustrated


* Standardize (z-score) variables
foreach var in j129_r j138_r j139 j155 j158_r j174 j177_r ///
               j135_r j142_r j159 j160_r j163_r j134_r j170_r ///
               j178a j148 ///
               j128 j131 j137 j143 j145 j147 j149 j156 j175_r j178b ///
               j127 j157  {
    egen z_`var' = std(`var')
}


**# Personality measures
 

* Conscientiousness proxy
egen conscientiousness = rowmean(z_j129_r z_j138_r z_j139 z_j155 z_j158_r z_j174 z_j177_r)

* Agreeableness proxy
egen agreeableness = rowmean(z_j135_r z_j142_r z_j159 z_j160_r z_j163_r z_j134_r z_j170_r)

* Extraversion proxy
* Standardize the selected Extraversion-related variables
egen z_j061 = std(j061)   // Talkative with friends
egen z_j062 = std(j062)   // Talkative with teacher
egen z_j063 = std(j063)   // Talkative with friends
egen z_j064 = std(j064)   // Talkative with teacher
egen z_j122 = std(j122)   // Popularity with peers
egen z_j123 = std(j123)   // Friendships
egen z_j124 = std(j124)   // Boldness

* Create the Extraversion composite as the average of standardized items
egen extraversion = rowmean(z_j063 z_j064 z_j122 z_j123 z_j124 z_j061 z_j062 )

 

* Neuroticism proxy
egen neuroticism = rowmean(z_j128 z_j131 z_j137 z_j143 z_j145 z_j147 z_j149 z_j156 z_j175_r z_j178b)

* Openness proxy (very weak)
* Binary self-reports for our new measure of openness.

foreach var in m85 m86 m89 m90 m91 m92 m93 m94 m95 m96 m97 {
	replace `var'=. if `var' <0 | k078>3
	egen z_`var'= std(`var') 
}

foreach var in m235 m233 m243 m245 m247 m249 {
	replace `var'=. if `var'<0
	egen z_`var'= std(`var')
}

* Create Openness to Experience proxy as the average of standardized items
egen openness = rowmean(z_m85 z_m86 z_m89 z_m90 z_m91 z_m92 z_m93 z_m94 z_m95 z_m96 z_m97 z_m235 z_m233 z_m243 z_m245 z_m247 z_m249)




* Check internal consistency (Cronbach's alpha)
matrix results = J(5, 3, .)
matrix rownames results = Extroversion Openness Conscientiousness Agreeableness Neuroticism
matrix colnames results = Items AvgCov Alpha

alpha z_j063 z_j064 z_j122 z_j123 z_j124 z_j061 z_j062 	// Extroversion
matrix results[1,1] = r(k), r(cov), r(alpha)

alpha z_m85 z_m86 z_m89 z_m90 z_m91 z_m92 z_m93 z_m94 z_m95 z_m96 z_m97 z_m235 z_m233 z_m243 z_m245 z_m247 z_m249  // openness
matrix results[2,1] = r(k), r(cov), r(alpha)


alpha z_j129_r z_j138_r z_j139 z_j155 z_j158_r z_j174 z_j177_r     // Conscientiousness
matrix results[3,1] = r(k), r(cov), r(alpha)

alpha z_j135_r z_j142_r z_j159 z_j160_r z_j163_r z_j134_r z_j170_r // Agreeableness
matrix results[4,1] = r(k), r(cov), r(alpha)

alpha z_j128 z_j131 z_j137 z_j143 z_j145 z_j147 z_j149 z_j156 z_j175_r z_j178b // Neuroticism
matrix results[5,1] = r(k), r(cov), r(alpha)
* Display matrix
matrix list results

outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/alpha_summary",mat(results) replace  f(%9.3f)



* Generate indicator variables for each strong response (1 or 5)
foreach var in B9SCQ3A B9SCQ3B B9SCQ3C B9SCQ3D B9SCQ3E B9SCQ3F B9SCQ3G B9SCQ3H B9SCQ3I B9SCQ3J {
    gen strong_`var' = inlist(`var', 1, 5)
	replace strong_`var'=. if `var' <0 |`var' ==.
}

* Create a composite "strong views" count (0–10)
egen strong_views = rowtotal(strong_B9SCQ3A strong_B9SCQ3B strong_B9SCQ3C strong_B9SCQ3D strong_B9SCQ3E ///
                              strong_B9SCQ3F strong_B9SCQ3G strong_B9SCQ3H strong_B9SCQ3I strong_B9SCQ3J), missing

* Optional: drop intermediate indicators
drop strong_B9SCQ3*


/*
Political view variables
B9SCQ3A B9SCQ3B B9SCQ3C B9SCQ3D B9SCQ3E B9SCQ3F B9SCQ3G B9SCQ3H B9SCQ3I B9SCQ3J B9SCQ4 B9SCQ5A B9SCQ5B B9SCQ5C B9SCQ6 B9SCQ7 B9SCQ8A B9SCQ8B B9SCQ8C B9SCQ8D B9SCQ8E B9SCQ8F B9SCQ8G B9SCQ8H B9SCQ8I B9SCQ8J B9SCQ8K B9SCQ8L B9SCQ8M B9SCQ8N B9SCQ8O B9SCQ8P B9SCQ9 B9SCQ11A
*/

foreach var of varlist B9SCQ3A B9SCQ3B B9SCQ3C B9SCQ3D B9SCQ3E B9SCQ3F B9SCQ3G B9SCQ3H B9SCQ3I B9SCQ3J   {
	replace `var' =. if `var'<0	
}

cluster kmeans strong_views, k(4) name(cluster4) measure(L2) start(random) 

gen strong_views_clusters = .
replace strong_views_clusters = 1 if strong_views ==0 | strong_views==1
replace strong_views_clusters = 2 if strong_views ==2 | strong_views==3
replace strong_views_clusters = 3 if strong_views ==4 | strong_views==5| strong_views==6
replace strong_views_clusters = 4 if strong_views ==7 | strong_views==8| strong_views==9 | strong_views==10



foreach var of varlist B9SCQ3A B9SCQ3B B9SCQ3C B9SCQ3D B9SCQ3E B9SCQ3F B9SCQ3G B9SCQ3H B9SCQ3I B9SCQ3J {
    * Create new variable with _3cat suffix
    gen `var'_3cat = .
    * Assign new categories
    replace `var'_3cat = 1 if inlist(`var', 1, 2)
    replace `var'_3cat = 2 if inlist(`var', 3)
    replace `var'_3cat = 3 if inlist(`var', 4, 5)
    * Add labels for clarity
    label define `var'_3cat_lbl 1 "Agree (1-2)" 2 "neither (3)" 3 "disagree", replace
    label values `var'_3cat `var'_3cat_lbl
    * Verify new variable
    tab `var' `var'_3cat, m
}


* Interest in politics
replace B9SCQ4=. if  B9SCQ4<0

* Effective political activism 
/*
B9SCQ5A  // attended a public meeting
B9SCQ5B  // attended a protest
B9SCQ8A  // member of political party
B9SCQ8B  // member of trade union
B9SCQ8C  // member of environmental group
B9SCQ8L  // member of feminist group
*/
recode B9SCQ5A (2=0)
recode B9SCQ5B (2=0)
foreach var in B9SCQ5A B9SCQ5B B9SCQ8A B9SCQ8B B9SCQ8C B9SCQ8L {
    replace `var' = . if `var' < 0
}

egen activism = rowtotal(B9SCQ5A B9SCQ5B B9SCQ8A B9SCQ8B B9SCQ8C B9SCQ8L), missing
gen any_activism = activism > 0 if !missing(activism)


**# Analysis 

	
run_mlogit_weights , ///
	keyvars(pc1 conscientiousness agreeableness extraversion neuroticism openness) ///
    categorical(B9SCQ3A_3cat B9SCQ3B_3cat B9SCQ3C_3cat B9SCQ3D_3cat B9SCQ3E_3cat B9SCQ3F_3cat B9SCQ3G_3cat B9SCQ3H_3cat ///
	B9SCQ3I_3cat B9SCQ3J_3cat ) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6  )
matrix w = r(categorical)
matrix list w
outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/attitudes", mat(w) nobox format(%9.3f) replace

	
	
run_mlogit_weights , ///
	keyvars(pc1 conscientiousness agreeableness extraversion neuroticism openness) ///
    categorical(strong_views_clusters ) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6  )
matrix w = r(categorical)
matrix list w
outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/strongviews", mat(w) nobox format(%9.3f) replace

	
run_mlogit_weights , ///
	keyvars(pc1 conscientiousness agreeableness extraversion neuroticism openness) ///
    categorical(B9SCQ4) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6  )
matrix w = r(categorical)
matrix list w
outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/interestPols", mat(w) nobox format(%9.3f) replace
	
	
run_models_weights , ///
	keyvars(pc1 conscientiousness agreeableness extraversion neuroticism openness) ///
    binary(any_activism votedForBrexit) ///
    continuous(loginc) ///
    covariates($set1 $set2 $set3 $set4 $set5 $set6  )
	
matrix w = r(results)
matrix list w
outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/activism", mat(w) nobox format(%9.3f) replace

**# Panel analysis: You need to change the way standard errors are calculated to cluster(bcsid) in the estimating functions before running this. 

*Voting behaviour as a panel.

rename(votewho97 votewho01 votewho10 votewho05 votewho15 votewho17 votewho19)(votewho1997 votewho2001 votewho2010 votewho2005 votewho2015 votewho2017 votewho2019)

preserve

	keep 	bcsid pc1 conscientiousness agreeableness extraversion neuroticism openness  votewho* novoted* votedConservative* votedLabour* votedLibDem* votedLeft* votedRight*  interestPolitics* c_delivery_age_mum c_delivery_age_mum_2 c_region_birth c_mumsmoked c_abnormal_gest c_sexbirth c_birthweight c_stayedEd_mum c_stayedEd_dad c_parity c_prevDeathUnder7 c_prevDeathAfter7 c_prevStillbirth c_prevMiscarries  c_breastfed c_inpatientAt5 c_hearSeeSpeak_issueAt5 c_childReadAt_5 c_mum_nonqual c_mum_secondary c_mum_degree c_dad_nonqual c_dad_secondary c_dad_degree c_mumWorksFTAt5 c_mumWorksPTAt5   c_nonwhite c_noSiblingsAt5 c_childHhldAt10 c_inpatientBefore10  c_fatherUnemployed10 c_fatherStayhome10 c_motherUnemployed10 c_motherStayhome10 c_benefit c_incband* c_accTenureAt101 c_accTenureAt102 c_accTenureAt103 c_accTenureAt104 c_accTenureAt105 c_accTenureAt106 c_accTenureAt107 c_roomNumAt10   c_motherSmokesat10 c_fatherSmokesat10 c_eyesightIssue10 c_hearingIssue10 c_speechIssue10 c_heightAt10_sd c_weightAt10_sd c_disengagedAt10 loginc 

	

	reshape long novoted votedConservative votedLabour votedLibDem votedLeft votedRight interestPolitics votewho, i(bcsid) j(year)

	encode bcsid, generate(personid_num)
	
	bysort person : gen t = _n
	replace votewho =. if votewho>9 | votewho <0
	xtset person t
	bysort personid_num: gen diff = d.votewho
	replace diff = 1 if (diff <0 | diff >0 ) & diff !=.
	

	run_models_weights , ///
		keyvars(pc1 conscientiousness agreeableness extraversion neuroticism 		openness) ///
		binary(novoted votedConservative votedLabour votedLibDem votedLeft 		  votedRight diff) ///
		continuous(loginc  ) ///
		covariates($set1 $set2 $set3 $set4 $set5 $set6  )
			
matrix w = r(results)
matrix list w
outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/votingBehav", mat(w) nobox format(%9.3f) replace


	run_mlogit_weights , ///
		keyvars(pc1 conscientiousness agreeableness extraversion neuroticism 		openness) ///
		categorical(interestPolitics) ///
		covariates($set1 $set2 $set3 $set4 $set5 $set6  )
	matrix w = r(categorical)
	matrix list w
	outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/interestPols_panel", mat(w) nobox format(%9.3f) replace

	
	
	
	
	
	
*-------------------------------------------*
* 1. Attrition analysis (Did not vote)
*-------------------------------------------*
local var novoted
capture drop obs_`var'
gen obs_`var' = !missing(`var')

probit obs_`var' pc1 conscientiousness agreeableness extraversion neuroticism openness ///
      $set1 $set2 $set3 $set4 $set5 $set6

* Get marginal effects
margins, dydx(pc1 conscientiousness agreeableness extraversion neuroticism openness)

* Export to LaTeX
esttab using "/Users/user/Dropbox/Econometrics/politicalBehaviour/attrition_voting_skills.tex", ///
    cells("b(star fmt(3)) se(par fmt(3))") replace label title("Probit Marginal Effects: Attrition and Voting Behavior") ///
    keep(pc1 conscientiousness agreeableness extraversion neuroticism openness) ///
    mgroups("Attrition (Did not vote)" , pattern(1 0)) noobs

*-------------------------------------------*
* 2. Voting choice (Who voted)
*-------------------------------------------*

 


local var votewho
capture drop obs_`var'
gen obs_`var' = !missing(`var')

probit obs_`var' pc1 conscientiousness agreeableness extraversion neuroticism openness ///
      $set1 $set2 $set3 $set4 $set5 $set6 if novoted==0

* Get marginal effects
margins, dydx(pc1 conscientiousness agreeableness extraversion neuroticism openness)

* Export to the same LaTeX table (append)
esttab using "/Users/user/Dropbox/Econometrics/politicalBehaviour/attrition_voting_skills.tex", ///
    cells("b(star fmt(3)) se(par fmt(3))") append label ///
    keep(pc1 conscientiousness agreeableness extraversion neuroticism openness) ///
    mgroups("Voting choice (Who voted)" , pattern(1 0)) noobs

	
	
	
	
restore	



**# Descriptive statistics


* Install estout if not already installed

* Run summarize with detail and store results
estpost summarize pc1 conscientiousness agreeableness extraversion neuroticism openness novoted* votedConservative* votedLabour* votedLibDem* votedLeft* votedRight* interestPolitics* c_delivery_age_mum c_delivery_age_mum_2 c_region_birth c_mumsmoked c_abnormal_gest c_sexbirth c_birthweight c_stayedEd_mum c_stayedEd_dad c_parity c_prevDeathUnder7 c_prevDeathAfter7 c_prevStillbirth c_prevMiscarries c_breastfed c_inpatientAt5 c_hearSeeSpeak_issueAt5 c_childReadAt_5 c_mum_nonqual c_mum_secondary c_mum_degree c_dad_nonqual c_dad_secondary c_dad_degree c_mumWorksFTAt5 c_mumWorksPTAt5 c_nonwhite c_noSiblingsAt5 c_childHhldAt10 c_inpatientBefore10 c_fatherUnemployed10 c_fatherStayhome10 c_motherUnemployed10 c_motherStayhome10 c_benefit c_incband1 c_incband2 c_incband3 c_incband4 c_incband5 c_incband6 c_incband7 c_accTenureAt101 c_accTenureAt102 c_accTenureAt103 c_accTenureAt104 c_accTenureAt105 c_accTenureAt106 c_accTenureAt107 c_roomNumAt10 c_motherSmokesat10 c_fatherSmokesat10 c_eyesightIssue10 c_hearingIssue10 c_speechIssue10 c_heightAt10_sd c_weightAt10_sd c_disengagedAt10 loginc

* Export to LaTeX
esttab using "/Users/user/Dropbox/Econometrics/politicalBehaviour/descriptive_stats.tex", cells("count mean(fmt(3)) sd(fmt(3)) min(fmt(3)) max(fmt(3))") ///
  noobs label title("Descriptive Statistics") ///
  mtitles("N" "Mean" "SD" "Min" "Max") ///
  style(tex) replace
  
   
**# Personality and cognition long run correlates


* Personality at age 51

foreach v in b11q25f b11q25g b11q25h b11q25i b11q25j b11q25o b11q25p b11q25q b11q25r b11q25s b11q25t {
    gen `v'_r = 6 - `v'   // if scale is 1–5
    // use 8 - `v' if scale is 1–7
}

egen extraversion51 = rowmean(b11q25a b11q25k b11q25f_r b11q25p_r)
egen agreeableness51 = rowmean(b11q25b b11q25l b11q25g_r b11q25q_r)
egen conscientiousness51 = rowmean(b11q25c b11q25m b11q25h_r b11q25r_r)
egen neuroticism51 = rowmean(b11q25d b11q25n b11q25i_r b11q25s_r)
egen openness51 = rowmean(b11q25e b11q25j_r b11q25o_r b11q25t_r)



* Cognitive skills
* Keep valid ranges only (assuming 0–10 for recall, and plausible animal counts for fluency)
gen mem_immediate = b11cflisn if inrange(b11cflisn, 0, 10)
gen mem_delayed   = b11cflisd if inrange(b11cflisd, 0, 10)
gen fluency       = b11cfani  if inrange(b11cfani, 0, 100)  // adjust upper bound if needed

* Standardize each variable
egen z_immediate = std(mem_immediate)
egen z_delayed   = std(mem_delayed)
egen z_fluency   = std(fluency)

* Create an overall cognitive score (average of standardized components)
egen pc151 = rowmean(z_immediate z_delayed z_fluency)

 
 
* Define the six traits (including pc1)
local traits pc1 agreeableness extraversion neuroticism openness conscientiousness  

* Create a matrix to store results: 6 rows × 3 columns (rho, p-value, mean Δ-rank)
matrix results = J(6, 3, .)
local row = 1

foreach trait of local traits {

    di as text "Processing: `trait'"

    * Standardize if needed
    egen z_`trait'_young = std(`trait')
    egen z_`trait'_old   = std(`trait'51)

    * Spearman rank correlation (stability of rank)
    spearman z_`trait'_young z_`trait'_old
    matrix R = r(rho)
    local rho = R[1,1]
    local pval = r(p)

    * Create deciles for both waves
    xtile `trait'_young_decile = z_`trait'_young, n(10)
    xtile `trait'_old_decile   = z_`trait'_old, n(10)

    * Compute absolute decile rank change
    gen `trait'_rank_change = abs(`trait'_young_decile - `trait'_old_decile)

    quietly summarize `trait'_rank_change
    local mean_change = r(mean)

    * Store results in matrix
    matrix results[`row', 1] = `rho'
    matrix results[`row', 2] = `pval'
    matrix results[`row', 3] = `mean_change'

    * Clean up temporary vars
    drop z_`trait'_young z_`trait'_old `trait'_rank_change `trait'_young_decile `trait'_old_decile

    local ++row
}

* Label matrix rows and columns
matrix rownames results = pc1 agreeableness extraversion neuroticism openness conscientiousness
matrix colnames results = Spearman_rho p_value Mean_DeltaRank

* Display matrix
matrix list results

* Export results to LaTeX
outtable using "/Users/user/Dropbox/Econometrics/politicalBehaviour/personality_stability.tex", mat(results) replace ///
    format(%6.3f)

	
 
	