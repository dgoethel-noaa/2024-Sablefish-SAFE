#Model name
Model_1
#data file name
tem_2023_na_wh.dat
# data_reweight_switch; ==0 do not do data reweighting, ==1 do data reweighting (only used in combo with associated R file)
1
#number of growth time blocks
2
#growth block cutoffs
1960
1995
# number of weight time blocks
1
# weight block cutoffs
1960
# number of maturity time blocks
1
# maturity block cutoffs
1960
#SR type ==3, average; ==2, BH, ==1, Ricket
3
#styr for rec devs est
1930
#end yr for rec devs est (usually estimate until terminal year-1)...represents yearclass of terminal year-3 (recruitment is at age-2)
2022
#rec likelihood type
2
#bias_ramp (do bias adjustment ramp ==1)
1
#bmax    max value of bias adjustment
1.0
#b_a50  approx selectivity a50 for bias ramp (essentially offset for asc/desc limbs of ramp to account for cohorts that won't be observed in age comps hist/recent respectively)
5
#b_year_st    year bias ramp starts
1980
#b_year_end    year bias ramp reaches max
1990
#sigma_R_early_Switch
1
#sigma_R_early
0.4
#sigma_R_early_end
1975
#ph mean recruit
1
#ph rec_devs
2
#ph steepness
-6
#ph R zero
-3
#ph sigma_R
2
#ph M
2
#ph mdelta
-1
#ph Mdevs
-5
#ph Mdevs age
-5
#ph ave_F
2
#ph Fdevs
3
#ph F50
-6
#ph fish_sel
4
#ph fish2_sel
-4
#ph_fish_sel_delt     
4
#ph_fish_sel_delt alt     
-4
#ph surv_sel
4
#ph surv 2 sel
-4
#ph surv sel delt
4
#ph surv sel delt alt
-4
#ph_IFQ, sets different sel for IFQ fleets (i.e. Fish4)
4
#ph_IFQ_block2, sets different sel for IFQ fleets (i.e. Fish4) in recent years (due to recent high yc)
4
#ph_LL_block2, sets different sel for LL survey (i.e. Srvy3) in recent years (due to recent high yc)
4
#yr_sel_chg_fish, sets year for change in sel for IFQ fleets (i.e. Fish4) when ph_IFQ_block2>0
2016
#yr_sel_chg_srv3, sets year for change in sel for LL survey (i.e. Srvy3) when ph_LL_block2>0
2016
#ph surv1 q
1
#ph surv2 q
1
#ph surv3 q, uses survey 1 q
-3
#ph surv4 q, uses survey 2 q
-3
#ph surv5 q
-3
#ph surv6 q
1
#ph surv7 q
1
#ph surv8 q
1
#ph srv2 q2, est different catchability of srv2
-3
#ph_q_LL_srv_rec
-3
#ph_q_IFQ_rec
1
####### Selectivity Types ##############################################
#fish1 sel type
2
#fish2 selc type
2
#fish 3 sel type
3
#fish4 sel type
2
#fish5 sel type, post-IFQ recent time block
2
#srv1 sel type
2
#srv2 sel type
2
#srv7 sel type
3
#srv10 sel type, LL survey recent time block
2
#age of full sel fishery (age that stop est selectivity par for fishery)
0
#age of full sel survey (age that stop est selectivity par for survey)
0
#### Priors ############################################
#prior mean for M
0.1
#prior CV for M
0.1
#prior mean steepness
2
#prior CV steepnes
0.1
#prior mean sigma_R
-0.15
#prior CV sigma_R
0.2
#prior mean surv1 q
7.857
#prio CV surv1 q
0.33
#prior mean surv2 q
4.693
#prior CV surv2 q
0.242
#prior mean surv3 q
7.954994735
#prior CV surv3 q
5
#prior mean surv4 q
3.962994792
#prior CV surv4 q
5
#prior mean surv5 q
4.967
#prior CV surv5 q
0.328
#prior mean surv6 q
4.967
#prior CV surv6 q
0.328
#prior mean surv7 q
0.692
#prior CV surv67 q
0.295
#prior mean surv8 q
4.967
#prior CV surv8 q
0.328
#YR catchwt changes (I think for LL weighting)
1977
#catch LL weight fish1
50
#catch LL weight fish3
50
#wt surv1 RPW
0
#wt surv2 RPW
0
#wt DOM LL Srvy RPN
0.448
#wt JPN LL Srvy RPN
0.448
#wt surv5 DOM CPUE RPW
0.448
#wt surv6 JPN CPUE RPW
0.448
#wt GOA Trawl bio
0.448
#wt surv8 bio
0
#wt fish1 age comp
7.8
#wt surv1 age comp
7.95
#wt fish1 size comp
1
#wt surv1 size comp
1
#wt fish2 size comp
0
#wt surv2 size comp
1
#wt fish 3 size comp
4.1
#wt srv7 size comp
7.25
#wt fish4 size
0
#wt surv7 age
0
#wt surv extra size 
0
#wt surv5 size
0
#wt fish6 size
0
#wt srv6 size
0
#wt rec var penalty
1.5
0.797868466479416 #wt fish1 age comp iter
3.72382880611179 #wt surv1 age comp iter
1.27151115308662 #wt surv2 age comp iter
5.21629903889062 #wt fish1 size comp male iter
4.94489032547033 #wt fish1 size comp female iter
1.11519421735366 #wt surv1 size comp male iter
1.49981566417443 #wt surv1 size comp female iter
0.902165573415091 #wt surv2 size comp male iter
1.26841555181968 #wt surv2 size comp female iter
0.255150139457704 #wt fish3 size comp male iter
0.35007890900583 #wt fish3 size comp female iter
0.450498802823617 #wt srv7 size comp male iter
0.672532119442201 #wt srv7 size comp female iter
#wt sel fish1 regularity (i.e., smoother)
10
#wt sel fish2 regularity (i.e., smoother)
10
#wt sel fish3 regularity (i.e., smoother)
10
#wt sel fish4 regularity (i.e., smoother)
10
#wt sel fish1 dome (i.e., prevent strong dome)
10
#wt sel fish2 dome (i.e., prevent strong dome)
10
#wt sel fish3 dome (i.e., prevent strong dome)
10
#wt sel fish4 dome (i.e., prevent strong dome)
10
#wt sel surv1 regularity (i.e., smoother)
10
#wt sel surv2 regularity (i.e., smoother)
10
#wt sel surv1 dome (i.e., prevent strong dome)
10
#wt sel surv2 dome (i.e., prevent strong dome)
0
#wt F regularity
0.1
#wt M regularity
0.1
#wt_q_priors
0
#wt_M_priors
1
#wt_sigr_priors
0
#hist_hal_prop, add data for BS flag (no idea what this is)  
0.1
#yield ratio
0.757585413766995
### Expansion factor after October 1 for 2018 catch was 1.162208, dat file maker variable mean(ratios)
