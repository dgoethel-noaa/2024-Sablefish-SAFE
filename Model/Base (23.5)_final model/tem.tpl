//==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+
//
//  Statistical, separable age-structured population model for sablefish
//  Alaska Fisheries Science Center, October 2008
//  D. Hanselman:dana.hanselman@noaa.gov
//  UPDATED (and commented)  by D. Goethel: daniel.goethel@noaa.gov   (10/15/20)
//####################################################################

// ### File Structure ###########################################################################################
//  OLD_Input file:   tem.dat, as in tem(plate).....provides the settings for likelihoods, phases, data weights, SRR settings
//  Control file: tem.ctl.....same as tem.dat, but no comments...read this into tpl so need to make same changes here as to tem.dat
//  Data File: tem_2020_na_wh.dat....actual data inputs updated each year with new catch, age comp, length comp, etc.
//  Program file: tem.tpl
//  Output files: tem.rep, tem.std, proj.dat, newproj.dat
// ###############################################################################################

// ### Fleet/Area Structure #############################################################################################
// Single Area (combined across GOA, BS, AI)
// Split sexes, seperate weight at ages, age-length keys and unsexed
// All fishery catch include  catch and discards
// LL catch is adjusted (increased) for whale depredation (CPUE is not ajusted)
// DOM LL survey RPN/RPW is adjusted (increased) for whale depredation
// All fisheries and surveys are across entire model domain (GOA, BS, AI)

// srv1= Domestic LL survey RPW (not fit in model)                       (logistic selectivity)
// srv2= Japanese LL survey (not fit in model)                           (logistic selectivity)
// srv3= Domestic LL survey RPN (fit in model)                           (logistic selectivity)
// srv4 = Japanese LL survey RPN (fit in model)                          (logistic selectivity)
// srv5 = Domestic LL fishery RPW (pre-IFQ--before 1995 (fit in model)   (logistic selectivity)
// srv6 = Japanese LL fishery RPW (fit in model)                         (logistic selectivity)
// srv7 = NMFS bottom trawl survey (currently GOA only; fit in model)    (POWER FUNCTION selectivity)
// srv8 = Post-IFQ DOM LL fishery q par  RPW (fit in model)              (logistic selectivity)
// srv9 = JPN LL q par for 1990-1994 (not estimated just set to qsrv2, only used if qsrv2 is not estimated)
// srv10 = DOM LL survey RPN sel par recent time block
              // account for recent changes in selectivity and age/length comps due to strong 2014, 2016, 2017 yc

// fish1 = U.S. Longline fishery (fit in model)                           (logistic selectivity)
// fish2 = Japanese LL fishery (fit in model)                             (logistic selectivity)
// fish3 = U.S. trawl fishery (fit in model)                              (GAMMA FUNCTION selectivity)
// fish4 = Post IFQ LL fishery sel par                                    (logistic selectivity)                 // Dana has it listed as Japanese trawl fishery, but I think this is some sort of hold over no longer used...fish4 sel coefficients are clearly post_IFQ LL
// fish5 = Post IFQ LL fishery sel par second time block                  (logistic selectivity)
                        // account for recent changes in selectivity and age/length comps due to strong 2014, 2016, 2017 yc

// ############################################################################################

// ### Estimated Parameters Base Model #########################################################################################

//   RECRUITMENT
// log_mean_rec = mean recruitment estimate                                                                        (1 parameter)
// log_rec_dev = Recruitment yearly devs from mean Rec                                                             ((nages-2)+(nyrs-1)  parameters; est devs for each age at initial equil (except age 1 and plus group) and for each year (except terminal year))

//   NATURAL MORTALITY
// log_M          = Base Natural Mortality                                                                         (1 parameter)  
// log_M_devs     = Yearly Deviations from M                                                                       (NOT EST BASE; nyrs parameters)
// log_M_devs_age = Age Deviations from M                                                                          (NOT EST BASE; nages parameters)

//   FISHING MORTALITY
// log_avg_F_Fish1 = F multiplier for LL fishery                                                                   (1 parameter)  
// log_avg_F_fish3 = F multiplier for Trawl fishery                                                                (1 parameter)
// log_F_devs_fish1 = F yearly devs for LL fishery                                                                 (nyrs  parameters)
// log_F_devs_fish3 = F yearly devs for Trawl fishery                                                              (nyrs  parameters)

//   CATCHABILITY COEFFICIENTS
// log_q_srv1 = q DOM LL survey RPN/RPW                                                                            (1 parameter)
// log_q_srv2 = q JPN LL survey RPN/RPW                                                                            (1 parameter)
// log_q_srv5 = q DOM LL fishery RPW (pre-IFQ---pre-1995)                                                          (1 parameter)
// log_q_srv6 = q JPN LL fishery RPW                                                                               (1 parameter)
// log_q_srv7 = q DOM trawl survey                                                                                 (1 parameter)
// log_q_srv8 = q DOM LL fishery RPW post-IFQ  (1995+)                                                             (1 parameter)
// log_q_LL_srvy_recent = q DOM LL survey recent timeblock                                                         (1 parameter, not currently estimated)
// log_q_LL_fish_recent = q DOM LL fishery recent period (2016)                                                    (1 parameter)

//   LOGISTIC SELECTIVITY PARAMETERS FISHERY
// log_a50_fish1_f = DOM LL fishery age at 50% selection females                                                   (1 parameter)
// log_a50_fish1_m = DOM LL fishery age at 50% selection males                                                     (1 parameter)
// log_delta_fish1_f = DOM LL fishery diff in age at 50% and 90% selection females                                 (NOT EST BASE; USE FISH4 FEMALE EST; 1 parameter)
// log_delta_fish1_m = DOM LL fishery diff in age at 50% and 90% selection males                                   (NOT EST BASE; USE FISH4 FEMALE EST; 1 parameter)
// log_a50_fish2 = JPN LL fishery age at 50% selection                                                             (1 parameter)
// log_delta_fish2 = JPN LL fishery diff in age at 50% and 90% selection                                           (NOT EST BASE; USE FISH4 FEMALE EST; 1 parameter)
// log_a50_fish4_f = DOM LL fishery post-IFQ age at 50% selection females                                          (1 parameter)
// log_a50_fish4_m = DOM LL fishery post-IFQ age at 50% selection males                                            (1 parameter)
// log_delta_fish4_f = DOM LL fishery post-IFQ diff in age at 50% and 90% selection females                        (1 parameter; USED FOR ALL OTHER DELTAS)
// log_delta_fish4_m = DOM LL fishery post-IFQ diff in age at 50% and 90% selection males                          (NOT EST BASE; USE FISH4 FEMALE EST; 1 parameter)
// log_a50_fish5_f = DOM LL fishery post-IFQ recent time block age at 50% selection females                        (NOT EST BASE; EST 1 parameter when active)
// log_a50_fish5_m = DOM LL fishery post-IFQ recent time block age at 50% selection males                          (NOT EST BASE; EST 1 parameter when active)
// log_delta_fish5_f = DOM LL fishery post-IFQ recent time block diff in age at 50% and 90% selection females      (NOT EST BASE; USE FISH4 FEMALE EST WHEN ACTIVE; 1 parameter)
// log_delta_fish5_m = DOM LL fishery post-IFQ recent time block diff in age at 50% and 90% selection males        (NOT EST BASE; USE FISH4 FEMALE EST WHEN ACTIVE; 1 parameter)

//   GAMMA FUNCTION SELECTIVITY PARAMETERS FISHERY
// log_a50_fish3_f = DOM Trawl fishery age at 50% selection females                                                (1 parameter)
// log_a50_fish3_m = DOM Trawl fishery age at 50% selection males                                                  (1 parameter)
// log_delta_fish3_f = DOM Trawl fishery diff in age at 50% and 90% selection females                              (1 parameter)
// log_delta_fish3_m = DOM Trawl fishery diff in age at 50% and 90% selection males                                (1 parameter)

//   LOGISTIC SELECTIVITY PARAMETERS SURVEY
// log_a50_srv1_f = DOM LL survey age at 50% selection females                                                     (1 parameter)
// log_a50_srv1_m = DOM LL survey age at 50% selection males                                                       (1 parameter)
// log_delta_srv1_f = DOM LL survey diff in age at 50% and 90% selection females                                   (NOT EST BASE; USE SRV2 MALE EST; 1 parameter)
// log_delta_srv1_m = DOM LL survey diff in age at 50% and 90% selection males                                     (NOT EST BASE; USE SRV2 MALE EST;1 parameter)
// log_a50_srv2 = JPN LL survey age at 50% selection females                                                       (1 parameter)
// log_a50_srv2 = JPN LL survey age at 50% selection males                                                         (1 parameter)
// log_delta_srv2 = JPN LL survey diff in age at 50% and 90% selection females                                     (NOT EST BASE; USE SRV2 MALE EST;1 parameter)
// log_delta_srv2 = JPN LL survey diff in age at 50% and 90% selection males                                       (1 parameter; USED FOR ALL LL SURVEYS)
// log_a50_srv10_f = DOM LL survey recent time block age at 50% selection females                                  (NOT EST BASE; EST 1 parameter when active)
// log_a50_srv10_m = DOM LL survey recent time block age at 50% selection males                                    (NOT EST BASE; EST 1 parameter when active)
// log_delta_srv10_f = DOM LL survey recent time block diff in age at 50% and 90% selection females                (NOT EST BASE; USE SRV2 MALE EST WHEN ACTIVE; 1 parameter)
// log_delta_srv10_m = DOM LL survey recent time block diff in age at 50% and 90% selection males                  (NOT EST BASE; USE SRV2 MALE EST WHEN ACTIVE; 1 parameter)

//   POWER FUNCTION SELECTIVITY PARAMETERS SURVEY
// log_a50_srv7_f = DOM Trawl survey age at 50% selection females                                                  (1 parameter)
// log_a50_srv7_m = DOM Trawl survey age at 50% selection males                                                    (1 parameter)

//   SPR F VALUES
// mF50 = FSPR50 for  SPR calcs                                                                                    (1 parameter)
// mF40 = FSPR40 for  SPR calcs                                                                                    (1 parameter)
// mF35 = FSPR35 for  SPR calcs                                                                                    (1 parameter)

/////////////////
// estimable but phase<0 for Base
//////////////////////

// steepness               (no SRR used)
// log_Rzero               (no SRR used)
// sigr                    (recruit sigma input)
// mdelta                  (M set equal for males and females; this male offset not estimated)
// selectivity parameters  (4 selectivity variations coded, only use those parameters noted as estimated above)

// ############################################################################################


// ############################################################################################
//### Changes in August 2016 ##################################################################################### 
// added SARA outputs, new SDNR summaries at end
// standard deviation for SSB
// Fixed hessian problem on tot_bio_proj(endyr+1)
// revised sable-r-report.cxx, and mhp-func.cpp to reflect new C++ protocols
// Added executive material sections and new report sable.rep with SDs
// Cleaned up a couple of messes like srv22, srv2q
//########################################################################################

//### Changes in October 2020 #####################################################################################
// add switches and  estimation for  selectivity and survey parameters in a new recent time block to account for recent changes in targeting/age-length comp in fishery and LL survey due to large 2014, 2016, and 2017 ycs
// add ability to estimate yearly Mdevs
// #####################################################################################

//==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+

GLOBALS_SECTION
  #include "admodel.h"
  #include "statsLib.h"
  #include "qfclib.h"
  #include "mhp-s-funcs.cpp"                // Include S-compatible output functions (needs preceding)
  #include <contrib.h>
  #define EOUT(var) cout <<#var<<" "<<var<<endl;
  adstring model_name;
  adstring data_file;

DATA_SECTION
  !!CLASS ofstream evalout("evalout.prj");
// Read data from the control file
  !! ad_comm::change_datafile_name("tem.ctl");    //Read in phases, penalties and priors from "tem.ctl"
  !! *(ad_comm::global_datafile) >>  model_name;  //define model name from .ctl file
  !! *(ad_comm::global_datafile) >>  data_file;   //define .dat file from .ctl file (tem_2020_na_wh.dat...replace 2020 with current year)

  init_int    data_reweight_switch   //switch for running data reweighting, ==0 no reweighting, ==1 do data reweighting
  init_int    growth_blocks         // number of growth time blocks to accomodate (i.e., # size-age transition matrices to read in)
  init_vector growth_cutoffs(1,growth_blocks)        // break points for growth (dim is #growth_blocks-1)
  init_int    weight_blocks         // number of weight time blocks to accomodate (i.e., WAA vectors to read in)
  init_vector weight_cutoffs(1,weight_blocks)        // break points for weight
  init_int    maturity_blocks         // number of weight time blocks to accomodate (i.e., WAA vectors to read in)
  init_vector maturity_cutoffs(1,maturity_blocks)        // break points for weight
  init_int    SrType                // 3=average, 2=Bholt, 1=Ricker
  init_int    styr_rec_est
  init_int    endyr_rec_est         //typical estimate up to terminal year-1
  int nrecs_est;
  !! nrecs_est = endyr_rec_est-styr_rec_est+1;
  init_int    rec_like_type         //Type of likelihood specified
  init_int    bias_ramp             // perform SRR bias adjustment using a bias ramp as per Methot and Taylor (2011)
  init_number bmax                  // max value of bias adjustment
  init_number b_a50                 // approx selectivity a50 for setting bias ramp asc/desc limbs to account for cohorts that won't be observed in hist/rec age comp data
  init_number b_year_st             // year bias ramp starts
  init_number b_year_end            // year bias ramp reaches max
  init_int sigma_R_early_switch  // allow for reduced sigma_R in early years
  init_number sigma_R_early         // sigma_R in early years
  init_number sigma_R_early_end     // sigma_R_early end year

// Phases that general parameter estimation begins
  init_int    ph_mean_rec           // Phase for mean recruitment
  init_int    ph_recdev             // Phase for estimating recruitment deviations
  init_int    ph_steepness          // Phase for steepness (not used base), for estimating S-R parameters
  init_number ph_Rzero;             // Phase for estimating Rzero (do not estimate for base-case, neg phase)
  init_int    ph_sigr               // Phase for recruiment sigma R (not est)
  
  init_int    ph_m                  //Phase for estimating natural mortality
  init_number ph_mdelta;            // Phase for estimating change in M for males (base sets to negative; i.e., male M same as female)
  init_int    ph_Mdevs              // Phase for natural mortality deviations
  init_int    ph_Mdevs_age          // Phase for natural mortality age deviations 

  init_int    ph_avg_F              // Phase for estimating average fishing mortality
  init_int    ph_Fdev               // Phase for fishing mortality deviations
  init_int    ph_F50                // Phase for SPR F estimation
  init_int    ph_fish_sel           // Phase for estimating fishing selectivity
  init_int    ph_fish2_sel          // Phase for estimating fishing selectivity in fishery 2
  init_int    ph_fish_sel_delt      // Phase for estimating the fishery selectivity delta parameters (diff in a50 and a95)
  init_int    ph_fish_sel_delt_alt  // Phase for estimating the fishery selectivity delta parameters (diff in a50 and a95) if want to est all deltas (not just 1 par)
  init_int    ph_srv1_sel           // Phase for estimating survey 1 (and others) selectivity
  init_int    ph_srv2_sel           // Phase for estimating survey 2 selectivity
  init_int    ph_srv_sel_delt       // Phase for estimating the survey selectivity delta parameters (diff in a50 and a95)
  init_int    ph_srv_sel_delt_alt   // Phase for estimating the survey selectivity delta parameters (diff in a50 and a95) if want to est all deltas (not just 1 par)
  init_number ph_ifq;               // Phase for estimating selectivity change at onset of IFQ in LL fishery (fishery 1, use fihsery 4 selectivity parameter inputs for post-IFQ fishery)
  init_number ph_ifq_block2;        // Phase for estimating selectivity change for recent years of IFQ LL fishery (fishery 1, use fishery 5 selectivity parameter inputs)
  init_number ph_LL_block2;         // Phase for estimating selectivity change for recent years of  LL survey (survey 1, use survey 8 selectivity parameter inputs)
  init_number yr_sel_chg_fish       // Year that selectivity changes for the LL fishery
  init_number yr_sel_chg_srv1       // Year that selectivity changes for the LL survey (typical make same as for fishery)

  init_int    ph_q_srv1             // Phase for estimating catchability
  init_int    ph_q_srv2             // Phase for estimating catchability
  init_int    ph_q_srv3             // Phase for estimating catchability
  init_int    ph_q_srv4             // Phase for estimating catchability
  init_int    ph_q_srv5             // Phase for estimating catchability
  init_int    ph_q_srv6             // Phase for estimating catchability
  init_int    ph_q_srv7             // Phase for estimating catchability
  init_int    ph_q_srv8             // Phase for estimating catchability
  init_number ph_srv2_q2;           // Phase for estimating catchability of Japanese LL survey (survey 2) for 1990-1994; PAR Section sets phase to opp of this input, so only estimated if input is negative
  init_int    ph_q_LL_srv_rec       // phase for estimating catchability of LL survey in recent time block
  init_int    ph_q_IFQ_rec          // phase for estimating catchabiliyt of LL fishery in recent time block

// selectivity type for each fishery and survey
// ==1, coefficients for each age (requires LL component to ensure smooth transitions among ages)
// ==2, 2 parameter logistic (used for fixed gear and LL survey)
// ==3, 6 parameter double logistic
// ==4, exponential logistic 
  init_int    fsh1_sel_opt
  init_int    fsh2_sel_opt
  init_int    fsh3_sel_opt
  init_int    fsh4_sel_opt
  init_int    fsh5_sel_opt  //post-IFQ recent time block
  
  init_int    srv1_sel_opt
  init_int    srv2_sel_opt
  init_int    srv7_sel_opt
  init_int    srv10_sel_opt  // recent time block

// Age of full selection for fishery and survey
  init_int    n_fish_sel_ages       // Age that fishery selectivity stops estimating
  init_int    n_srv1_sel_ages       // Age that survey selectivity stops estimating

// Priors
  init_number mprior                // Prior mean for natural mortality
  number log_mprior
  !! log_mprior = log(mprior);
  init_number cvmprior              // Prior CV for natural mortality

  init_number steep_prior           // Prior mean for steepness
  init_number cv_steep_prior        // Prior CV for steepness
 
  init_number sigrprior             // Prior mean for recruitment deviations
  init_number cvsigrprior           // Prior CV for recruitment deviations

  init_number q_srv1prior           // Prior mean for catchability coefficient
  init_number cvq_srv1prior         // Prior CV for catchability coefficient
  number log_q1prior
  !! log_q1prior = log(q_srv1prior);
  
  init_number q_srv2prior           // Prior mean for catchability coefficient
  init_number cvq_srv2prior         // Prior CV for catchability coefficient
  number log_q2prior
  !! log_q2prior = log(q_srv2prior);

// ############### survey 3 and 4 use q from survey 1 and 2 respectively #############################
  init_number q_srv3prior           // Prior mean for catchability coefficient
  init_number cvq_srv3prior         // Prior CV for catchability coefficient
  number log_q3prior
  !! log_q3prior = log(q_srv3prior);
  
  init_number q_srv4prior           // Prior mean for catchability coefficient
  init_number cvq_srv4prior         // Prior CV for catchability coefficient
  number log_q4prior
  !! log_q4prior = log(q_srv4prior);
// ###########################################################################################################################

  init_number q_srv5prior           // Prior mean for catchability coefficient
  init_number cvq_srv5prior         // Prior CV for catchability coefficient
  number log_q5prior
  !! log_q5prior = log(q_srv5prior);
  
  init_number q_srv6prior           // Prior mean for catchability coefficient
  init_number cvq_srv6prior         // Prior CV for catchability coefficient
  number log_q6prior
  !! log_q6prior = log(q_srv6prior);
  
  init_number q_srv7prior           // Prior mean for catchability coefficient
  init_number cvq_srv7prior         // Prior CV for catchability coefficient
  number log_q7prior
  !! log_q7prior = log(q_srv7prior);
  
  init_number q_srv8prior           // Prior mean for catchability coefficient
  init_number cvq_srv8prior         // Prior CV for catchability coefficient
  number log_q8prior
  !! log_q8prior = log(q_srv8prior);
  
  init_int    yr_catchwt            // year catch-wt changes...... 
  init_number wt_ssqcatch_fish1           // Weight for catch estimation
  init_number wt_ssqcatch_fish3        // Weight for catch estimation
  init_number wt_srv1               // Weight for survey biomass estimation
  init_number wt_srv2               // Weight for survey biomass estimation
  init_number wt_srv3               // Weight for survey biomass estimation
  init_number wt_srv4               // Weight for survey biomass estimation
  init_number wt_srv5               // Weight for survey biomass estimation
  init_number wt_srv6               // Weight for survey biomass estimation
  init_number wt_srv7               // Weight for survey biomass estimation
  init_number wt_srv8               // Weight for survey biomass estimation
  init_number wt_fish1_age           // Weight for fishery age compositions
  init_number wt_srv1_age           // Weight for survey age compositions
  init_number wt_fish1_size          // Weight for fishery size compositions
  init_number wt_srv1_size          // Weight for survey size compostiions
  init_number wt_fish2_size          // Weight for fishery size compositions
  init_number wt_srv2_size          // Weight for survey size compostiions
  init_number wt_fish3_size          // Weight for fishery size compositions
  init_number wt_srv7_size          // Weight for survey size compostiions
  init_number wt_fish4_size          // Weight for fishery size compositions
  init_number wt_srv7_age          // Weight for survey size compostiions
  init_number wt_srv_extra_size          // Weight for fishery size compositions
  init_number wt_srv5_size          // Weight for survey size compostiions
  init_number wt_fish6_size          // Weight for fishery size compositions
  init_number wt_srv6_size          // Weight for survey size compostiions
  init_number wt_rec_var            // Weight for estimation recruitment variatiions penalty

 //############### for iterative reweighting runs
  init_number wt_fish1_age_iter           // Weight for fishery age compositions
  init_number wt_srv1_age_iter           // Weight for survey age compositions
  init_number wt_srv2_age_iter           // Weight for survey age compositions

  init_number wt_fish1_size_male_iter          // Weight for fishery size compositions
  init_number wt_fish1_size_female_iter          // Weight for fishery size compositions
  init_number wt_srv1_size_male_iter          // Weight for survey size compostiions
  init_number wt_srv1_size_female_iter          // Weight for survey size compostiions
  init_number wt_srv2_size_male_iter          // Weight for survey size compostiions
  init_number wt_srv2_size_female_iter          // Weight for survey size compostiions
  init_number wt_fish3_size_male_iter          // Weight for fishery size compositions
  init_number wt_fish3_size_female_iter          // Weight for fishery size compositions
  init_number wt_srv7_size_male_iter          // Weight for fishery size compositions
  init_number wt_srv7_size_female_iter          // Weight for fishery size compositions
 //########################################################
 
 // selectivity likelihoods do not appear to be added to the LL function so these do not have any effect
  init_number wt_sel_reg_fish1       // Weight on fishery selectivity regularity penalty
  init_number wt_sel_reg_fish2       // Weight on fishery selectivity regularity penalty
  init_number wt_sel_reg_fish3       // Weight on fishery selectivity regularity penalty
  init_number wt_sel_reg_fish4       // Weight on fishery selectivity regularity penalty
  init_number wt_sel_dome_fish1     // Weight on fishery selectivity dome-shape penalty   
  init_number wt_sel_dome_fish2     // Weight on fishery selectivity dome-shape penalty   
  init_number wt_sel_dome_fish3     // Weight on fishery selectivity dome-shape penalty   
  init_number wt_sel_dome_fish4     // Weight on fishery selectivity dome-shape penalty   
  init_number wt_sel_reg_srv1       // Weight on survey selectivity regularity penalty
  init_number wt_sel_reg_srv2       // Weight on survey selectivity regularity penalty
  init_number wt_sel_dome_srv1      // Weight on survey selectivity dome-shape penalty
  init_number wt_sel_dome_srv2      // Weight on survey selectivity dome-shape penalty
 //#################################################################################################
 
  init_number wt_fmort_reg          // Weight on fishing mortality regularity
  init_number wt_M_reg              // weight for penalizing Mdevs
  init_number wt_q_priors           // weight for catchability priors
  init_number wt_M_prior           // weight for catchability priors
  init_number wt_sigr_prior           // weight for catchability priors

  init_number hist_hal_prop;          // additional data for BS flag
  init_number yieldratio;
  
    !! cout<<yieldratio<<endl;

    !! cout<<"done with ctl"<<endl;
    !! ad_comm::change_datafile_name(data_file);    // Read data from the data file (i.e., tem_2020_na_wh.dat)
    
// Start and end years, recruitment age, number of age and length bins
  init_int      styr
  init_int      endyr
  init_int      recage
  init_int      nages
  init_int      nlenbins
  init_vector   len_bin_labels(1,nlenbins)

  int styr_rec
  int styr_sp
  int endyr_sp
  int endyr_rec
  int nyrs
  !!  nyrs = endyr - styr + 1;
  !!  styr_rec = (styr - nages) + 1;     // First year of recruitment
  !! styr_sp  = styr_rec - recage ;     // First year of spawning biomass  
  !! endyr_sp = endyr   - recage - 1;   // endyr year of (main) spawning biomass
  !! endyr_rec= endyr_rec_est;  // 
  vector yy(styr,endyr);
  !! yy.fill_seqadd(styr,1) ;
  vector aa(1,nages);
  !! aa.fill_seqadd(recage,1) ;

  init_number spawn_fract; // Spawning Month
  !! spawn_fract = (spawn_fract - 1) / 12;

// Natural mortality, proportion mature and weight at age
  init_vector   p_mature1(1,nages)
  init_vector   p_mature2(1,nages)
  init_vector   p_mature3(1,nages)
  init_vector   wt_m1(1,nages)
  init_vector   wt_f1(1,nages)
  init_vector   wt_m2(1,nages)
  init_vector   wt_f2(1,nages)
  init_vector   wt_all(1,nages)  //not used
  init_vector   wt_old(1,nages)  // not used
  init_vector prop_m(styr,endyr);  // proportion of males in RPN
  init_vector prop_m2(styr,endyr);  // proportion of males in RPN
// Observed catches
  init_vector   obs_catch_fish1(styr,endyr)
  init_vector   obs_catch_fish3(styr,endyr)


// Survey biomass estimates
  init_int      nyrs_srv1               // number of years of survey biomass estimates
  init_ivector  yrs_srv1(1,nyrs_srv1)         // years survey conducted in
  init_vector   obs_srv1_biom(1,nyrs_srv1)      // mean estimate of biomass
  init_vector   obs_srv1_se(1,nyrs_srv1)        // standard error of survey biomass estimates
  init_vector   obs_srv1_lci(1,nyrs_srv1)       // lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv1_uci(1,nyrs_srv1)       // upper confidence interval

  init_int      nyrs_srv2               // number of years of survey biomass estimates
  init_ivector  yrs_srv2(1,nyrs_srv2)         // years survey conducted in
  init_vector   obs_srv2_biom(1,nyrs_srv2)      // mean estimate of biomass
  init_vector   obs_srv2_se(1,nyrs_srv2)        // standard error of survey biomass estimates
  init_vector   obs_srv2_lci(1,nyrs_srv2)       // lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv2_uci(1,nyrs_srv2)       // upper confidence interval
 !! cout<<"Number of srv1 rpn is:"<<nyrs_srv2<<endl;

  init_int      nyrs_srv3               // number of years of survey biomass estimates
  init_ivector  yrs_srv3(1,nyrs_srv3)         // years survey conducted in
  init_vector   obs_srv3_biom(1,nyrs_srv3)      // mean estimate of biomass
  init_vector   obs_srv3_se(1,nyrs_srv3)        // standard error of survey biomass estimates
  init_vector   obs_srv3_lci(1,nyrs_srv3)       // lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv3_uci(1,nyrs_srv3)       // upper confidence interval

  init_int      nyrs_srv4               // number of years of survey biomass estimates
  init_ivector  yrs_srv4(1,nyrs_srv4)         // years survey conducted in
  init_vector   obs_srv4_biom(1,nyrs_srv4)      // mean estimate of biomass
  init_vector   obs_srv4_se(1,nyrs_srv4)        // standard error of survey biomass estimates
  init_vector   obs_srv4_lci(1,nyrs_srv4)       // lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv4_uci(1,nyrs_srv4)       // upper confidence interval
  
  init_int      nyrs_srv5               // number of years of survey biomass estimates
  init_ivector  yrs_srv5(1,nyrs_srv5)         // years survey conducted in
  init_vector   obs_srv5_biom(1,nyrs_srv5)      // mean estimate of biomass
  init_vector   obs_srv5_se(1,nyrs_srv5)        // standard error of survey biomass estimates
  init_vector   obs_srv5_lci(1,nyrs_srv5)       // lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv5_uci(1,nyrs_srv5)       // upper confidence interval

  init_int      nyrs_srv6               // number of years of survey biomass estimates
  init_ivector  yrs_srv6(1,nyrs_srv6)         // years survey conducted in
  init_vector   obs_srv6_biom(1,nyrs_srv6)      // mean estimate of biomass
  init_vector   obs_srv6_se(1,nyrs_srv6)        // standard error of survey biomass estimates
  init_vector   obs_srv6_lci(1,nyrs_srv6)       // lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv6_uci(1,nyrs_srv6)       // upper confidence interval

  init_int      nyrs_srv7               // number of years of survey biomass estimates
  init_ivector  yrs_srv7(1,nyrs_srv7)         // years survey conducted in
  init_vector   obs_srv7_biom(1,nyrs_srv7)      // mean estimate of biomass
  init_vector   obs_srv7_se(1,nyrs_srv7)        // standard error of survey biomass estimates
  init_vector   obs_srv7_lci(1,nyrs_srv7)       // lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv7_uci(1,nyrs_srv7)       // upper confidence interval
 !! cout<<"Number of trawl years:"<<nyrs_srv7<<endl;
 
  // Fishery age composition data
  init_int      nyrs_fish1_age            // number of years of fishery age compos
  init_vector  yrs_fish1_age(1,nyrs_fish1_age)    //the years of age comps
  init_vector   nsamples_fish1_age(1,nyrs_fish1_age)  // some measure of relative sample size for each age comp.
  init_matrix   oac_fish1(1,nyrs_fish1_age,1,nages)   // the actual year by year age comps
  number sdnr_fish1_age
// Survey age composition data
  init_int      nyrs_srv1_age             // number of years of survey age compositions
  init_ivector  yrs_srv1_age(1,nyrs_srv1_age)     // the years of survey age comps
  init_vector   nsamples_srv1_age(1,nyrs_srv1_age)  // some measure of relative sample size for each age comp.
  init_matrix   oac_srv1(1,nyrs_srv1_age,1,nages)   // the year by year age survey age comps
  number sdnr_srv1_age
  init_int      nyrs_srv2_age             // number of years of survey age compositions
  init_ivector  yrs_srv2_age(1,nyrs_srv2_age)     // the years of survey age comps
  init_vector   nsamples_srv2_age(1,nyrs_srv2_age)  // some measure of relative sample size for each age comp.
  init_matrix   oac_srv2(1,nyrs_srv2_age,1,nages)   // the year by year age survey age comps
  number sdnr_srv2_age
  !! cout<<"Number of ages srv2:"<<nyrs_srv2_age<<endl;
  
// Fishery size composition data
  init_int      nyrs_fish1_size             // number of years of fish1ery size comps
  init_ivector  yrs_fish1_size(1,nyrs_fish1_size)   // the years of fish1ery size comps
  init_vector   nsamples_fish1_size(1,nyrs_fish1_size)// some measure of relative sample size for each fish1ery comp
  init_matrix   osc_fish1_m(1,nyrs_fish1_size,1,nlenbins)// year by year fishery size comps
  init_matrix   osc_fish1_f(1,nyrs_fish1_size,1,nlenbins)// year by year fishery size comps
  number sdnr_fish1_size
  
  init_int      nyrs_fish2_size             // number of years of fish2ery size comps
  init_ivector  yrs_fish2_size(1,nyrs_fish2_size)   // the years of fishery size comps
  init_vector   nsamples_fish2_size(1,nyrs_fish2_size)// some measure of relative sample size for each fish2ery comp
  init_matrix   osc_fish2(1,nyrs_fish2_size,1,nlenbins)// year by year fishery size comps
   number sdnr_fish2_size
  !! cout<<"Number of fish 2 sizes:"<< nyrs_fish2_size <<endl;
  
   
  init_int      nyrs_fish3_size             // number of years of fish3ery size comps
  init_ivector  yrs_fish3_size(1,nyrs_fish3_size)   // the years of fish3ery size comps
  init_vector   nsamples_fish3_size(1,nyrs_fish3_size)// some measure of relative sample size for each fish3ery comp
  init_matrix   osc_fish3_m(1,nyrs_fish3_size,1,nlenbins)// year by year fishery size comps
  init_matrix   osc_fish3_f(1,nyrs_fish3_size,1,nlenbins)// year by year fishery size comps
  number sdnr_fish3_size
  !! cout<<"Number of fish 3 sizes:"<< nyrs_fish3_size <<endl;
  !! cout<<"Years fish 4 sizes:"<< yrs_fish3_size <<endl;

  
  init_int      nyrs_fish4_size             // number of years of fish3ery size comps
  init_ivector  yrs_fish4_size(1,nyrs_fish4_size)   // the years of fishery size comps
  init_vector   nsamples_fish4_size(1,nyrs_fish4_size)// some measure of relative sample size for each fish3ery comp
  init_matrix   osc_fish4(1,nyrs_fish4_size,1,nlenbins)// year by year fishery size comps
  number sdnr_fish4_size
 
  !! cout<<"Number of fish 4 sizes:"<< nyrs_fish4_size <<endl;
// Survey size composition data
  init_int      nyrs_srv1_size            // number of years of survey size comps
  init_ivector  yrs_srv1_size(1,nyrs_srv1_size)   // the years of survey size comps
  init_vector   nsamples_srv1_size(1,nyrs_srv1_size)//some measure of relative sample size for each survey size comp
  init_matrix   osc_srv1_m(1,nyrs_srv1_size,1,nlenbins)//year by year size comps
  init_matrix   osc_srv1_f(1,nyrs_srv1_size,1,nlenbins)//year by year size comps
  number sdnr_srv1_size
 
   !! cout<<"Number of srv1 size is:"<<nyrs_srv1_size<<endl;

    init_int      nyrs_srv2_size            // number of years of survey size comps
  init_ivector  yrs_srv2_size(1,nyrs_srv2_size)   // the years of survey size comps
  init_vector   nsamples_srv2_size(1,nyrs_srv2_size)//some measure of relative sample size for each survey size comp
  init_matrix   osc_srv2_m(1,nyrs_srv2_size,1,nlenbins)//year by year size comps
   init_matrix   osc_srv2_f(1,nyrs_srv2_size,1,nlenbins)//year by year size comps
  number sdnr_srv2_size
 
    init_int      nyrs_srv7_size            // number of years of survey size comps
  init_ivector  yrs_srv7_size(1,nyrs_srv7_size)   // the years of survey size comps
  init_vector   nsamples_srv7_size(1,nyrs_srv7_size)//some measure of relative sample size for each survey size comp
  init_matrix   osc_srv7_m(1,nyrs_srv7_size,1,nlenbins)//year by year size comps
  init_matrix   osc_srv7_f(1,nyrs_srv7_size,1,nlenbins)//year by year size comps
 number sdnr_srv7_size
 
    init_int      nyrs_srv7_age             // number of years of survey age compositions
  init_ivector  yrs_srv7_age(1,nyrs_srv7_age)     // the years of survey age comps
  init_vector   nsamples_srv7_age(1,nyrs_srv7_age)  // some measure of relative sample size for each age comp.
  init_matrix   oac_srv7(1,nyrs_srv7_age,1,nages)   // the year by year age survey age comps
 number sdnr_srv7_age
 
  !! cout<<"Number of survey age years is "<<nyrs_srv7_age<<endl;

  int phase_selcoff_fsh1 
  int phase_selcoff_fsh2 
  int phase_selcoff_fsh3 
  int phase_selcoff_fsh4
  int phase_selcoff_fsh5
  int phase_logist_fsh1
  int phase_logist_fsh2
  int phase_logist_fsh3
  int phase_logist_fsh4
  int phase_logist_fsh5
  int phase_dlogist_fsh1
  int phase_dlogist_fsh2
  int phase_dlogist_fsh3
  int phase_dlogist_fsh4
  int phase_dlogist_fsh5
  int phase_selcoff_srv1 
  int phase_selcoff_srv2 
  int phase_selcoff_srv7
  int phase_selcoff_srv10
  int phase_logist_srv1 
  int phase_logist_srv2 
  int phase_logist_srv7
  int phase_logist_srv10
  int phase_dlogist_srv1 
  int phase_dlogist_srv2 
  int phase_dlogist_srv7
  int phase_dlogist_srv10

// Size-age transition matrix:  proportion at size given age
  init_matrix   sizeage1_m(1,nages,1,nlenbins)     //size comp #1
  init_matrix   sizeage1_f(1,nages,1,nlenbins)     //lets you add another size-age matrix, remove this here and in Get_Predicted section to use only one, or just have two identical
  init_matrix   sizeage2_m(1,nages,1,nlenbins)     //size comp #1
  init_matrix   sizeage2_f(1,nages,1,nlenbins)     //lets you add another size-age matrix, remove this here and in Get_Predicted section to use only one, or just have two identical
  init_matrix   sizeage3_m(1,nages,1,nlenbins)     //size comp #1
  init_matrix   sizeage3_f(1,nages,1,nlenbins)     //lets you add another size-age matrix, remove this here and in Get_Predicted section to use only one, or just have two identical
  init_matrix   sizeage_all(1,nages,1,nlenbins)     //lets you add another size-age matrix, remove this here and in Get_Predicted section to use only one, or just have two identical

// Ageing error transition matrix:  proportion at reader age given true age
  init_matrix   ageage(1,nages,1,nages)       // ageing error matrix
  init_number eof

     !! cout<<"The universal answer is "<<eof<<endl;

// Initialize some counting variables
  int iyr
  int i
  int j
 int ii
 int l
 int y
  vector offset(1,16);                                    // Multinomial "offset"

 LOCAL_CALCS


  if(rec_like_type>0) styr_rec=styr-nages+1;
   if(wt_rec_var==0)                
   {
     if (ph_sigr>0)                 
     {
       cout << "Warning, wt_rec_var is zero, so can't estimate sigr!@"<<endl;
       cout << "turning sigr off "<<endl;
       ph_sigr =-4;
       cout << "hit any key, then enter to continue"<<endl;
       char  xxx;
       cin >> xxx;
     }
   }

  // define phases for selectivity
  //////////////////////////////////////////////////////////////////////////////////////////////
  // NOTE THAT MANY OF THESE ARE HARDWIRED IN THE PAR SECTION WHEN SEL PAR DEFINED
  /////////////////////////////////////////////////////////////////////////////////////
   switch (srv1_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_srv1 = ph_srv1_sel;
       phase_logist_srv1  = -1;
       phase_dlogist_srv1 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_srv1 = -1; 
       phase_logist_srv1  = ph_srv1_sel;
       phase_dlogist_srv1 = -1;
      }
     break;
     case 3 : // Double logistic
     {
       phase_selcoff_srv1 = -1; 
       phase_logist_srv1  = ph_srv1_sel;
       phase_dlogist_srv1 = ph_srv1_sel;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_srv1 = -1; 
       phase_logist_srv1  = ph_srv1_sel;
       phase_dlogist_srv1 = ph_srv1_sel;
     }
     break;
   }
// -----------
   
   switch (srv2_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_srv2 = ph_srv1_sel;
       phase_logist_srv2  = -1;
       phase_dlogist_srv2 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_srv2 = -1; 
       phase_logist_srv2  = ph_srv1_sel;
       phase_dlogist_srv2 = -1;
      }
     break;
     case 3 : // Double logistic
     {
       phase_selcoff_srv2 = -1; 
       phase_logist_srv2  = ph_srv1_sel;
       phase_dlogist_srv2 = ph_srv1_sel;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_srv2 = -1; 
       phase_logist_srv2  = ph_srv1_sel;
       phase_dlogist_srv2 = ph_srv1_sel;
     }
     break;
   }
// -----------
   
   switch (srv7_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_srv7 = ph_srv1_sel;
       phase_logist_srv7  = -1;
       phase_dlogist_srv7 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_srv7 = -1; 
       phase_logist_srv7  = ph_srv1_sel;
       phase_dlogist_srv7 = -1;
      }
     break;
     case 3 : // POWER FUNCTION
     {
       phase_selcoff_srv7 = -1; 
       phase_logist_srv7  = ph_srv1_sel;
       phase_dlogist_srv7 = -1;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_srv7 = -1; 
       phase_logist_srv7  = ph_srv1_sel;
       phase_dlogist_srv7 = ph_srv1_sel;
     }
     break;
     case 5 : // Exponential logistic
     {
       phase_selcoff_srv7 = -1; 
       phase_logist_srv7  = -1;
       phase_dlogist_srv7 = -1;
     }
     break;
    if(wt_srv7==0) {
       phase_selcoff_srv7=-1;
       phase_logist_srv7=-1;
       phase_dlogist_srv7 = -1; }
     }

   switch (srv10_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_srv10 = ph_LL_block2;
       phase_logist_srv10  = -1;
       phase_dlogist_srv10 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_srv10 = -1; 
       phase_logist_srv10  = ph_LL_block2;
       phase_dlogist_srv10 = -1;
      }
     break;
     case 3 : // Double logistic
     {
       phase_selcoff_srv10 = -1; 
       phase_logist_srv10  = ph_LL_block2;
       phase_dlogist_srv10 = ph_LL_block2;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_srv10 = -1; 
       phase_logist_srv10  = ph_LL_block2;
       phase_dlogist_srv10 = ph_LL_block2;
     }
     break;
   }
// -----------
// -----------
   
   switch (fsh1_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_fsh1 = ph_fish_sel;
       phase_logist_fsh1  = -1;
       phase_dlogist_fsh1 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_fsh1 = -1; 
       phase_logist_fsh1  = ph_fish_sel;
       phase_dlogist_fsh1 = -1;
      }
     break;
     case 3 : // Double logistic
     {
       phase_selcoff_fsh1 = -1; 
       phase_logist_fsh1  = ph_fish_sel;
       phase_dlogist_fsh1 = ph_fish_sel;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_fsh1 = -1; 
       phase_logist_fsh1  = ph_fish_sel;
       phase_dlogist_fsh1 = ph_fish_sel;
     }
     break;
     case 5 : // Exponential logistic
     {
       phase_selcoff_fsh1 = -1; 
       phase_logist_fsh1  = -1;
       phase_dlogist_fsh1 = -1;
     }
     break;
   }
// -----------
   
   switch (fsh2_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_fsh2 = ph_fish_sel;
       phase_logist_fsh2  = -1;
       phase_dlogist_fsh2 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_fsh2 = -1; 
       phase_logist_fsh2  = ph_fish_sel;
       phase_dlogist_fsh2 = -1;
      }
     break;
     case 3 : // Double logistic
     {
       phase_selcoff_fsh2 = -1; 
       phase_logist_fsh2  = ph_fish_sel;
       phase_dlogist_fsh2 = ph_fish_sel;
     }
     break;
      case 5 : // Exponential logistic
     {
       phase_selcoff_fsh2 = -1; 
       phase_logist_fsh2  = -1;
       phase_dlogist_fsh2 = -1;
     }
  }// -----------
   
   switch (fsh3_sel_opt) 
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_fsh3 = ph_fish_sel;
       phase_logist_fsh3  = -1;
       phase_dlogist_fsh3 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_fsh3 = -1; 
       phase_logist_fsh3  = ph_fish_sel;
       phase_dlogist_fsh3 = -1;
      }
     break;
     case 3 : // GAMMA FUNCTION
     {
       phase_selcoff_fsh3 = -1; 
       phase_logist_fsh3  = ph_fish_sel;
       phase_dlogist_fsh3 = -1;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_fsh3 = -1; 
       phase_logist_fsh3  = ph_fish_sel;
       phase_dlogist_fsh3 = ph_fish_sel;
     }
     break;
   }
//----
   switch (fsh4_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_fsh4 = ph_fish_sel;
       phase_logist_fsh4  = -1;
       phase_dlogist_fsh4 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_fsh4 = -1; 
       phase_logist_fsh4  = ph_fish_sel;
       phase_dlogist_fsh4 = -1;
      }
     break;
     case 3 : // Double logistic
     {
       phase_selcoff_fsh4 = -1; 
       phase_logist_fsh4  = ph_fish_sel;
       phase_dlogist_fsh4 = ph_fish_sel;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_fsh4 = -1; 
       phase_logist_fsh4  = ph_fish_sel;
       phase_dlogist_fsh4 = ph_fish_sel;
     }
    if(ph_ifq<1) {
       phase_selcoff_fsh4=-1;
       phase_logist_fsh4 = -1;
       phase_dlogist_fsh4 = -1; }
     break;
   }
   
//----

   switch (fsh5_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_fsh5 = ph_ifq_block2;
       phase_logist_fsh5  = -1;
       phase_dlogist_fsh5 = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_fsh5 = -1; 
       phase_logist_fsh5  = ph_ifq_block2;
       phase_dlogist_fsh5 = -1;
      }
     break;
     case 3 : // Double logistic
     {
       phase_selcoff_fsh5 = -1; 
       phase_logist_fsh5  = ph_ifq_block2;
       phase_dlogist_fsh5 = ph_ifq_block2;
     }
     break;
     case 4 : // Exponential logistic
     {
       phase_selcoff_fsh5 = -1; 
       phase_logist_fsh5  = ph_ifq_block2;
       phase_dlogist_fsh5 = ph_ifq_block2;
     }
    if(ph_ifq_block2<1) {
       phase_selcoff_fsh5=-1;
       phase_logist_fsh5 = -1;
       phase_dlogist_fsh5 = -1; }
     break;
   }
   
// -----------
// Calculate "offset" for multinomials - survey age, fishery size, survey size
//   "Offset" value lets the multinomial likelihood equal zero when the observed and
//     predicted are equal as in Fournier (1990) "robustifies"
//   First step is to ensure that the data are expressed as proportions
  for (i=1; i<=nyrs_fish1_age; i++)
  {
   oac_fish1(i)/=sum(oac_fish1(i));
   offset(1) -= nsamples_fish1_age(i) *((oac_fish1(i) + 0.001)*log(oac_fish1(i) + 0.001)); 
  }

  for (i=1; i<=nyrs_srv1_age; i++)
  {
   oac_srv1(i)/=sum(oac_srv1(i));
   offset(2) -= nsamples_srv1_age(i)*((oac_srv1(i) + 0.001)*log(oac_srv1(i) + 0.001));
  }

  for (i=1; i<=nyrs_fish1_size; i++)
  {
   osc_fish1_f(i)/=sum(osc_fish1_f(i));
   offset(3) -= nsamples_fish1_size(i)*((osc_fish1_f(i) + 0.001)*log(osc_fish1_f(i) + 0.001));
   osc_fish1_m(i)/=sum(osc_fish1_m(i));
   offset(4) -= nsamples_fish1_size(i)*((osc_fish1_m(i) + 0.001)*log(osc_fish1_m(i) + 0.001));
  }
  for (i=1; i<=nyrs_fish2_size; i++)
  {
   osc_fish2(i)/=sum(osc_fish2(i));
   offset(5) -= nsamples_fish2_size(i)*((osc_fish2(i) + 0.001)*log(osc_fish2(i) + 0.001));
  }
  for (i=1; i<=nyrs_fish3_size; i++)
  {
   osc_fish3_f(i)/=sum(osc_fish3_f(i));
   offset(6) -= nsamples_fish3_size(i)*((osc_fish3_f(i) + 0.001)*log(osc_fish3_f(i) + 0.001));
   osc_fish3_m(i)/=sum(osc_fish3_m(i));
   offset(7) -= nsamples_fish3_size(i)*((osc_fish3_m(i) + 0.001)*log(osc_fish3_m(i) + 0.001));

     }
  for (i=1; i<=nyrs_fish4_size; i++)
  {
   osc_fish4(i)/=sum(osc_fish4(i));
   offset(8) -= nsamples_fish4_size(i)*((osc_fish4(i) + 0.001)*log(osc_fish4(i) + 0.001));
  }

  for (i=1; i<=nyrs_srv1_size; i++)
  {
   osc_srv1_f(i)/=sum(osc_srv1_f(i));
   offset(9) -= nsamples_srv1_size(i)*((osc_srv1_f(i) + 0.001)*log(osc_srv1_f(i) + 0.001));
    osc_srv1_m(i)/=sum(osc_srv1_m(i));
   offset(10) -= nsamples_srv1_size(i)*((osc_srv1_m(i) + 0.001)*log(osc_srv1_m(i) + 0.001));
  }
  for (i=1; i<=nyrs_srv2_size; i++)
  {
   osc_srv2_f(i)/=sum(osc_srv2_f(i));
   offset(11) -= nsamples_srv2_size(i)*((osc_srv2_f(i) + 0.001)*log(osc_srv2_f(i) + 0.001));
    osc_srv2_m(i)/=sum(osc_srv2_m(i));
   offset(12) -= nsamples_srv2_size(i)*((osc_srv2_m(i) + 0.001)*log(osc_srv2_m(i) + 0.001));
  }
  for (i=1; i<=nyrs_srv7_size; i++)
  {
   osc_srv7_f(i)/=sum(osc_srv7_f(i));
   offset(13) -= nsamples_srv7_size(i)*((osc_srv7_f(i) + 0.001)*log(osc_srv7_f(i) + 0.001));
    osc_srv7_m(i)/=sum(osc_srv7_m(i));
   offset(14) -= nsamples_srv7_size(i)*((osc_srv7_m(i) + 0.001)*log(osc_srv7_m(i) + 0.001));
  }
  for (i=1; i<=nyrs_srv7_age; i++)
  {
   oac_srv7(i)/=sum(oac_srv7(i));
   offset(15) -= nsamples_srv7_age(i)*((oac_srv7(i) + 0.001)*log(oac_srv7(i) + 0.001)); }
  for (i=1; i<=nyrs_srv2_age; i++)
  {
   oac_srv2(i)/=sum(oac_srv2(i));
   offset(16) -= nsamples_srv2_age(i)*((oac_srv2(i) + 0.001)*log(oac_srv2(i) + 0.001));
  }

 END_CALCS
 
INITIALIZATION_SECTION

 log_mean_rec        2.0
 log_sigr           sigrprior
 logm               -2.0
 log_avg_F_fish1    -1.0
 log_avg_F_fish3    -1.0
 log_mF50           -2.0
 log_mF40           -2.0
 log_mF35           -2.0
 log_q_srv1          0.5
 log_q_srv2          0.5
 log_q_srv3          0.5
 log_q_srv4          0.5
 log_q_srv5          0.5
 log_q_srv6          0.5
 log_q_srv7          0.5
 log_q_srv8         -1.0
 log_q_LL_fish_recent   -1.0

 log_a50_fish1_f      1.5 
 log_a50_fish1_m      1.5
 log_a50_fish2        1.5
 log_a50_fish3_m      1.0
 log_a50_fish3_f      1.0
 log_a50_fish4_m      1.5
 log_a50_fish4_f      1.5
 log_a50_fish5_m      1.5
 log_a50_fish5_f      1.5

 log_a50_srv1_f       1.5
 log_a50_srv1_m       1.5
 log_a50_srv2_f       1.5
 log_a50_srv2_m       1.5
 log_a50_srv7_f       1.0
 log_a50_srv7_m       1.0
 log_a50_srv10_f      1.5
 log_a50_srv10_m      1.5

 log_delta_fish1_f    1.5
 log_delta_fish1_m    1.5
 log_delta_fish2      1.5
 log_delta_fish3_f    1.5
 log_delta_fish3_m    1.5
 log_delta_fish4_f    1.5
 log_delta_fish4_m    1.5
 log_delta_fish5_f    1.5
 log_delta_fish5_m    1.5

 log_delta_srv1_f     1.5
 log_delta_srv1_m     1.5
 log_delta_srv2_f     1.5
 log_delta_srv2_m     1.5
 log_delta_srv10_f    1.5
 log_delta_srv10_m    1.5
 
 PARAMETER_SECTION

  // Key parameters
 // ########################################################################################################
 
 // Stock-recruitment
// ########################################################################################################
 // steepness and R0 not estimated in base model, par have neg phases
  init_bounded_number  steepness(0.2001,0.999,ph_steepness)   // Stock recruitment steepness (not estimated)
  init_bounded_number  log_Rzero(1,5,ph_Rzero);                     // Unfish equil recruitment (logged) (not estimated)
// ########################################################################################################
  init_bounded_dev_vector  log_rec_dev(styr-nages+2,endyr_rec_est,-10,10,ph_recdev);  // Recruitment deviations from before the asssessment starts to present
  init_bounded_number      log_mean_rec(-1,10,ph_mean_rec);          // Unfish equil recruitment (logged) (estimated)
  init_bounded_number      log_sigr(-2,1,ph_sigr);           // Recruitment sdev parameter  (not estimated)

 // Parameters for computing SPR rates
// ########################################################################################################
  init_bounded_number   log_mF50(-5,0.5,ph_F50)            // Estimated F50 (phase hardwired in tpl data section to 6)
  init_bounded_number   log_mF40(-5,0.5,ph_F50)            // Estimated F40 (phase hardwired in tpl data section to 6)
  init_bounded_number   log_mF35(-5,0.5,ph_F50)            // Estimated F35 (phase hardwired in tpl data section to 6)

 // Natural Mortality
//###########################################################################################################
  init_bounded_number   logm(-4,0.5,ph_m);                            // Estimate log natural mortality (estimated)
  init_bounded_vector   log_M_devs(styr,endyr,-5,5,ph_Mdevs);  // Annual natural mortality deviations
  init_bounded_vector   log_M_devs_age(1,nages,-5,5,ph_Mdevs_age);  // Annual natural mortality deviations
  init_bounded_number   mdelta(-0.05,0.05,ph_mdelta);          // Use parameter to have a male natural mortality (not estimated)
  number                natmort;                               //natural mortality when constant
  vector                M_year(styr,endyr);                    // Natural mortality by year
  vector                M_age(1,nages);                    // Natural mortality by year
  matrix                M(styr,endyr,1,nages);                 // Natural mortality by year and age
  vector                natmortv(1,nages);                     // Create a vector of natural mortalities for proj.dat

  // Fishing mortality
 // ########################################################################################################
  init_bounded_number       log_avg_F_fish1(-10,2,ph_avg_F);                      // Log average fishing mortality
  init_bounded_dev_vector   log_F_devs_fish1(styr,endyr,-5,5,ph_Fdev);  // Annual fishing mortality deviations
  init_bounded_number       log_avg_F_fish3(-10,2,ph_avg_F);                      // Log average fishing mortality
  init_bounded_dev_vector   log_F_devs_fish3(styr+3,endyr,-5,5,ph_Fdev);  // Annual fishing mortality deviations
  vector                Fmort_fish1(styr,endyr);                        // Fishing mortality by year
  vector                Fmort_fish3(styr,endyr);                        // Fishing mortality by year
  matrix                F_fish1_f(styr,endyr,1,nages);                    // Fishing mortality by year and age
  matrix                F_fish1_m(styr,endyr,1,nages);                    // Fishing mortality by year and age
  matrix                F_fish3_f(styr,endyr,1,nages);                    // Fishing mortality by year and age
  matrix                F_fish3_m(styr,endyr,1,nages);                    // Fishing mortality by year and age
  number                hist_hal_F; // Option of adding historical proportion of current average hook and line catch
  
 // Catchability
// ########################################################################################################
 //    srv1 and srv2 (NPWs) are not fit, but catchability is used by srvy 3 and srvy 4 (RPNs) that are fit in model
 //    qsrvy 9 is only estimated if qsrvy2 is not estimated 
  init_bounded_number           log_q_srv1(-5,5,ph_q_srv1);           // Estimate Log survey catchability for DOM LL SRV RPN/RPW        (estimated)
  init_bounded_number           log_q_srv2(-5,5,ph_q_srv2);           // Estimate Log survey catchability for JPN LL SRV RPN/RPW        (estimated)
  init_bounded_number           log_q_srv3(-5,5,ph_q_srv3);                  // Estimate Log survey catchability use same value as q1          (not estimated)
  init_bounded_number           log_q_srv4(-5,5,ph_q_srv4);                  // Estimate Log survey catchability use same value as q2          (not estimated)
  init_bounded_number           log_q_srv5(-5,5,ph_q_srv5);           // Estimate Log survey catchability DOM LL fishery RPW            (estimated)
  init_bounded_number           log_q_srv6(-5,5,ph_q_srv6);           // Estimate Log survey catchability JPN LL fishery RPW            (estimated)
  init_bounded_number           log_q_srv7(-5,5,ph_q_srv7);           // Estimate Log survey catchability DOM trawl survey              (estimated)
  init_bounded_number           log_q_srv8(-15,5,ph_q_srv8);           // Estimate Log survey catchability DOM LL fishery RPW post-IFQ   (estimated)
  init_bounded_number           log_q_srv9(-5,5,ph_srv2_q2);         // Estimate Log survey catchability JPN LL SRV RPN/RPW 1990-1994  (set to qsrv2 unless qsrv2 not estimated)
  init_bounded_number           log_q_LL_srvy_recent(-5,5,ph_q_LL_srv_rec);  // estimate log survey catchability for DOM LL SRVY RPN recent time block
  init_bounded_number           log_q_LL_fish_recent(-15,5,ph_q_IFQ_rec);  // estimate log fishery catchability for DOM LL fishery RPW recent time block

  // Fishery selectivity
 // ########################################################################################################
    // note hardwiring of many phases and repeat names for alternate sel forms (e.g., power and gamma functions for trawl data use parameters named from logistic selectivity function)
    // NEED TO CLEAN THIS UP BIG TIME
    
  init_vector       log_fish1_sel_coffs_f(1,n_fish_sel_ages,phase_selcoff_fsh1); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish1_f(-1,4,ph_fish_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish1_f(-5,4,ph_fish_sel_delt);                 // age between 50% selection and 95% selection....
  number          a50_fish1_f;                    // age at 50% selection                                                   
  init_number       d50_fish1_f(phase_dlogist_fsh1);
  number        delta_fish1_f;                    // age between 50% selection and 95% selection....
  init_number       gamma_fish1_f(phase_dlogist_fsh1);
  vector        log_fish1_sel_f(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish1_sel_f(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish1sel_f;                               // average fishery selectivity

  init_vector       log_fish1_sel_coffs_m(1,n_fish_sel_ages,phase_selcoff_fsh1); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish1_m(-1,4,ph_fish_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish1_m(-5,4,ph_fish_sel_delt_alt);                 // age between 50% selection and 95% selection....
  number          a50_fish1_m;                    // age at 50% selection                                                   
  init_number       d50_fish1_m(phase_dlogist_fsh1);
  number        delta_fish1_m;                    // age between 50% selection and 95% selection....
  init_number       gamma_fish1_m(phase_dlogist_fsh1);
  vector        log_fish1_sel_m(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish1_sel_m(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish1sel_m;                               // average fishery selectivity

  init_vector       log_fish2_sel_coffs(1,n_fish_sel_ages,phase_selcoff_fsh2); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number       log_a50_fish2(-1,2.5,ph_fish2_sel);                 // age at 50% selection                                                   
  init_bounded_number       log_delta_fish2(-5,4,ph_fish_sel_delt_alt);                 // age between 50% selection and 95% selection....
  number          a50_fish2;                    // age at 50% selection                                                   
  number        delta_fish2;                    // age between 50% selection and 95% selection....
  vector        log_fish2_sel(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish2_sel(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish2sel;                               // average fishery selectivity

  init_vector       log_fish3_sel_coffs_f(1,n_fish_sel_ages,phase_selcoff_fsh3); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish3_f(-1,4,ph_fish_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish3_f(-5,4,ph_fish_sel_delt);                 // age between 50% selection and 95% selection....
  number          a50_fish3_f;                    // age at 50% selection                                                   
  init_number       d50_fish3_f(-1);
  number        delta_fish3_f;                    // age between 50% selection and 95% selection....
  init_bounded_number log_gamma_fish3_f(-0.00001,-0.00000000001,-4);
  vector        log_fish3_sel_f(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish3_sel_f(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish3sel_f;                               // average fishery selectivity
  number        gamma_fish3_f;

  init_vector       log_fish3_sel_coffs_m(1,n_fish_sel_ages,phase_selcoff_fsh3); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish3_m(-1,4,ph_fish_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish3_m(-5,4,ph_fish_sel_delt_alt);                 // age between 50% selection and 95% selection....
  number          a50_fish3_m;                    // age at 50% selection                                                   
  init_number       d50_fish3_m(-1);
  number        delta_fish3_m;                    // age between 50% selection and 95% selection....
  init_bounded_number log_gamma_fish3_m(-0.18,-0.173,-4);
  vector        log_fish3_sel_m(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish3_sel_m(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish3sel_m;                               // average fishery selectivity
  number          gamma_fish3_m;

  init_vector       log_fish4_sel_coffs_f(1,n_fish_sel_ages,phase_selcoff_fsh4); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish4_f(-1,4,ph_ifq);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish4_f(-5,4,ph_fish_sel_delt);                 // age between 50% selection and 95% selection....
  number          a50_fish4_f;                    // age at 50% selection                                                   
  init_number       d50_fish4_f(phase_dlogist_fsh4);
  number        delta_fish4_f;                    // age between 50% selection and 95% selection....
  init_number       gamma_fish4_f(phase_dlogist_fsh4);
  vector        log_fish4_sel_f(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish4_sel_f(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish4sel_f;                               // average fishery selectivity

  init_vector       log_fish4_sel_coffs_m(1,n_fish_sel_ages,phase_selcoff_fsh4); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish4_m(-1,4,ph_ifq);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish4_m(-5,4,ph_fish_sel_delt);                 // age between 50% selection and 95% selection....
  number          a50_fish4_m;                    // age at 50% selection                                                   
  init_number       d50_fish4_m(phase_dlogist_fsh4);
  number        delta_fish4_m;                    // age between 50% selection and 95% selection....
  init_number       gamma_fish4_m(-4);
  vector        log_fish4_sel_m(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish4_sel_m(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish4sel_m;                               // average fishery selectivity


  init_vector       log_fish5_sel_coffs_f(1,n_fish_sel_ages,phase_selcoff_fsh5); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish5_f(-1,4,ph_ifq_block2);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish5_f(-5,4,ph_ifq_block2);                 // age between 50% selection and 95% selection....
  number          a50_fish5_f;                    // age at 50% selection                                                   
  init_number       d50_fish5_f(phase_dlogist_fsh5);
  number        delta_fish5_f;                    // age between 50% selection and 95% selection....
  init_number       gamma_fish5_f(phase_dlogist_fsh5);
  vector        log_fish5_sel_f(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish5_sel_f(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish5sel_f;                               // average fishery selectivity

  init_vector       log_fish5_sel_coffs_m(1,n_fish_sel_ages,phase_selcoff_fsh5); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_fish5_m(-1,4,ph_ifq_block2);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_fish5_m(-5,4,ph_ifq_block2);                 // age between 50% selection and 95% selection....
  number          a50_fish5_m;                    // age at 50% selection                                                   
  init_number       d50_fish5_m(phase_dlogist_fsh5);
  number        delta_fish5_m;                    // age between 50% selection and 95% selection....
  init_number       gamma_fish5_m(-4);
  vector        log_fish5_sel_m(1,nages);                      // vector of fishery selectivy log parameters including those not estimated
  vector        fish5_sel_m(1,nages);                          // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgfish5sel_m;                               // average fishery selectivity

  vector        sel_rep_proj_f(1,nages);                      //for reporting single sex rescaled selectivity from fish1   
  vector        sel_rep_proj_m(1,nages);                      //for reporting single sex rescaled selectivity from fish1   

// ########################################################################################################

 // Survey selectivities
// ########################################################################################################

  init_vector       log_srv1_sel_coffs_f(1,n_srv1_sel_ages,phase_selcoff_srv1);   // vector of survey selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv1_f(-1,4,ph_srv1_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv1_f(-5,4,ph_srv_sel_delt);                 // age between 50% selection and 95% selection....
  number          a50_srv1_f;                     // age at 50% selection                                                   
  init_number       d50_srv1_f(phase_dlogist_srv1);
  number        delta_srv1_f;                     // age between 50% selection and 95% selection....
  init_number       gamma_srv1_f(phase_dlogist_srv1);
  vector        log_srv1_sel_f(1,nages);              // vector of survey selectivy log parameters including those not estimated
  vector        srv1_sel_f(1,nages);                  // vectory of survey selectivty parameters on arithmetic scale
  number        log_avgsrv1sel_f;                     // average survey selectivity

  init_vector       log_srv1_sel_coffs_m(1,n_srv1_sel_ages,phase_selcoff_srv1);   // vector of survey selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv1_m(-1,4,ph_srv1_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv1_m(-5,4,ph_srv_sel_delt);                 // age between 50% selection and 95% selection....
  number          a50_srv1_m;                     // age at 50% selection                                                   
  init_number       d50_srv1_m(phase_dlogist_srv1);
  number        delta_srv1_m;                     // age between 50% selection and 95% selection....
  init_number       gamma_srv1_m(phase_dlogist_srv1);
  vector        log_srv1_sel_m(1,nages);              // vector of survey selectivy log parameters including those not estimated
  vector        srv1_sel_m(1,nages);                  // vectory of survey selectivty parameters on arithmetic scale
  number        log_avgsrv1sel_m;                     // average survey selectivity

  init_vector       log_srv2_sel_coffs_f(1,n_srv1_sel_ages,phase_selcoff_srv2);   // vector of survey selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv2_f(-1,4,ph_srv2_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv2_f(-5,4,ph_srv_sel_delt_alt);                 // age between 50% selection and 95% selection....
  number          a50_srv2_f;                     // age at 50% selection                                                   
  init_number       d50_srv2_f(phase_dlogist_srv2);
  init_number       gamma_srv2_f(phase_dlogist_srv2);
  number        delta_srv2_f;                     // age between 50% selection and 95% selection....
  vector        log_srv2_sel_f(1,nages);              // vector of survey selectivy log parameters including those not estimated
  vector        srv2_sel_f(1,nages);                  // vectory of survey selectivty parameters on arithmetic scale
  number        log_avgsrv2sel_f;                     // average survey selectivity

  init_vector       log_srv2_sel_coffs_m(1,n_srv1_sel_ages,phase_selcoff_srv2);   // vector of survey selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv2_m(-1,4,ph_srv2_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv2_m(-5,4,ph_srv_sel_delt_alt);                 // age between 50% selection and 95% selection....
  number          a50_srv2_m;                     // age at 50% selection                                                   
  init_number       d50_srv2_m(phase_dlogist_srv2);
  init_number       gamma_srv2_m(phase_dlogist_srv2);
  number        delta_srv2_m;                     // age between 50% selection and 95% selection....
  vector        log_srv2_sel_m(1,nages);              // vector of survey selectivy log parameters including those not estimated
  vector        srv2_sel_m(1,nages);                  // vectory of survey selectivty parameters on arithmetic scale
  number        log_avgsrv2sel_m;                     // average survey selectivity

  init_vector       log_srv7_sel_coffs_f(1,n_srv1_sel_ages,phase_selcoff_srv7); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv7_f(-3,5,ph_srv1_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv7_f(-4,4,-4);        //never use this for power fxn so leave phase neg         // age between 50% selection and 95% selection....
  number          a50_srv7_f;                     // age at 50% selection                                                   
  init_number       log_d50_srv7_f(phase_dlogist_srv7);
  number        delta_srv7_f;                     // age between 50% selection and 95% selection....
  init_bounded_number log_gamma_srv7_f(-10,-0.0000000000000001,-4);
  vector        log_srv7_sel_f(1,nages);                       // vector of fishery selectivy log parameters including those not estimated
  vector        srv7_sel_f(1,nages);                           // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgsrv7sel_f;                              // average fishery selectivity
  number        gamma_srv7_f;
  number        d50_srv7_f;

  init_vector       log_srv7_sel_coffs_m(1,n_srv1_sel_ages,phase_selcoff_srv7); // vector of fishery selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv7_m(-3,5,ph_srv1_sel);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv7_m(-4,4,-4);         //never use this for the power fxn so leave phase neg         // age between 50% selection and 95% selection....
  number          a50_srv7_m;                     // age at 50% selection                                                   
  init_number       log_d50_srv7_m(phase_dlogist_srv7);
  number        delta_srv7_m;                     // age between 50% selection and 95% selection....
  init_bounded_number log_gamma_srv7_m(-10,-0.0000000000000001,-4);
  vector        log_srv7_sel_m(1,nages);                       // vector of fishery selectivy log parameters including those not estimated
  vector        srv7_sel_m(1,nages);                           // vectory of fishery selectivty parameters on arithmetic scale
  number        log_avgsrv7sel_m;                              // average fishery selectivity
  number        gamma_srv7_m;
  number        d50_srv7_m;

  init_vector       log_srv10_sel_coffs_f(1,n_srv1_sel_ages,phase_selcoff_srv10);   // vector of survey selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv10_f(-1,4,ph_LL_block2);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv10_f(-5,4,ph_srv_sel_delt_alt);                 // age between 50% selection and 95% selection....
  number          a50_srv10_f;                     // age at 50% selection                                                   
  init_number       d50_srv10_f(phase_dlogist_srv10);
  number        delta_srv10_f;                     // age between 50% selection and 95% selection....
  init_number       gamma_srv10_f(phase_dlogist_srv10);
  vector        log_srv10_sel_f(1,nages);              // vector of survey selectivy log parameters including those not estimated
  vector        srv10_sel_f(1,nages);                  // vectory of survey selectivty parameters on arithmetic scale
  number        log_avgsrv10sel_f;                     // average survey selectivity

  init_vector       log_srv10_sel_coffs_m(1,n_srv1_sel_ages,phase_selcoff_srv10);   // vector of survey selectivy log parameters up until they are constant
  init_bounded_number   log_a50_srv10_m(-1,4,ph_LL_block2);                 // age at 50% selection                                                   
  init_bounded_number   log_delta_srv10_m(-5,4,ph_srv_sel_delt_alt);                 // age between 50% selection and 95% selection....
  number          a50_srv10_m;                     // age at 50% selection                                                   
  init_number       d50_srv10_m(phase_dlogist_srv10);
  number        delta_srv10_m;                     // age between 50% selection and 95% selection....
  init_number       gamma_srv10_m(phase_dlogist_srv10);
  vector        log_srv10_sel_m(1,nages);              // vector of survey selectivy log parameters including those not estimated
  vector        srv10_sel_m(1,nages);                  // vectory of survey selectivty parameters on arithmetic scale
  number        log_avgsrv10sel_m;                     // average survey selectivity

 // ########################################################################################################
 // ########################################################################################################
 // ########################################################################################################
 // ########################################################################################################
 // ########################################################################################################

// Mortality and Survivorship Matrices

  number        F50;                      // Standard deviation report for F50
  sdreport_number       F40;                      // " " " F40
  number        F35;                      // " " " F35
  number        mF50;
  number        mF40;
  number        mF35;
  number        SB0                         // Spawning biomass at no fishing
  number        SBF50                         // " " at F50
  number        SBF40                         // " " at F40
  number        SBF35                         // " " at F35
  number        sprpen                          // Likelihood penalty to make ADMB estimate spr rates
  matrix        Nspr(1,4,1,nages)                   // Matrix of number of spawners at age at each fishing mortality level
  
  matrix                Z_f(styr,endyr,1,nages);                    // Total mortality by year and age
  matrix                Z_m(styr,endyr,1,nages);                    // Total mortality by year and age
  matrix                S_f(styr,endyr,1,nages);                    // Survivorship by year and age
  matrix                S_m(styr,endyr,1,nages);                    // Survivorship by year and age
  matrix                S_f_mid(styr,endyr,1,nages);                    // Survivorship by year and age
  matrix                S_m_mid(styr,endyr,1,nages);                    // Survivorship by year and age
  
  sdreport_number         q_srv1;                                  // " " for Survey catchability
  sdreport_number         q_srv2;
  number                  q_srv3;  //not used now that RPWs are taken out
  number                  q_srv4;
  number                  q_srv5;
  sdreport_number         q_srv6;
  sdreport_number         q_srv7;
  sdreport_number         q_srv8;
  number                  q_srv9;
  sdreport_number         q_LL_fish_recent;
  number                  q_LL_srvy_recent;
  sdreport_number         M_est;

 // Stock-Recruit
  vector                sam_rec(styr_rec,endyr)              // As estimated by assessment model
  vector                srm_rec(styr_rec,endyr)              // As estimated by stock-recruitment curve
  number                sigrsq                 // Recruitment variance parameter
  number                alpha;                  // alpha parameter for B-H
  number                beta;                 // beta parameter for B-H
  number                Bzero;                  // Virgin spawner biomass
  number                Rzero;                  // Virgin recruitment
  number                phizero;               // SPR
  number                log_Rztemp;       // temporary logRzero
  sdreport_vector       pred_rec(styr,endyr);         // " " for predicted recruitments
  number                avg_rec;                                 // " " for Average recruitment
  number                b_yr_end;                       //bias ramp desc limb starting year
  vector                b(styr-nages+2,endyr_rec_est);  // bias ramp adjustment
  vector                rec_like_bias_adj(styr-nages+2,endyr_rec_est);  // bias ramp adjustment
  number                sigr;

// Numbers at age and biomass
  matrix                natage_m(styr,endyr,1,nages);         // Matrix of numbers at age from start year to end year
  matrix                natage_f(styr,endyr,1,nages);         // Matrix of numbers at age from start year to end year
  matrix                num_len_m(styr,endyr,1,nlenbins);         // Matrix of numbers at length from start year to end year
  matrix                num_len_f(styr,endyr,1,nlenbins);         // Matrix of numbers at length from start year to end year

  vector                tot_biom(styr,endyr);
  number                spbiom_trend;                            // " " of Trend in biomass over last 2 years (B(t)/B(t-1); t=endyr)
  number                Depletion;                               // " " for Depletion
  vector                spawn_biom(styr,endyr);         // " " for spawning biomass vector
  vector                Sp_Biom(styr_sp,endyr)
  sdreport_vector       ssbsd(styr,endyr);
  
// Catch at age
  matrix                catage_fish1_f(styr,endyr,1,nages)            // Matrix of predicted catch at age from start year to endyear
  matrix                catage_fish1_m(styr,endyr,1,nages)            // Matrix of predicted catch at age from start year to endyear
  matrix                catage_fish3_f(styr,endyr,1,nages)            // Matrix of predicted catch at age from start year to endyear
  matrix                catage_fish3_m(styr,endyr,1,nages)            // Matrix of predicted catch at age from start year to endyear
  vector                pred_catch_fish1(styr,endyr)        // Vector of predicted catches
  vector                pred_catch_fish3(styr,endyr)
  vector                pred_srv1(styr,endyr);                  // Predicted survey
  vector                pred_srv2(styr,endyr);                  // Predicted survey
  vector                pred_srv3(styr,endyr);                  // Predicted survey
  vector                pred_srv4(styr,endyr);                  // Predicted survey
  vector                pred_srv5(styr,endyr);                  // Predicted survey
  vector                pred_srv6(styr,endyr);                  // Predicted survey
  vector                pred_srv7(styr,endyr);                  // Predicted survey
  vector                pred_srv8(styr,endyr);                  // Predpreedicted survey

  matrix                eac_fish1(1,nyrs_fish1_age,1,nages)        // Expected proportion at age in fish
  matrix                eac_srv1(1,nyrs_srv1_age,1,nages)        // Expected proportion at age in survey
  matrix                eac_srv2(1,nyrs_srv2_age,1,nages)        // Expected proportion at age in survey
  matrix                esc_fish1_m(1,nyrs_fish1_size,1,nlenbins)    // Expected proportion at size in fishery
  matrix                esc_fish1_f(1,nyrs_fish1_size,1,nlenbins)    // Expected proportion at size in fishery
  matrix                esc_fish2(1,nyrs_fish2_size,1,nlenbins)    // Expected proportion at size in fishery
  matrix                esc_fish3_m(1,nyrs_fish3_size,1,nlenbins)    // Expected proportion at size in fishery
  matrix                esc_fish3_f(1,nyrs_fish3_size,1,nlenbins)    // Expected proportion at size in fishery
  matrix                esc_fish4(1,nyrs_fish4_size,1,nlenbins)    // Expected proportion at size in fishery
  matrix                esc_srv1_m(1,nyrs_srv1_size,1,nlenbins)    // Expected proportion at size in survey
  matrix                esc_srv1_f(1,nyrs_srv1_size,1,nlenbins)    // Expected proportion at size in survey
  matrix                esc_srv2_m(1,nyrs_srv2_size,1,nlenbins)    // Expected proportion at size in survey
  matrix                esc_srv2_f(1,nyrs_srv2_size,1,nlenbins)    // Expected proportion at size in survey
  matrix                esc_srv7_m(1,nyrs_srv7_size,1,nlenbins)    // Expected proportion at size in survey
  matrix                esc_srv7_f(1,nyrs_srv7_size,1,nlenbins)    // Expected proportion at size in survey
  matrix                eac_srv7(1,nyrs_srv7_age,1,nages)    // Expected proportion at size in survey

// Likelihoods and penalty functions
  vector        surv_like(1,8);   // Likelihood values for survey biomasses, allowance for up to 3 surveys
  vector        age_like(1,16);     // Likelihood values for age and size compositions allowance for up 6 comps
  vector        sel_like(1,12);     // LIkelihood values for selectivities with alowance for up to 6 selectivities
  number        rec_like;         // Likelihood value for recruitments
  number        ssqcatch;         // Likelihood value for catch estimation
  number        F_mort_regularity;  // Penalty value for fishing mortality regularity
  number        M_mort_regularity;  // Penalty value for fishing mortality regularity
  number        avg_sel_penalty;    // Penalty value for selectivity regularity penalty

 // Priors
  vector        priors(1,13);       // Prior penalty values for sigr,q,natural mortality,steepness
  
// Define an objective function
  number        Like;             // Likelihood for data fits
  objective_function_value obj_fun;       // Total likelihood for objective function value


  vector        xdum2(styr,endyr);                // Dummy variable for use in pop-report.cxx
  vector        pred_catch(styr,endyr);
  number        fratio;
//  number        B40;
  number        ABC3;
  
///////////////////////////////////////////
/// Population projection Hanselsiginelli
//////////////////////////////////////////

  matrix        N_proj_f(endyr+1,endyr+15,1,nages);
  matrix        N_proj_m(endyr+1,endyr+15,1,nages);
  number        FABC_proj;
  vector        FABC_tot_proj_f(1,nages);
  vector        FABC_tot_proj_m(1,nages);
  number        FOFL_proj;
  vector        FOFL_tot_proj_f(1,nages);
  vector        FOFL_tot_proj_m(1,nages);
  sdreport_number ABC;
  sdreport_number B40;
  number        OFL;
  vector        Z_proj_f(1,nages);
  vector        Z_proj_m(1,nages);
  vector        ZOFL_proj_f(1,nages);
  vector        ZOFL_proj_m(1,nages);

  vector        S_proj_f(1,nages);
  vector        S_proj_m(1,nages);
  matrix        catage_proj_f(endyr+1,endyr+15,1,nages);
  matrix        catage_proj_m(endyr+1,endyr+15,1,nages);
  matrix        catage_proj_OFL_f(endyr+1,endyr+15,1,nages);
  matrix        catage_proj_OFL_m(endyr+1,endyr+15,1,nages);
  vector        pred_catch_proj_OFL_f(endyr+1,endyr+15);
  vector        pred_catch_proj_OFL_m(endyr+1,endyr+15);
  sdreport_vector         spawn_biom_proj(endyr+1,endyr+15);
  sdreport_vector         tot_biom_proj(endyr+1,endyr+15);
  sdreport_vector         pred_catch_proj(endyr+1,endyr+15);
  sdreport_vector         pred_catch_proj_OFL(endyr+1,endyr+15);
  number        stdev_rec;
  number        FOFL;
  number        FABC;
  number        FOFL2;
  number        FABC2; 


  3darray size_age_f(styr,endyr,1,nages,1,nlenbins)
  3darray size_age_m(styr,endyr,1,nages,1,nlenbins)
  matrix weight_f(styr,endyr,1,nages)
  matrix weight_m(styr,endyr,1,nages)
  matrix maturity(styr,endyr,1,nages)
  matrix weight_maturity_prod_f(styr,endyr,1,nages);                  // Weight of mature fish vector at age


 PROCEDURE_SECTION
  l=l+1; // Initiate counter for random seeds in projection
  
  switch (SrType) {
     case 3:
     if (rec_like_type==1) log_Rztemp=(log_mean_rec);
      default:
      log_Rztemp=log_Rzero; }
     Get_Biologicals();
     Get_Selectivity();                     // Call function to get selectivities
       Get_Mortality_Rates();                 // Call function to get fishing and natural mortality
   //  Get_Bzero();                           // OjO
       Get_Numbers_At_Age();                    // Call function to get numbers at age per year
       Get_Catch_at_Age();                      // Call function to get catch at age per year
       Get_Predicted_Values();                  // Get predicted values for catch, survbio, age and size comps
       Calc_priors();                         // Solve for priors
 //  if (last_phase())
 //  {
        
       Get_Dependent_Vars();                  // Solve for dependent variables like total bio, recruitment etc.
       ssbsd = spawn_biom;
       compute_spr_rates();                 // Compute f40 etc.
       Get_Population_Projection();
 //  }
 switch (SrType) {
     case 3:
      if (rec_like_type==1)log_Rztemp=(log_mean_rec);
      default:
      log_Rztemp=log_Rzero; }
  Evaluate_Objective_Function();          // Minimize objective function value
  
    if (mceval_phase())                     // For outputting MCMC simulations in text format 
  
    {
     evalout<<log_mean_rec<<" "<<sigr<<" "<<q_srv1<<" "<<q_srv2<<" "<<q_srv3<<" "<<q_srv4<<" "<<q_srv5<<" "<<q_srv6<<" "<<q_srv7<<" "<<q_srv8<<" "<<F40<<" "<<natmort<<" "<<obj_fun<<" "<<tot_biom<<" "<<spawn_biom<<" "<<pred_rec<<" "<<B40<<" "<<B40<<" "<<spawn_biom_proj<<" "<<pred_catch_proj<<" "<<endl;
     }

FUNCTION dvar_vector SRecruit(const dvar_vector& Stmp)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector RecTmp(Stmp.indexmin(),Stmp.indexmax());
      // dvariable R_alpha;
      // dvariable R_beta;
  switch (SrType)
  {
    case 1:
      RecTmp = elem_prod((Stmp / phizero) , mfexp( alpha * ( 1. - Stmp / Bzero ))) ; //Ricker form from Dorn
      // R_alpha = exp(alpha)/phizero;
      // R_alpha  = 0.036494;
      // R_beta  = alpha/Bzero;
      // R_beta  = 0.0037207;
      // RecTmp = elem_prod( R_alpha*Stmp ,mfexp(-R_beta*Stmp)) ; // Ricker model

      break;
    case 2:
      // cout<<"BH a b "<<alpha<<" "<<beta<<endl;
      RecTmp = elem_prod(Stmp , 1. / ( alpha + beta * Stmp));        //Beverton-Holt form
      break;
    case 3:
       RecTmp = mfexp(log_mean_rec);                    //Avg recruitment
      if(rec_like_type==1) log_Rztemp=log_mean_rec;
      break;
    case 4:
      RecTmp = elem_prod(Stmp , mfexp( alpha  - Stmp * beta)) ; //Old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION dvariable SRecruit(const double& Stmp)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
  switch (SrType)
  {
    case 1:
      RecTmp = (Stmp / phizero) * mfexp( alpha * ( 1. - Stmp / Bzero )) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = Stmp / ( alpha + beta * Stmp);        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(log_mean_rec);                    //Avg recruitment
      if(rec_like_type==1) log_Rztemp=log_mean_rec;
      break;
    case 4:
      RecTmp = Stmp * mfexp( alpha  - Stmp * beta) ; //old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION dvariable SRecruit(CONST dvariable& Stmp)
  RETURN_ARRAYS_INCREMENT();
  dvariable RecTmp;
      // dvariable R_alpha;
      // dvariable R_beta;
  switch (SrType)
  {
    case 1:
      RecTmp = (Stmp / phizero) * mfexp( alpha * ( 1. - Stmp / Bzero )) ; //Ricker form from Dorn
      break;
    case 2:
      RecTmp = Stmp / ( alpha + beta * Stmp);        //Beverton-Holt form
      break;
    case 3:
      RecTmp = mfexp(log_mean_rec );                    //Avg recruitment
      break;
    case 4:
      RecTmp = Stmp * mfexp( alpha  - Stmp * beta) ; //old Ricker form
      break;
  }
  RETURN_ARRAYS_DECREMENT();
  return RecTmp;

FUNCTION Get_Biologicals

  sigr=mfexp(log_sigr);

 //allow for multiple growth blocks
 for (i=styr;i<=endyr;i++)
  {
   if(growth_blocks==1)
    {
      size_age_f(i)=sizeage1_f;
      size_age_m(i)=sizeage1_m;   
    }
   if(growth_blocks==2)
    {
     if(i<growth_cutoffs(2))
      {
       size_age_f(i)=sizeage1_f;
       size_age_m(i)=sizeage1_m;
      }
     if(i>=growth_cutoffs(2))
      {
       size_age_f(i)=sizeage2_f;
       size_age_m(i)=sizeage2_m;
      }      
    }
   if(growth_blocks==3)
    {
     if(i<growth_cutoffs(2))
      {
       size_age_f(i)=sizeage1_f;
       size_age_m(i)=sizeage1_m;
      }
     if(i>=growth_cutoffs(2) && i<growth_cutoffs(3) )
      {
       size_age_f(i)=sizeage2_f;
       size_age_m(i)=sizeage2_m;
      }      
     if(i>=growth_cutoffs(3))
      {
       size_age_f(i)=sizeage3_f;
       size_age_m(i)=sizeage3_m;
      }
    }
  }

 //allow for multiple weight blocks
 for (i=styr;i<=endyr;i++)
  {
   if(weight_blocks==1)
    {
     weight_f(i)=wt_f1;
     weight_m(i)=wt_m1;   
    }
   if(weight_blocks==2)
    {
     if(i<weight_cutoffs(2))
      {
       weight_f(i)=wt_f1;
       weight_m(i)=wt_m1; 
      }
     if(i>=growth_cutoffs(2))
      {
       weight_f(i)=wt_f2;
       weight_m(i)=wt_m2; 
      }      
    }
  }

 //allow for multiple maturity blocks
 for (i=styr;i<=endyr;i++)
  {
   if(maturity_blocks==1)
    {
     maturity(i)=p_mature1;
    }
   if(maturity_blocks==2)
    {
     if(i<maturity_cutoffs(2))
      {
       maturity(i)=p_mature1;
      }
     if(i>=maturity_cutoffs(2))
      {
       maturity(i)=p_mature2;
      }      
    }
   if(maturity_blocks==3)
    {
     if(i<maturity_cutoffs(2))
      {
       maturity(i)=p_mature1;
      }
     if(i>=maturity_cutoffs(2) && i<maturity_cutoffs(3))
      {
       maturity(i)=p_mature2;
      }
     if(i>=maturity_cutoffs(3))
      {
       maturity(i)=p_mature3;
      } 
    }
  }

     weight_maturity_prod_f = elem_prod(weight_f,maturity);

FUNCTION Get_Bzero
  Bzero.initialize();
  Rzero    =  mfexp(log_Rztemp); 
  sigrsq   = sigr*sigr;

  dvariable survtmp = exp(-natmort);
  dvariable spawn_adj=pow(survtmp,spawn_fract) ;

  dvar_matrix natagetmp(styr_rec,styr,1,nages);
  natagetmp(styr_rec,1) = Rzero;
  for (j=2; j<=nages; j++)
    natagetmp(styr_rec,j) = natagetmp(styr_rec,j-1) * survtmp;

  natagetmp(styr_rec,nages) /= (1.-survtmp); 

  Bzero = weight_maturity_prod_f(styr) * spawn_adj *natagetmp(styr_rec) ;
  phizero = Bzero/Rzero;

  switch (SrType)
  {
    case 1:
      alpha = log(-4.*(steepness)/(steepness-1.));
      break;
    case 2:
      alpha  =  Bzero * (1. - (steepness - 0.2) / (0.8*steepness) ) / Rzero;
      beta   = (5. * steepness - 1.) / (4. * steepness * Rzero);
      break;
    case 4:
    //R = S * EXP(alpha - beta * S))
      beta  = log(5.*(steepness))/(0.8*Bzero) ;
      alpha = log((Rzero)/Bzero)+beta*Bzero;
      break;
  }

  Sp_Biom.initialize();
  Sp_Biom(styr_sp,styr_rec-1) = Bzero;
  for (i=styr_rec;i<styr;i++)
  {
    natagetmp(i,1)          = mfexp(log_rec_dev(i) + log_Rztemp);
    Sp_Biom(i) = natagetmp(i)*spawn_adj * weight_maturity_prod_f(styr); 
    natagetmp(i+1)(2,nages) = ++(natagetmp(i)(1,nages-1)*mfexp(-natmort ));
    natagetmp(i+1,nages)   += natagetmp(i,nages)*mfexp(-natmort);
  }
  natagetmp(styr,1)   = mfexp(log_rec_dev(styr) + log_Rztemp);
  sam_rec(styr_rec,styr) = column(natagetmp,1);
  natage_f(styr)  = natagetmp(styr)/2; // OjO
  natage_m(styr)  = natagetmp(styr)/2; // OjO
  Sp_Biom(styr) = natagetmp(styr)*spawn_adj * weight_maturity_prod_f(styr); 

 // OjO to here...
  
FUNCTION Get_Selectivity

//   Selectivity does not change for ages greater than n_fish_sel_ages

// ###################################################################################################
// #####################################################################################################
// ####  FISHERY SELECTIVITY ####################
// #####################################################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ###################################################################################################
// #####################################################################################################
// ####  BASE MODEL FISHERY SELECTIVITY PARAMETERS ####################
// #####################################################################################################
// ###################################################################################################

// LOG-TRANSFORM FISHERY SELECTIVITY ESTIMATED PARAMETERS
// ###################################################################################################

 // All fishery fleets have their own sex specific a50
  a50_fish1_f=mfexp(log_a50_fish1_f);
  a50_fish1_m=mfexp(log_a50_fish1_m);
  a50_fish2=mfexp(log_a50_fish1_f);  //link because keeps bounding 
  a50_fish3_f=mfexp(log_a50_fish3_f); // FOR GAMMA FXN
  a50_fish3_m=mfexp(log_a50_fish3_m); // FOR GAMMA FXN
  a50_fish4_f=mfexp(log_a50_fish4_f);
  a50_fish4_m=mfexp(log_a50_fish4_m);
  a50_fish5_f=mfexp(log_a50_fish5_f);
  a50_fish5_m=mfexp(log_a50_fish5_m);
  
  delta_fish1_f=mfexp(log_delta_fish1_f);    // linked because only 5 years length comp data (1990-1995) to est historic sel
  delta_fish1_m=mfexp(log_delta_fish1_f);    // linked because only 5 years length comp data (1990-1995) to est historic sel
  delta_fish2=mfexp(log_delta_fish1_f);     // linked because only 5 years length comp data (1990-1995) to est historic sel
  delta_fish3_f=mfexp(log_delta_fish3_f);
  delta_fish3_m=mfexp(log_delta_fish3_f);    // linked bec parameters not well est prob due to limited sex-specific data, does not appear much diff in delta betw sexes so fixed by sex
  delta_fish4_f=mfexp(log_delta_fish4_f);   
  delta_fish4_m=mfexp(log_delta_fish4_m);   
  delta_fish5_f=mfexp(log_delta_fish5_f);   
  delta_fish5_m=mfexp(log_delta_fish5_m);   

// ###################################################################################################

// CALCULATE FISHERY SELECTIVITY BASED ON EST PARAMETERS
// ###################################################################################################
  for (j=1;j<=nages;j++)
   {
   
    if(fsh1_sel_opt==2)
     {
      fish1_sel_f(j)=1/ (1+mfexp(-delta_fish1_f*(j-a50_fish1_f))); 
      fish1_sel_m(j)=1/ (1+mfexp(-delta_fish1_m*(j-a50_fish1_m)));
     }

    if(fsh2_sel_opt==2) fish2_sel(j)=1/ (1+mfexp(-delta_fish2*(double(j)-a50_fish2))); 

    if(fsh3_sel_opt==3) // Punt et. al 1996 gamma parameterization
     { 
      fish3_sel_f(j)=(pow(j/a50_fish3_f,a50_fish3_f/(0.5*(sqrt(square(a50_fish3_f)+4*square(delta_fish3_f))-a50_fish3_f)))*mfexp((a50_fish3_f-j)/(0.5*(sqrt(square(a50_fish3_f)+4*square(delta_fish3_f))-a50_fish3_f))));
      fish3_sel_m(j)=(pow(j/a50_fish3_m,a50_fish3_m/(0.5*(sqrt(square(a50_fish3_m)+4*square(delta_fish3_m))-a50_fish3_m)))*mfexp((a50_fish3_m-j)/(0.5*(sqrt(square(a50_fish3_m)+4*square(delta_fish3_m))-a50_fish3_m))));
     }
 
    if(fsh4_sel_opt==2)
     {
      fish4_sel_f(j)=1/ (1+mfexp(-delta_fish4_f*(j-a50_fish4_f)));  
      fish4_sel_m(j)=1/ (1+mfexp(-delta_fish4_m*(j-a50_fish4_m)));
     }

    if(fsh5_sel_opt==2)
     {
      fish5_sel_f(j)=1/ (1+mfexp(-delta_fish5_f*(j-a50_fish5_f)));  
      fish5_sel_m(j)=1/ (1+mfexp(-delta_fish5_m*(j-a50_fish5_m)));
     }

   } // END j SUBSCRIPT (NAGES) FOR FISHERY SELECTIVITY
  
   fish2_sel=fish2_sel/max(fish2_sel);       //standardize logistic FXN to max of 1
   fish3_sel_f=fish3_sel_f/max(fish3_sel_f); //standardize gamma FXN to max of 1
   fish3_sel_m=fish3_sel_m/max(fish3_sel_m); //standardize gamma FXN to max of 1
 
// #####################################################################################################
// #####################################################################################################
// #####################################################################################################
// #####################################################################################################



// ###################################################################################################
// #####################################################################################################
// ####  SURVEY SELECTIVITY ####################
// #####################################################################################################
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ###################################################################################################
// #####################################################################################################
// ####  BASE MODEL SURVEY SELECTIVITY PARAMETERS ####################
// #####################################################################################################
// ###################################################################################################

// LOG-TRANSFORM SURVEY SELECTIVITY ESTIMATED PARAMETERS
// ###################################################################################################
     // all surveys have sex specific a50s
      a50_srv1_f=mfexp(log_a50_srv1_f);
      a50_srv1_m=mfexp(log_a50_srv1_m);
      a50_srv2_f=mfexp(log_a50_srv2_f);
      a50_srv2_m=mfexp(log_a50_srv2_m);
      a50_srv7_f=mfexp(log_a50_srv7_f);
      a50_srv7_m=mfexp(log_a50_srv7_m);
      a50_srv10_f=mfexp(log_a50_srv10_f);
      a50_srv10_m=mfexp(log_a50_srv10_m);

      delta_srv1_f=mfexp(log_delta_srv1_f);    //all deltas linked due to trouble estimating and unrealistic resultant sel curves, appears to be more variation in delta by sex than year based on exploratory runs
      delta_srv1_m=mfexp(log_delta_srv1_m);
      delta_srv2_f=mfexp(log_delta_srv1_f);
      delta_srv2_m=mfexp(log_delta_srv1_m);    
      delta_srv10_f=mfexp(log_delta_srv1_f);
      delta_srv10_m=mfexp(log_delta_srv1_m);
      
// ###################################################################################################

// CALCULATE SURVEY SELECTIVITY BASED ON EST PARAMETERS
// ###################################################################################################
 for (j=1;j<=nages;j++)
  {
  
   if(srv1_sel_opt==2)
    { 
     srv1_sel_f(j)=1/ (1+mfexp(-delta_srv1_f*(j-a50_srv1_f)));
     srv1_sel_m(j)=1/ (1+mfexp(-delta_srv1_m*(j-a50_srv1_m)));
    }

   if(srv2_sel_opt==2)
    { 
     srv2_sel_f(j)=1/ (1+mfexp(-delta_srv2_f*(j-a50_srv2_f)));
     srv2_sel_m(j)=1/ (1+mfexp(-delta_srv2_m*(j-a50_srv2_m)));
    }

   if(srv7_sel_opt==3) // POWER FUNCTION FOR TRAWL SURVEY
    { 
      srv7_sel_f(j)=(1/pow(j,a50_srv7_f));
      srv7_sel_m(j)=(1/pow(j,a50_srv7_m));
    }

   if(srv10_sel_opt==2)
    { 
     srv10_sel_f(j)=1/ (1+mfexp(-delta_srv10_f*(j-a50_srv10_f)));
     srv10_sel_m(j)=1/ (1+mfexp(-delta_srv10_m*(j-a50_srv10_m)));
    }
    
  }  // END j SUBSCRIPT (NAGES) FOR SURVEY SELECTIVITY CALCS

    srv7_sel_f /= max(srv7_sel_f);   // standardize power function to max of 1.0
    srv7_sel_m /= max(srv7_sel_m);   // standardize power function to max of 1.0

// #####################################################################################################
// #####################################################################################################
// #####################################################################################################
// #####################################################################################################






// ###################################################################################################
// #####################################################################################################
// ####  OTHER FISHERY SELECTIVITY OPTIONS ####################
// #####################################################################################################
// #####################################################################################################

  gamma_fish3_m=mfexp(log_gamma_fish3_m); // used for ALT GAMMA FXN, SEL OPT = 4
  gamma_fish3_f=mfexp(log_gamma_fish3_f); // used for ALT GAMMA FXN, SEL OPT = 4
  
  if (fsh1_sel_opt==1)
   {
    for (j=1;j<=n_fish_sel_ages;j++)
    {
      log_fish1_sel_f(j) = log_fish1_sel_coffs_f(j);
      log_fish1_sel_m(j) = log_fish1_sel_coffs_m(j);
      log_fish1_sel_f(j) = log_fish1_sel_f(j-1);
      log_fish1_sel_m(j) = log_fish1_sel_m(j-1); 
      log_avgfish1sel_f = log(mean(mfexp(log_fish1_sel_coffs_f)));
      log_fish1_sel_f  -= log(mean(mfexp(log_fish1_sel_f)));
      log_avgfish1sel_m = log(mean(mfexp(log_fish1_sel_coffs_m)));
      log_fish1_sel_m  -= log(mean(mfexp(log_fish1_sel_m)));
      fish1_sel_f       = mfexp(log_fish1_sel_f)/mfexp(max(log_fish1_sel_f));  // keeping maximum fish selectivity at 1
      fish1_sel_m       = mfexp(log_fish1_sel_m)/mfexp(max(log_fish1_sel_m));  // keeping maximum fish selectivity at 1
    }
   }
  
  if (fsh2_sel_opt==1)
   {
    for (j=1;j<=n_fish_sel_ages;j++) log_fish2_sel(j) = log_fish2_sel_coffs(j);
      log_avgfish2sel = log(mean(mfexp(log_fish2_sel_coffs)));
      log_fish2_sel  -= log(mean(mfexp(log_fish2_sel)));
      fish2_sel       = mfexp(log_fish2_sel)/mfexp(max(log_fish2_sel));
   } // keeping maximum fish selectivity at 1

   if (fsh3_sel_opt==1) 
    {
     for (j=1;j<=n_fish_sel_ages;j++)
      {
       log_fish3_sel_f(j) = log_fish3_sel_coffs_f(j);
       log_fish3_sel_m(j) = log_fish3_sel_coffs_m(j);
       log_fish3_sel_f(j) = log_fish3_sel_f(j-1);
       log_fish3_sel_m(j) = log_fish3_sel_m(j-1); 
       log_avgfish3sel_f = log(mean(mfexp(log_fish3_sel_coffs_f)));
       log_fish3_sel_f  -= log(mean(mfexp(log_fish3_sel_f)));
       log_avgfish3sel_m = log(mean(mfexp(log_fish3_sel_coffs_m)));
       log_fish3_sel_m  -= log(mean(mfexp(log_fish3_sel_m)));
       fish3_sel_f       = mfexp(log_fish3_sel_f)/mfexp(max(log_fish3_sel_f));  // keeping maximum fish selectivity at 1
       fish3_sel_m       = mfexp(log_fish3_sel_m)/mfexp(max(log_fish3_sel_m));  // keeping maximum fish selectivity at 1
      }
    }
    
  if (fsh4_sel_opt==1)
   {
    for (j=1;j<=n_fish_sel_ages;j++)
     {
      log_fish4_sel_f(j) = log_fish4_sel_coffs_f(j);
      log_fish4_sel_m(j) = log_fish4_sel_coffs_m(j);
      log_fish4_sel_f(j) = log_fish4_sel_f(j-1);
      log_fish4_sel_m(j) = log_fish4_sel_m(j-1);
     }
      log_avgfish4sel_f = log(mean(mfexp(log_fish4_sel_coffs_f)));
      log_fish4_sel_f  -= log(mean(mfexp(log_fish4_sel_f)));
      log_avgfish4sel_m = log(mean(mfexp(log_fish4_sel_coffs_m)));
      log_fish4_sel_m  -= log(mean(mfexp(log_fish4_sel_m)));
      fish4_sel_f       = mfexp(log_fish4_sel_f)/mfexp(max(log_fish4_sel_f));  // keeping maximum fish selectivity at 1
      fish4_sel_m       = mfexp(log_fish4_sel_m)/mfexp(max(log_fish4_sel_m));
   }  // keeping maximum fish selectivity at 1
   
// #####################################################################################################
// #####################################################################################################

  for (j=1;j<=nages;j++)
   {
   
    if(fsh1_sel_opt==3)
     { 
      fish1_sel_f(j) = (1/(1+mfexp(-delta_fish1_f*(double(j)-a50_fish1_f))))*(1-(1/(1+mfexp(-gamma_fish1_f*(j-d50_fish1_f)))));
      fish1_sel_m(j) = (1/(1+mfexp(-delta_fish1_m*(double(j)-a50_fish1_m))))*(1-(1/(1+mfexp(-gamma_fish1_m*(j-d50_fish1_m)))));
     }
     
    if(fsh1_sel_opt==4)
     { 
      fish1_sel_f(j) = (1/(1-gamma_fish1_f))*pow(((1-gamma_fish1_f)/gamma_fish1_f),gamma_fish1_f)*(mfexp(delta_fish1_f*gamma_fish1_f*(a50_fish1_f-j))/(1+mfexp(delta_fish1_f*(a50_fish1_f-j))));
      fish1_sel_m(j) = (1/(1-gamma_fish1_m))*pow(((1-gamma_fish1_m)/gamma_fish1_m),gamma_fish1_m)*(mfexp(delta_fish1_m*gamma_fish1_m*(a50_fish1_m-j))/(1+mfexp(delta_fish1_m*(a50_fish1_m-j))));
     }
     
    if(fsh1_sel_opt==5)
     {
      fish1_sel_f(j)=1/ (1+mfexp(-delta_fish1_f*(j-a50_fish1_f))); 
      fish1_sel_m(j)=1/ (1+mfexp(-delta_fish1_m*(j-a50_fish1_m)));
     }

    if(fsh2_sel_opt==3)  fish2_sel(j)=(pow(j/a50_fish2,a50_fish2/(0.5*(sqrt(square(a50_fish2)+4*square(delta_fish2))-a50_fish2)))*mfexp((a50_fish2-j)/(0.5*(sqrt(square(a50_fish2)+4*square(delta_fish2))-a50_fish2))));

    if(fsh2_sel_opt==5) fish2_sel(j)=1/ (1+mfexp(-delta_fish2*(double(j)-a50_fish2)));
    
    if(fsh3_sel_opt==2)
     {
      fish3_sel_f(j)=1/ (1+mfexp(-delta_fish3_f*(j-a50_fish3_f))); 
      fish3_sel_m(j)=1/ (1+mfexp(-delta_fish3_m*(j-a50_fish3_m)));
     }
     
    if(fsh3_sel_opt==4)
     { 
      fish3_sel_f(j) = (1/(1-gamma_fish3_f))*pow(((1-gamma_fish3_f)/gamma_fish3_f),gamma_fish3_f)*(mfexp(delta_fish3_f*gamma_fish3_f*(a50_fish3_f-j))/(1+mfexp(delta_fish3_f*(a50_fish3_f-j))));
      fish3_sel_m(j) = (1/(1-gamma_fish3_m))*pow(((1-gamma_fish3_m)/gamma_fish3_m),gamma_fish3_m)*(mfexp(delta_fish3_m*gamma_fish3_m*(a50_fish3_m-j))/(1+mfexp(delta_fish3_m*(a50_fish3_m-j))));
     }

    if(fsh4_sel_opt==3)
     {
      fish3_sel_f(j)=(pow(j/a50_fish3_f,a50_fish3_f/(0.5*(sqrt(square(a50_fish3_f)+4*square(delta_fish3_f))-a50_fish3_f)))*mfexp((a50_fish3_f-j)/(0.5*(sqrt(square(a50_fish3_f)+4*square(delta_fish3_f))-a50_fish3_f))));
      fish3_sel_m(j)=(pow(j/a50_fish3_m,a50_fish3_m/(0.5*(sqrt(square(a50_fish3_m)+4*square(delta_fish3_m))-a50_fish3_m)))*mfexp((a50_fish3_m-j)/(0.5*(sqrt(square(a50_fish3_m)+4*square(delta_fish3_m))-a50_fish3_m))));
     }
     
    if(fsh4_sel_opt==4)
     { 
      fish4_sel_f(j) = (1/(1-gamma_fish4_f))*pow(((1-gamma_fish4_f)/gamma_fish4_f),gamma_fish4_f)*(mfexp(delta_fish4_f*gamma_fish4_f*(a50_fish4_f-j))/(1+mfexp(delta_fish4_f*(a50_fish4_f-j))));
      fish4_sel_m(j) = (1/(1-gamma_fish4_m))*pow(((1-gamma_fish4_m)/gamma_fish4_m),gamma_fish4_m)*(mfexp(delta_fish4_m*gamma_fish4_m*(a50_fish4_m-j))/(1+mfexp(delta_fish4_m*(a50_fish4_m-j))));
     }

   } // END j SUBSCRIPT (NAGES) FOR FISHERY SELECTIVITY

   fish2_sel=fish2_sel/max(fish2_sel);
   fish3_sel_f=fish3_sel_f/max(fish3_sel_f);
   fish3_sel_m=fish3_sel_m/max(fish3_sel_m);



// ###################################################################################################
// #####################################################################################################
// #### OTHER SURVEY SELECTIVITY OPTIONS ####################
// #####################################################################################################
// ###################################################################################################

      delta_srv7_f=mfexp(log_delta_srv7_f); // used for GAMMA FXN or LOGISTIC FXN, SEL OPT = 2 or 4 or 5
      delta_srv7_m=mfexp(log_delta_srv7_m); // used for GAMMA FXN or LOGISTIC FXN, SEL OPT = 2 or 4 or 5
      
      gamma_srv7_f=mfexp(log_gamma_srv7_f); // used for GAMMA FXN, SEL OPT = 4 or 5
      gamma_srv7_m=mfexp(log_gamma_srv7_m); // used for GAMMA FXN, SEL OPT = 4 or 5


    if (srv1_sel_opt==1)
     {
      for (j=1;j<=n_srv1_sel_ages;j++)
       {
        log_srv1_sel_f(j) = log_srv1_sel_coffs_f(j);
        log_srv1_sel_m(j) = log_srv1_sel_coffs_m(j); 
        log_srv1_sel_f(j) = log_srv1_sel_f(j-1);
        log_avgsrv1sel_f    = log(mean(mfexp(log_srv1_sel_coffs_f)));
        log_srv1_sel_f     -= log(mean(mfexp(log_srv1_sel_f)));
        srv1_sel_f          = mfexp(log_srv1_sel_f)/mfexp(max(log_srv1_sel_f));  //keeping max survey selectiviy at 1
        log_srv1_sel_m(j) = log_srv1_sel_m(j-1);
       }
        log_avgsrv1sel_m    = log(mean(mfexp(log_srv1_sel_coffs_m)));
        log_srv1_sel_m     -= log(mean(mfexp(log_srv1_sel_m)));
        srv1_sel_m          = mfexp(log_srv1_sel_m)/mfexp(max(log_srv1_sel_m));  //keeping max survey selectiviy at 1
      }
      
    if (srv2_sel_opt==1)
     {
      for (j=1;j<=n_srv1_sel_ages;j++)
       {
        log_srv2_sel_m(j) = log_srv2_sel_coffs_m(j);
        log_srv2_sel_f(j) = log_srv2_sel_coffs_f(j); 
        log_srv2_sel_f(j) = log_srv2_sel_f(j-1);
        log_avgsrv2sel_f    = log(mean(mfexp(log_srv2_sel_coffs_f)));
        log_srv2_sel_f     -= log(mean(mfexp(log_srv2_sel_f)));
        srv2_sel_f          = mfexp(log_srv2_sel_f)/mfexp(max(log_srv2_sel_f));  //keeping max survey selectiviy at 1
        log_srv2_sel_m(j) = log_srv2_sel_m(j-1);
       }
        log_avgsrv2sel_m    = log(mean(mfexp(log_srv2_sel_coffs_m)));
        log_srv2_sel_m     -= log(mean(mfexp(log_srv2_sel_m)));
        srv2_sel_m          = mfexp(log_srv2_sel_m)/mfexp(max(log_srv2_sel_m));  //keeping max survey selectiviy at 1
      }
      
    if (srv7_sel_opt==1)
     {
      for (j=1;j<=n_srv1_sel_ages;j++)
       {
        log_srv7_sel_f(j) = log_srv7_sel_coffs_f(j); 
        log_srv7_sel_m(j) = log_srv7_sel_coffs_m(j); 
        log_srv7_sel_f(j) = log_srv7_sel_f(j-1);
        log_avgsrv7sel_f    = log(mean(mfexp(log_srv7_sel_coffs_f)));
        log_srv7_sel_f     -= log(mean(mfexp(log_srv7_sel_f)));
        srv7_sel_f          = mfexp(log_srv7_sel_f)/mfexp(max(log_srv7_sel_f));  //keeping max survey selectiviy at 1
        log_srv7_sel_m(j) = log_srv7_sel_m(j-1);
       }
        log_avgsrv7sel_m    = log(mean(mfexp(log_srv7_sel_coffs_m)));
        log_srv7_sel_m     -= log(mean(mfexp(log_srv7_sel_m)));
        srv7_sel_m          = mfexp(log_srv7_sel_m)/mfexp(max(log_srv7_sel_m));  //keeping max survey selectiviy at 1
     }
     
// #####################################################################################################
// #####################################################################################################

 for (j=1;j<=nages;j++)
  {
  
   if(srv1_sel_opt==3)
    { 
      srv1_sel_f(j) = (1/(1+mfexp(-delta_srv1_f*(double(j)-a50_srv1_f))))*(1-(1/(1+mfexp(-gamma_srv1_f*(j-d50_srv1_f)))));
      srv1_sel_m(j) = (1/(1+mfexp(-delta_srv1_m*(double(j)-a50_srv1_m))))*(1-(1/(1+mfexp(-gamma_srv1_m*(j-d50_srv1_m)))));
    }

   if(srv1_sel_opt==4)
    { 
      srv1_sel_f(j) = (1/(1-gamma_srv1_f))*pow(((1-gamma_srv1_f)/gamma_srv1_f),gamma_srv1_f)*(mfexp(delta_srv1_f*gamma_srv1_f*(a50_srv1_f-j))/(1+mfexp(delta_srv1_f*(a50_srv1_f-j))));
      srv1_sel_m(j) = (1/(1-gamma_srv1_m))*pow(((1-gamma_srv1_m)/gamma_srv1_m),gamma_srv1_m)*(mfexp(delta_srv1_m*gamma_srv1_m*(a50_srv1_m-j))/(1+mfexp(delta_srv1_m*(a50_srv1_m-j))));
    }

   if(srv2_sel_opt==3)
    { 
      srv2_sel_f(j) = (1/(1+mfexp(-delta_srv2_f*(double(j)-a50_srv2_f))))*(1-(1/(1+mfexp(-gamma_srv2_f*(j-d50_srv2_f)))));
      srv2_sel_m(j) = (1/(1+mfexp(-delta_srv2_m*(double(j)-a50_srv2_m))))*(1-(1/(1+mfexp(-gamma_srv2_m*(j-d50_srv2_m)))));
    }
    
   if(srv2_sel_opt==4)
    { 
      srv2_sel_f(j) = (1/(1-gamma_srv2_f))*pow(((1-gamma_srv2_f)/gamma_srv2_f),gamma_srv2_f)*(mfexp(delta_srv2_f*gamma_srv2_f*(a50_srv2_f-j))/(1+mfexp(delta_srv2_f*(a50_srv2_f-j))));
      srv2_sel_m(j) = (1/(1-gamma_srv2_m))*pow(((1-gamma_srv2_m)/gamma_srv2_m),gamma_srv2_m)*(mfexp(delta_srv2_m*gamma_srv2_m*(a50_srv2_m-j))/(1+mfexp(delta_srv2_m*(a50_srv2_m-j))));
    }
    
   if(srv7_sel_opt==2)
    { 
      srv7_sel_f(j)=1/ (1+mfexp(-delta_srv7_f*(j-a50_srv7_f)));
      srv7_sel_m(j)=1/ (1+mfexp(-delta_srv7_m*(j-a50_srv7_m)));
    }

   if(srv7_sel_opt==4)
    { 
      srv7_sel_f(j) = (1/(1-gamma_srv7_f))*pow(((1-gamma_srv7_f)/gamma_srv7_f),gamma_srv7_f)*(mfexp(delta_srv7_f*gamma_srv7_f*(a50_srv7_f-j))/(1+mfexp(delta_srv7_f*(a50_srv7_f-j))));
      srv7_sel_m(j) = (1/(1-gamma_srv7_m))*pow(((1-gamma_srv7_m)/gamma_srv7_m),gamma_srv7_m)*(mfexp(delta_srv7_m*gamma_srv7_m*(a50_srv7_m-j))/(1+mfexp(delta_srv7_m*(a50_srv7_m-j))));
    }
    
   if(srv7_sel_opt==5)
    { 
      srv7_sel_f(j) = (1/(1-gamma_srv7_f))*pow(((1-gamma_srv7_f)/gamma_srv7_f),gamma_srv7_f)*(mfexp(delta_srv7_f*gamma_srv7_f*(a50_srv7_f-j))/(1+mfexp(delta_srv7_f*(a50_srv7_f-j))));
      srv7_sel_m(j) = (1/(1-gamma_srv7_m))*pow(((1-gamma_srv7_m)/gamma_srv7_m),gamma_srv7_m)*(mfexp(delta_srv7_m*gamma_srv7_m*(a50_srv7_m-j))/(1+mfexp(delta_srv7_m*(a50_srv7_m-j))));
    }
  }

    srv7_sel_f /= max(srv7_sel_f);   // standardize power function to max of 1.0
    srv7_sel_m /= max(srv7_sel_m);   // standardize power function to max of 1.0
    
// #####################################################################################################
// #####################################################################################################
// #####################################################################################################
// #####################################################################################################


// ###################################################################################################
// #####################################################################################################
// ####  UNUSED CODE ####################
// #####################################################################################################
// ###################################################################################################

       //  if(fsh2_sel_opt==3) fish2_sel(j)=1/ (1+mfexp(-delta_fish2*(double(j)-a50_fish4_f)));
       
       //  fish3_sel_f(j) = (1/(1+mfexp(-delta_fish3_f*(double(j)-a50_fish3_f))))*(1-(1/(1+mfexp(-gamma_fish3_f*(j-d50_fish3_f)))));
       //  fish3_sel_m(j) = (1/(1+mfexp(-delta_fish3_m*(double(j)-a50_fish3_m))))*(1-(1/(1+mfexp(-gamma_fish3_m*(j-d50_fish3_m))))); }

       //  d50_srv7_m=mfexp(log_d50_srv7_m);
       //  d50_srv7_f=mfexp(log_d50_srv7_f);     
       //  if(srv7_sel_opt==3) { 
       //     srv7_sel_f(j) = (1/(1+mfexp(-delta_srv7_f*(double(j)-a50_srv7_f))))*(1-(1/(1+mfexp(-gamma_srv7_f*(j-d50_srv7_f)))));
       //     srv7_sel_m(j) = (1/(1+mfexp(-delta_srv7_m*(double(j)-a50_srv7_m))))*(1-(1/(1+mfexp(-gamma_srv7_m*(j-d50_srv7_m))))); }


       //               srv1_sel_f=srv1_sel_f/max(srv1_sel_f);
       //              srv1_sel_m=srv1_sel_m/max(srv1_sel_m);
       //            srv2_sel_f=srv2_sel_f/max(srv2_sel_f);
       //             srv2_sel_m=srv2_sel_m/max(srv2_sel_m);
   
       /*    srv7_sel_f(j)= -(pow(j,a50_srv7_f)); // more flexible way to get convex or concave domeshaped with one parameter
             srv7_sel_m(j)= -(pow(j,a50_srv7_m)); } 
             srv7_sel_f= 1- srv7_sel_f/min(srv7_sel_f);
             srv7_sel_m= 1- srv7_sel_m/min(srv7_sel_m);
             srv7_sel_f /= max(srv7_sel_f);
             srv7_sel_m /= max(srv7_sel_m);  */

// #####################################################################################################
// #####################################################################################################
// #####################################################################################################
// #####################################################################################################
      
FUNCTION Get_Mortality_Rates

// Calculate mortality rates
  M_est = mfexp(logm);
 // Natural Mortality
 // ########################################################################################
  for (iyr=styr; iyr<=endyr; iyr++) // for all years prior to IFQ
   {
    for (j = 1 ; j<= nages; j++)
     {
      if((ph_Mdevs<1) && (ph_Mdevs_age<1))              // if no yearly/age M being estimated, then set M to base
       {
         natmort        = mfexp(logm);                   // setting natural mortality to arithmetic scale
         M_age(j)       = mfexp(logm);                   
         M_year(iyr)    = mfexp(logm);
         M(iyr,j)       = natmort;
       }
      if((ph_Mdevs>0) && (ph_Mdevs_age<1))              // if yearly M being estimated, then set M to base+devs
       {
         natmort        = mfexp(logm);
         M_age(j)       = mfexp(logm);                  
         M_year(iyr)    = mfexp(logm+log_M_devs(iyr));
         M(iyr,j)       = M_year(iyr);
       }
      if((ph_Mdevs<1) && (ph_Mdevs_age>0))              // if age M being estimated, then set M to base+age_devs
       {
         natmort        = mfexp(logm);
         M_age(j)       = mfexp(logm+log_M_devs_age(j));                
         M_year(iyr)    = mfexp(logm);
         M(iyr,j)       = M_age(j);
       }
      if((ph_Mdevs>0) && (ph_Mdevs_age>0))              // if age+year M being estimated, then set M to base+age_devs+year_devs
       {
         natmort        = mfexp(logm);
         M_age(j)       = mfexp(logm+log_M_devs_age(j));                
         M_year(iyr)    = mfexp(logm+log_M_devs(iyr));
         M(iyr,j)       = mfexp(logm+log_M_devs_age(j)+log_M_devs(iyr));
       }
     }
   }

 // Fishing Mortality
 // ########################################################################################
  Fmort_fish1    = mfexp(log_avg_F_fish1+log_F_devs_fish1);          //setting fishing mortaltiy to arithmetic scale
  
 for (iyr=styr; iyr<=endyr; iyr++) // for all years prior to IFQ
   {
    if(iyr<1963)
     {
      Fmort_fish3(iyr)=0.0;
     }
    if(iyr>=1963)
     {
      Fmort_fish3(iyr)= mfexp(log_avg_F_fish3 +log_F_devs_fish3(iyr));         //trawl mortality starts in 1963 (hence +3)
   }
  }

  hist_hal_F     = hist_hal_prop*mfexp(log_avg_F_fish1);             // optional historical fishing mortality for initial age comps
   
  for (iyr=styr; iyr<=1994; iyr++) // for all years prior to IFQ
   {
    for (j = 1 ; j<= nages; j++)
     {
      F_fish1_f(iyr,j) = Fmort_fish1(iyr) * fish1_sel_f(j);           // Getting fully selected fishing mortality
      F_fish1_m(iyr,j) = Fmort_fish1(iyr) * fish1_sel_m(j);           // Getting fully selected fishing mortality
      F_fish3_f(iyr,j) = Fmort_fish3(iyr) * fish3_sel_f(j);
      F_fish3_m(iyr,j) = Fmort_fish3(iyr) * fish3_sel_m(j);
     }
    }
    
  for (iyr=1995; iyr<=endyr; iyr++) //for years after IFQ implementation
   {
    for (j = 1 ; j<= nages; j++)
     {
      if((ph_ifq>0) && (ph_ifq_block2<1))  // if IFQ years and no recent time block sel est, then use fish4_sel for LL fleet
       {
        F_fish1_f(iyr,j) =  Fmort_fish1(iyr)*fish4_sel_f(j);
        F_fish1_m(iyr,j) =  Fmort_fish1(iyr)*fish4_sel_m(j);
       }
      if((ph_ifq>0) && (ph_ifq_block2>0))  // if IFQ years and  recent time block sel est, then use fish4_sel for LL fleet prior to start of recent block, fish5_sel after
       {
        if(iyr<yr_sel_chg_fish) // post-IFQ, pre recent time block
         {
          F_fish1_f(iyr,j) =  Fmort_fish1(iyr)*fish4_sel_f(j);
          F_fish1_m(iyr,j) =  Fmort_fish1(iyr)*fish4_sel_m(j);
         }
        if(iyr>=yr_sel_chg_fish) // post-IFQ, during recent time block
         {
          F_fish1_f(iyr,j) =  Fmort_fish1(iyr)*fish5_sel_f(j);
          F_fish1_m(iyr,j) =  Fmort_fish1(iyr)*fish5_sel_m(j);
         }
       }
      if((ph_ifq<1) && (ph_ifq_block2<1))     // if no post-IFQ sel est, then use fish1_sel
       {
        F_fish1_f(iyr,j) = Fmort_fish1(iyr)*fish1_sel_f(j);
        F_fish1_m(iyr,j) = Fmort_fish1(iyr)*fish1_sel_m(j);
       } 
        F_fish3_f(iyr,j) = Fmort_fish3(iyr) * fish3_sel_f(j);
        F_fish3_m(iyr,j) = Fmort_fish3(iyr) * fish3_sel_m(j);
      }
     }

 // Total Mortality
 // ########################################################################################
  Z_f            = F_fish1_f +F_fish3_f + M;                      // Fully selected total mortality
  Z_m            = F_fish1_m +F_fish3_m + M + mdelta;               // Fully selected total mortality
  S_f            = mfexp(-1.0*Z_f);                                    // Fully selected survival
  S_m            = mfexp(-1.0*Z_m);                                    // Fully selected survival
  S_f_mid        = mfexp(-0.5*Z_f);
  S_m_mid        = mfexp(-0.5*Z_m);
  
FUNCTION Get_Numbers_At_Age

  // bias ramp calc from Methot and Taylor (2011)

 b_yr_end=endyr-b_a50;
 
 for(y=(styr-nages+2);y<=endyr_rec_est;y++)
 {
 if(bias_ramp==1)
  {
     if(y<b_year_st || y==endyr_rec_est)
     {
      b(y)=0.0;
     }
     if(y>b_year_st && y<b_year_end)
     {
      b(y)=bmax*((y-b_year_st)/(b_year_end-b_year_st));
     }
     if(y>=b_year_end && y<=b_yr_end)
     {
      b(y)=bmax;
     }
     if(y>b_yr_end && y<endyr_rec_est)
     {
      b(y)=bmax*((endyr_rec_est-y)/(endyr_rec_est-b_yr_end));
     }
    }
 if(bias_ramp==0)
  {
   b(y)=1.0;
  }
 }

// Next two sections are based on Baranov catch equations
// Calculate Numbers at age
// Start year

 switch (SrType)
  {
    case 3: 
      int itmp;
  
    if(rec_like_type==3)
     {
      natage_f(styr,1)=mfexp(log_Rztemp)*mfexp((-sigr*sigr)/2)/2;
      natage_m(styr,1)=mfexp(log_Rztemp)*mfexp((-sigr*sigr)/2)/2;
     }

    if(rec_like_type<3)
     {
     if(sigma_R_early_switch==1 && styr<sigma_R_early_end)
      {
       natage_f(styr,1)=mfexp(log_mean_rec+log_rec_dev(styr)-b(styr)*sigma_R_early*sigma_R_early/2)/2;
       natage_m(styr,1)=mfexp(log_mean_rec+log_rec_dev(styr)-b(styr)*sigma_R_early*sigma_R_early/2)/2;       
      }
     else
      {
       natage_f(styr,1)=mfexp(log_mean_rec+log_rec_dev(styr)-b(styr)*sigr*sigr/2)/2;
       natage_m(styr,1)=mfexp(log_mean_rec+log_rec_dev(styr)-b(styr)*sigr*sigr/2)/2;
      }

    // ################ Using the rec_devs in years<styr to determine initial abundance
    for (j=2;j<nages;j++)
     {
      itmp = styr+1-j;
     if(sigma_R_early_switch==1 && itmp<sigma_R_early_end)
      {
       natage_f(styr,j) =( mfexp(log_mean_rec  - (M(styr,j)+hist_hal_F*fish1_sel_f(j)) * double(j-1)+ log_rec_dev(itmp)-b(itmp)*sigma_R_early*sigma_R_early/2))/2; 
       natage_m(styr,j) = (mfexp(log_mean_rec  - (M(styr,j)+mdelta+hist_hal_F*fish1_sel_m(j)) * double(j-1)+ log_rec_dev(itmp)-b(itmp)*sigma_R_early*sigma_R_early/2))/2;      
      }
     else
      {
       natage_f(styr,j) =( mfexp(log_mean_rec  - (M(styr,j)+hist_hal_F*fish1_sel_f(j)) * double(j-1)+ log_rec_dev(itmp)-b(itmp)*sigr*sigr/2))/2; 
       natage_m(styr,j) = (mfexp(log_mean_rec  - (M(styr,j)+mdelta+hist_hal_F*fish1_sel_m(j)) * double(j-1)+ log_rec_dev(itmp)-b(itmp)*sigr*sigr/2))/2; 
     }
     }
     
    // ########## calc for plus group
      natage_f(styr,nages)      = (mfexp(log_mean_rec - (M(styr,nages-1)+hist_hal_F*fish1_sel_f(nages-1)) * (nages-1))/ (1. - exp(-(M(styr,nages-1)+hist_hal_F*fish1_sel_f(nages-1)))))/2;
      natage_m(styr,nages)      = (mfexp(log_mean_rec - (M(styr,nages-1)+mdelta+hist_hal_F*fish1_sel_m(nages-1)) * (nages-1))/ (1. - exp(-(M(styr,nages-1)+mdelta+hist_hal_F*fish1_sel_m(nages-1)))))/2;
     }
    break;
   }

 // FIll in abundance matrix until last year of recruit devs
  for ( i=styr;i <= endyr_rec;i++)
   {
     if(sigma_R_early_switch==1 && i<sigma_R_early_end)
      {
       natage_f(i,1)           = mfexp(log_rec_dev(i) + log_mean_rec-b(i)*sigma_R_early*sigma_R_early/2 )/2;  //fill in recruitment in current year
       natage_m(i,1)           = mfexp(log_rec_dev(i) + log_mean_rec-b(i)*sigma_R_early*sigma_R_early/2 )/2;      
      }
     if(sigma_R_early_switch==0 || i>=sigma_R_early_end)
      {
       natage_f(i,1)           = mfexp(log_rec_dev(i) + log_mean_rec-b(i)*sigr*sigr/2 )/2;  //fill in recruitment in current year
       natage_m(i,1)           = mfexp(log_rec_dev(i) + log_mean_rec-b(i)*sigr*sigr/2 )/2;
      }
     sam_rec(i)              = natage_f(i,1)+natage_m(i,1); // OjO
     natage_f(i+1)(2,nages)  = ++elem_prod(natage_f(i)(1,nages-1),S_f(i)(1,nages-1));       // Fill in N in following year, Survival from prev since uses i and N used i+1
     natage_m(i+1)(2,nages)  = ++elem_prod(natage_m(i)(1,nages-1),S_m(i)(1,nages-1));       // Fill in N in following year, Survival from prev since uses i and N used i+1
     natage_f(i+1,nages)    += natage_f(i,nages)*S_f(i,nages);
     natage_m(i+1,nages)    += natage_m(i,nages)*S_m(i,nages);                              // Plus group calc, add fish in plus group from last year that survive
     Sp_Biom(i)  = elem_prod(natage_f(i),pow(S_f(i),spawn_fract)) * weight_maturity_prod_f(i); 
   }
  
 // FIll in abundance matrix for years with no rec_devs estimated, set recruit equal to Rave
   for(i=endyr_rec+1;i<endyr;i++)
    {
     natage_f(i,1)=mfexp(log_mean_rec-sigr*sigr/2)/2;
     natage_m(i,1)=mfexp(log_mean_rec-sigr*sigr/2)/2;
     sam_rec(i)            = natage_f(i,1)+natage_m(i,1); // OjO
     natage_f(i+1)(2,nages)  = ++elem_prod(natage_f(i)(1,nages-1),S_f(i)(1,nages-1));       // Following year
     natage_m(i+1)(2,nages)  = ++elem_prod(natage_m(i)(1,nages-1),S_m(i)(1,nages-1));       // Following year
     natage_f(i+1,nages)    += natage_f(i,nages)*S_f(i,nages);
     natage_m(i+1,nages)    += natage_m(i,nages)*S_m(i,nages);
     Sp_Biom(i)  = elem_prod(natage_f(i),pow(S_f(i),spawn_fract)) * weight_maturity_prod_f(i); 
    }
    
  if(rec_like_type==1)
   {
     natage_f(endyr,1)=(mfexp(log_mean_rec+log_rec_dev(endyr)))/2;
     natage_m(endyr,1)=(mfexp(log_mean_rec+log_rec_dev(endyr)))/2;
   }
  else
   {
    natage_f(endyr,1)         = mfexp(log_mean_rec )/2; 
    natage_m(endyr,1)         = mfexp(log_mean_rec )/2;
   }
   
  Sp_Biom(endyr)  = elem_prod(natage_f(endyr),pow(S_f(endyr),spawn_fract)) * weight_maturity_prod_f(i); 


  //Get numbers at length

  for (i=styr;i<=endyr;i++)
   {                     
    num_len_m(i)   = (natage_m(i))* size_age_m(i);                                              
    num_len_f(i)   = (natage_f(i))* size_age_f(i);
   }                                            

   
FUNCTION Get_Catch_at_Age
// Calculate catch at age
  pred_catch_fish1.initialize();
  pred_catch_fish3.initialize();
  for (iyr=styr; iyr<=endyr; iyr++)
  {
    catage_fish1_m(iyr) = elem_div(elem_prod(elem_prod(natage_m(iyr),F_fish1_m(iyr)),(1.-S_m(iyr))),Z_m(iyr));
    catage_fish1_f(iyr) = elem_div(elem_prod(elem_prod(natage_f(iyr),F_fish1_f(iyr)),(1.-S_f(iyr))),Z_f(iyr));
    pred_catch_fish1(iyr) = elem_div(elem_prod(elem_prod(natage_f(iyr),F_fish1_f(iyr)),(1.-S_f(iyr))),Z_f(iyr))*weight_f(iyr)+elem_div(elem_prod(elem_prod(natage_m(iyr),F_fish1_m(iyr)),(1.-S_m(iyr))),Z_m(iyr))*weight_m(iyr);
    catage_fish3_m(iyr) = elem_div(elem_prod(elem_prod(natage_m(iyr),F_fish3_m(iyr)),(1.-S_m(iyr))),Z_m(iyr));
    catage_fish3_f(iyr) = elem_div(elem_prod(elem_prod(natage_f(iyr),F_fish3_f(iyr)),(1.-S_f(iyr))),Z_f(iyr));
    pred_catch_fish3(iyr) = elem_div(elem_prod(elem_prod(natage_f(iyr),F_fish3_f(iyr)),(1.-S_f(iyr))),Z_f(iyr))*weight_f(iyr)+elem_div(elem_prod(elem_prod(natage_m(iyr),F_fish3_m(iyr)),(1.-S_m(iyr))),Z_m(iyr))*weight_m(iyr);
  }
  
    pred_catch=(pred_catch_fish1+pred_catch_fish3);

FUNCTION Get_Dependent_Vars
  for (i=styr;i<=endyr;i++)
  {
    pred_rec(i) = natage_f(i,1)+natage_m(i,1);                  // Setting up results based on estimated paramters
    tot_biom(i) = weight_f(i) * natage_f(i)+weight_m(i)*natage_m(i);          // Total biomass results
    spawn_biom(i) = weight_maturity_prod_f(i)*natage_f(i) ;                     // Spawning biomass result
  }

  avg_rec        = mean(pred_rec);
  Depletion      = spawn_biom(endyr)/spawn_biom(styr);                         // Depletion
  spbiom_trend   = spawn_biom(endyr)/spawn_biom(endyr-1);
  
FUNCTION Get_Predicted_Values

 // Transform estimated catchabilities
 // ###############################################################################################################################################################################
 
   q_srv1         = mfexp(log_q_srv1);                                        // Survey catchability at arithmetic scale
   q_srv2         = mfexp(log_q_srv2);                                        // Survey catchability at arithmetic scale
   q_srv3         = mfexp(log_q_srv3);                                        // Survey catchability at arithmetic scale
   q_srv4         = mfexp(log_q_srv4);                                        // Survey catchability at arithmetic scale
   q_srv5         = mfexp(log_q_srv5);                                        // Survey catchability at arithmetic scale
   q_srv6         = mfexp(log_q_srv6);                                        // Survey catchability at arithmetic scale
   q_srv7         = mfexp(log_q_srv7);                                        // Survey catchability at arithmetic scale
   q_srv8         = mfexp(log_q_srv8);                                        // Survey catchability at arithmetic scale
   q_srv9         = mfexp(log_q_srv2);                                        // Survey catchability at arithmetic scale
   q_LL_fish_recent = mfexp(log_q_LL_fish_recent);
   q_LL_srvy_recent = mfexp(log_q_LL_srvy_recent);
 // Survey and CPUE Calcs
 // ###############################################################################################################################################################################

  for (i=styr;i<=endyr;i++)
   { 
    pred_srv6(i) =    q_srv6 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(weight_f(i),fish2_sel)+ elem_prod(S_m_mid(i),natage_m(i))*elem_prod(weight_m(i),fish2_sel));   // Predicted Survey biomass
    pred_srv7(i) = 2*(q_srv7 * (1-prop_m2(i))*elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv7_sel_f,weight_f(i))+ q_srv7*prop_m2(i)*elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv7_sel_m,weight_m(i)));   // Predicted Survey biomass
    pred_srv4(i) = 2*(q_srv2 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv2_sel_f)+q_srv2 *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv2_sel_m));
   }
   
  for (i=1979;i<=1989;i++)
   {
    pred_srv2(i) = 2*(q_srv2* (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv2_sel_f,weight_f(i)))+q_srv2 * prop_m(i)*(elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv2_sel_m,weight_m(i)))); 
    pred_srv4(i) = 2*(q_srv2* (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv2_sel_f)+q_srv2 *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv2_sel_m));
   } 
       
  for (i=1990;i<=1994;i++)
   {
    pred_srv2(i) = 2*(q_srv9 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv2_sel_f,weight_f(i)))+q_srv9 *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv2_sel_m,weight_m(i))));    // Predicted Survey biomass
    pred_srv4(i) = 2*(q_srv9 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv2_sel_f)+q_srv9 *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv2_sel_m));
   }   

 // Domestic LL Survey Calcs
 // ###############################################################################################################################################################################

  for (i=styr;i<=endyr;i++)
   {
    if(i<yr_sel_chg_srv1)
     {   
       pred_srv1(i) = 2*(q_srv1 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv1_sel_f,weight_f(i)))+q_srv1*prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv1_sel_m,weight_m(i))));   // Predicted Survey biomass
       pred_srv3(i) = 2*(q_srv1 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv1_sel_f)+q_srv1 *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv1_sel_m));   // Predicted Survey biomass
     } 
    if(i>=yr_sel_chg_srv1)
     {
      if(ph_LL_block2>0 && ph_q_LL_srv_rec>0)  // if estimating new selectivity for recent time block then use srv10 selectivity AND NEW Q
       {   
        pred_srv1(i) = 2*(q_LL_srvy_recent * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv10_sel_f,weight_f(i)))+q_LL_srvy_recent*prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv10_sel_m,weight_m(i))));   // Predicted Survey biomass
        pred_srv3(i) = 2*(q_LL_srvy_recent * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv10_sel_f)+q_LL_srvy_recent *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv10_sel_m));   // Predicted Survey biomass
       }
      if(ph_LL_block2>0 && ph_q_LL_srv_rec<1)  // if estimating new selectivity for recent time block then use srv10 selectivity AND old Q
       {   
        pred_srv1(i) = 2*(q_srv1 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv10_sel_f,weight_f(i)))+q_srv1*prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv10_sel_m,weight_m(i))));   // Predicted Survey biomass
        pred_srv3(i) = 2*(q_srv1 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv10_sel_f)+q_srv1 *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv10_sel_m));   // Predicted Survey biomass
       } 
      if(ph_LL_block2<1 && ph_q_LL_srv_rec<1) // if not estimating new selectivity for recent time block then use srv1 selectivity
       {   
        pred_srv1(i) = 2*(q_srv1 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv1_sel_f,weight_f(i)))+q_srv1*prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv1_sel_m,weight_m(i))));   // Predicted Survey biomass
        pred_srv3(i) = 2*(q_srv1 * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv1_sel_f)+q_srv1 *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv1_sel_m));   // Predicted Survey biomass
       }
      if(ph_LL_block2<1 && ph_q_LL_srv_rec>0) // if not estimating new selectivity for recent time block then use srv1 selectivity
       {   
        pred_srv1(i) = 2*(q_LL_srvy_recent * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*elem_prod(srv1_sel_f,weight_f(i)))+q_LL_srvy_recent*prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(srv1_sel_m,weight_m(i))));   // Predicted Survey biomass
        pred_srv3(i) = 2*(q_LL_srvy_recent * (1-prop_m(i))*(elem_prod(S_f_mid(i),natage_f(i))*srv1_sel_f)+q_LL_srvy_recent *prop_m(i)* (elem_prod(S_m_mid(i),natage_m(i))*srv1_sel_m));   // Predicted Survey biomass
       }
     }
   }
   
 // Domestic LL Fishery CPUE Calcs
 // ###############################################################################################################################################################################

  // Pre-IFQ
  for (i=styr;i<=1994;i++) // pre-IFQ, use q_srv5
   {
    pred_srv5(i) = q_srv5 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish1_sel_f,weight_f(i)))+q_srv5 * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish1_sel_m,weight_m(i)));    // Predicted Survey biomass
   }
  // Post-IFQ
  for (i=1995;i<=endyr;i++)  // post-IFQ use q_srv8
   { 
      if((ph_ifq>0) && (ph_ifq_block2<1) && ph_q_IFQ_rec<1)  // if IFQ years and no recent time block sel est, then use fish4_sel for LL fleet
       {
         pred_srv5(i) = q_srv8 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish4_sel_f,weight_f(i)))+q_srv8 * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish4_sel_m,weight_m(i)));     // Predicted Survey biomass
       }


      if((ph_ifq>0) && (ph_ifq_block2>0) && ph_q_IFQ_rec>0)  // if IFQ years and  recent time block sel+q est, then use fish4_sel and q_srv8 for LL fleet prior to start of recent block, fish5_sel and q_LL_fish_recent after
       {
        if(i<yr_sel_chg_fish) // post-IFQ, pre recent time block
         {
          pred_srv5(i) = q_srv8 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish4_sel_f,weight_f(i)))+q_srv8 * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish4_sel_m,weight_m(i)));     // Predicted Survey biomass
         }
        if(i>=yr_sel_chg_fish) // post-IFQ, during recent time block
         {
          pred_srv5(i) = q_LL_fish_recent * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish5_sel_f,weight_f(i)))+q_LL_fish_recent * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish5_sel_m,weight_m(i)));     // Predicted Survey biomass
         }
       }

      if((ph_ifq>0) && (ph_ifq_block2>0) && ph_q_IFQ_rec<1)  // if IFQ years and  recent time block sel est only, then use fish4_sel and q_srv8 for LL fleet prior to start of recent block, fish5_sel and q_srv8 after
       {
        if(i<yr_sel_chg_fish) // post-IFQ, pre recent time block
         {
          pred_srv5(i) = q_srv8 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish4_sel_f,weight_f(i)))+q_srv8 * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish4_sel_m,weight_m(i)));     // Predicted Survey biomass
         }
        if(i>=yr_sel_chg_fish) // post-IFQ, during recent time block
         {
          pred_srv5(i) = q_srv8 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish5_sel_f,weight_f(i)))+q_srv8 * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish5_sel_m,weight_m(i)));     // Predicted Survey biomass
         }
       }

      if((ph_ifq>0) && (ph_ifq_block2<1) && ph_q_IFQ_rec>0)  // if IFQ years and  recent time block q est only, then use fish4_sel and q_srv8 for LL fleet prior to start of recent block, fish4_sel and q_LL_fish_recent after
       {
        if(i<yr_sel_chg_fish) // post-IFQ, pre recent time block
         {
          pred_srv5(i) = q_srv8 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish4_sel_f,weight_f(i)))+q_srv8 * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish4_sel_m,weight_m(i)));     // Predicted Survey biomass
         }
        if(i>=yr_sel_chg_fish) // post-IFQ, during recent time block
         {
          pred_srv5(i) = q_LL_fish_recent * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish4_sel_f,weight_f(i)))+q_LL_fish_recent * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish4_sel_m,weight_m(i)));     // Predicted Survey biomass
         }
       }


      if((ph_ifq<1) && (ph_ifq_block2<1) && ph_q_IFQ_rec<1)     // if no post-IFQ sel est, then use fish1_sel and q_srv5
       {
         pred_srv5(i) = q_srv5 * (elem_prod(S_f_mid(i),natage_f(i))*elem_prod(fish1_sel_f,weight_f(i)))+q_srv5 * (elem_prod(S_m_mid(i),natage_m(i))*elem_prod(fish1_sel_m,weight_m(i)));
       }
    }

 // ###############################################################################################################################################################################
 

 // Fishery Age Comp Calcs
 // ###############################################################################################################################################################################
  for (i=1;i<=nyrs_fish1_age;i++)
   {
    eac_fish1(i)  = ((catage_fish1_m(yrs_fish1_age(i))/sum(catage_fish1_m(yrs_fish1_age(i))))+(catage_fish1_f(yrs_fish1_age(i))/sum(catage_fish1_f(yrs_fish1_age(i)))))/2* ageage;                                                // Predicted Fishery age comps
    eac_fish1(i) /=sum(eac_fish1(i));
   }

 // Fishery Size Comp Calcs
 // ###############################################################################################################################################################################                                                           // Second Predicted Fishery size comps for 80s and 90s

  for (i=1;i<=nyrs_fish1_size;i++)
   {                    
    esc_fish1_m(i)  = catage_fish1_m(yrs_fish1_size(i))/sum(catage_fish1_m(yrs_fish1_size(i)))* size_age_m(yrs_fish1_size(i));    
    esc_fish1_f(i)  = catage_fish1_f(yrs_fish1_size(i))/sum(catage_fish1_f(yrs_fish1_size(i)))* size_age_f(yrs_fish1_size(i));
   }
   
  for (i=1;i<=nyrs_fish3_size;i++)
   {
    esc_fish3_m(i)  = (catage_fish3_m(yrs_fish3_size(i))/sum(catage_fish3_m(yrs_fish3_size(i))))* size_age_m(yrs_fish3_size(i));                                              // Second Predicted Fishery size comps for 80s and 90s
    esc_fish3_f(i)  = (catage_fish3_f(yrs_fish3_size(i))/sum(catage_fish3_f(yrs_fish3_size(i))))* size_age_f(yrs_fish3_size(i));
   } 
                                           // Second Predicted Fishery size comps for 80s and 90s
   
  for (i=1;i<=nyrs_fish2_size;i++)
   {                     // Lets you use a second matrix for part of it
    esc_fish2(i)  = ((1-prop_m(yrs_fish2_size(i)))*elem_prod(fish2_sel,natage_f(yrs_fish2_size(i)))+prop_m(yrs_fish2_size(i))*elem_prod(fish2_sel,natage_m(yrs_fish2_size(i))))* sizeage_all;                          // Predicted Survey age comps
    esc_fish2(i) /=sum(esc_fish2(i));
   }


 // Survey Age Comp Calcs
 // ###############################################################################################################################################################################
  for (i=1;i<=nyrs_srv1_age;i++)  // account for potential selectivity change in the DOM LL Survey
   {
    if(yrs_srv1_age(i)<yr_sel_chg_srv1)
     {   
      eac_srv1(i)  = ((1-prop_m(yrs_srv1_age(i)))*elem_prod(srv1_sel_f,natage_f(yrs_srv1_age(i)))+prop_m(yrs_srv1_age(i))*elem_prod(srv1_sel_m,natage_m(yrs_srv1_age(i))))* ageage;                         // Predicted Survey age comps
      eac_srv1(i) /=sum(eac_srv1(i));
     } 
    if(yrs_srv1_age(i)>=yr_sel_chg_srv1)
     {
      if(ph_LL_block2>0)  // if estimating new selectivity for recent time block then use srv10 selectivity
       {   
        eac_srv1(i)  = ((1-prop_m(yrs_srv1_age(i)))*elem_prod(srv10_sel_f,natage_f(yrs_srv1_age(i)))+prop_m(yrs_srv1_age(i))*elem_prod(srv10_sel_m,natage_m(yrs_srv1_age(i))))* ageage;                         // Predicted Survey age comps
        eac_srv1(i) /=sum(eac_srv1(i));
       } 
      if(ph_LL_block2<1) // if not estimating new selectivity for recent time block then use srv1 selectivity
       {   
        eac_srv1(i)  = ((1-prop_m(yrs_srv1_age(i)))*elem_prod(srv1_sel_f,natage_f(yrs_srv1_age(i)))+prop_m(yrs_srv1_age(i))*elem_prod(srv1_sel_m,natage_m(yrs_srv1_age(i))))* ageage;                         // Predicted Survey age comps
        eac_srv1(i) /=sum(eac_srv1(i));
       }
     }
   }


  for (i=1;i<=nyrs_srv2_age;i++)
   {
    eac_srv2(i)  = (elem_prod(srv2_sel_f,natage_f(yrs_srv2_age(i)))+elem_prod(srv2_sel_m,natage_m(yrs_srv2_age(i))))*ageage;                        // Predicted Survey age comps
    eac_srv2(i) /=sum(eac_srv2(i));
   }
   
  for (i=1;i<=nyrs_srv7_age;i++)
   {
    eac_srv7(i)  = ((1-prop_m(yrs_srv7_age(i)))*elem_prod(srv7_sel_f,natage_f(yrs_srv7_age(i)))+prop_m(yrs_srv7_age(i))*elem_prod(srv7_sel_m,natage_m(yrs_srv7_age(i))))* ageage;                         // Predicted Survey age comps
    eac_srv7(i) /=sum(eac_srv7(i));
   }


 // Survey Size Comp Calcs
 // ###############################################################################################################################################################################                                                           // Second Predicted Fishery size comps for 80s and 90s
  for ( i=1;i<=nyrs_srv1_size;i++) //old length_age key
   {
    if(yrs_srv1_size(i)<yr_sel_chg_srv1)
     {   
      esc_srv1_m(i)  = elem_prod(srv1_sel_m,natage_m(yrs_srv1_size(i)))* size_age_m(yrs_srv1_size(i));        // Predicted Survey size comps (not used in POP model)
      esc_srv1_m(i)  /=sum(esc_srv1_m(i));
      esc_srv1_f(i)  = elem_prod(srv1_sel_f,natage_f(yrs_srv1_size(i))) * size_age_f(yrs_srv1_size(i));        // Predicted Survey size comps (not used in POP model)
      esc_srv1_f(i)  /=sum(esc_srv1_f(i));
     } 
    if(yrs_srv1_size(i)>=yr_sel_chg_srv1)
     {
      if(ph_LL_block2>0)  // if estimating new selectivity for recent time block then use srv10 selectivity
       {   
        esc_srv1_m(i)  = elem_prod(srv10_sel_m,natage_m(yrs_srv1_size(i))) * size_age_m(yrs_srv1_size(i));        // Predicted Survey size comps (not used in POP model)
        esc_srv1_m(i)  /=sum(esc_srv1_m(i));
        esc_srv1_f(i)  = elem_prod(srv10_sel_f,natage_f(yrs_srv1_size(i))) * size_age_f(yrs_srv1_size(i));        // Predicted Survey size comps (not used in POP model)
        esc_srv1_f(i)  /=sum(esc_srv1_f(i));
       } 
      if(ph_LL_block2<1) // if not estimating new selectivity for recent time block then use srv1 selectivity
       {   
        esc_srv1_m(i)  = elem_prod(srv1_sel_m,natage_m(yrs_srv1_size(i)))* size_age_m(yrs_srv1_size(i));        // Predicted Survey size comps (not used in POP model)
        esc_srv1_m(i)  /=sum(esc_srv1_m(i));
        esc_srv1_f(i)  = elem_prod(srv1_sel_f,natage_f(yrs_srv1_size(i))) * size_age_f(yrs_srv1_size(i));        // Predicted Survey size comps (not used in POP model)
        esc_srv1_f(i)  /=sum(esc_srv1_f(i));
       }
     }
   }

  for ( i=1;i<=nyrs_srv2_size;i++)
   {
    esc_srv2_m(i)  = elem_prod(srv2_sel_m,natage_m(yrs_srv2_size(i)))*  size_age_m(yrs_srv2_size(i));        // Predicted Survey size comps (not used in POP model)
    esc_srv2_m(i)  /=sum(esc_srv2_m(i)); 
    esc_srv2_f(i)  = elem_prod(srv2_sel_f,natage_f(yrs_srv2_size(i))) * size_age_f(yrs_srv2_size(i));        // Predicted Survey size comps (not used in POP model)
    esc_srv2_f(i)  /=sum(esc_srv2_f(i));
   }
   
  for ( i=1;i<=nyrs_srv7_size;i++)
   {
    esc_srv7_m(i)  = elem_prod(srv7_sel_m,natage_m(yrs_srv7_size(i)))*size_age_m(yrs_srv7_size(i));        // Predicted Survey size comps (not used in POP model)
    esc_srv7_f(i)  = elem_prod(srv7_sel_f,natage_f(yrs_srv7_size(i)))* size_age_f(yrs_srv7_size(i));       // Predicted Survey size comps (not used in POP model)
    esc_srv7_f(i)  /=sum(esc_srv7_f(i)); 
    esc_srv7_m(i)  /=sum(esc_srv7_m(i));
   }

FUNCTION compute_spr_rates

  //Compute SPR Rates and add them to the likelihood for Females 
  fratio = Fmort_fish1(endyr)/(Fmort_fish1(endyr)+Fmort_fish3(endyr));

    mF50=mfexp(log_mF50);
    mF40=mfexp(log_mF40);
    mF35=mfexp(log_mF35);
    
  // Scale F-spr rates to be on full-selected values
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
       F50  = mF50*max(fish4_sel_f);
       F40  = mF40*max(fish4_sel_f);
       F35  = mF35*max(fish4_sel_f);
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
      {
       F50  = mF50*max(fish5_sel_f);
       F40  = mF40*max(fish5_sel_f);
       F35  = mF35*max(fish5_sel_f);
      }   
   if(ph_ifq<1) // no post-IFQ sel est
     {
      F50  = mF50*max(fish1_sel_f);
      F40  = mF40*max(fish1_sel_f);
      F35  = mF35*max(fish1_sel_f);
     }
  
  SB0   = 0;
  SBF50 = 0;
  SBF40 = 0;
  SBF35 = 0;


  for (i=1;i<=4;i++) // i index is for the various SPR fractions (0, 0.5, 0.4, 0.35)
   {
    Nspr(i,1)=1.;
   }
  // ########################## Use average estimate M 
    for (j=2;j<nages;j++)
     {
      Nspr(1,j)=Nspr(1,j-1)*mfexp(-1.*natmort);
      
      if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
       {
        Nspr(2,j)=Nspr(2,j-1)*mfexp(-1.*(natmort+fratio*mF50*fish4_sel_f(j-1)+(1-fratio)*mF50*fish3_sel_f(j-1)));
        Nspr(3,j)=Nspr(3,j-1)*mfexp(-1.*(natmort+fratio*mF40*fish4_sel_f(j-1)+(1-fratio)*mF40*fish3_sel_f(j-1)));
        Nspr(4,j)=Nspr(4,j-1)*mfexp(-1.*(natmort+fratio*mF35*fish4_sel_f(j-1)+(1-fratio)*mF35*fish3_sel_f(j-1)));
       }
      if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
       {
        Nspr(2,j)=Nspr(2,j-1)*mfexp(-1.*(natmort+fratio*mF50*fish5_sel_f(j-1)+(1-fratio)*mF50*fish3_sel_f(j-1)));
        Nspr(3,j)=Nspr(3,j-1)*mfexp(-1.*(natmort+fratio*mF40*fish5_sel_f(j-1)+(1-fratio)*mF40*fish3_sel_f(j-1)));
        Nspr(4,j)=Nspr(4,j-1)*mfexp(-1.*(natmort+fratio*mF35*fish5_sel_f(j-1)+(1-fratio)*mF35*fish3_sel_f(j-1)));
       }   
      if(ph_ifq<1) // no post-IFQ sel est
       {
        Nspr(2,j)=Nspr(2,j-1)*mfexp(-1.*(natmort+fratio*mF50*fish1_sel_f(j-1)+(1-fratio)*mF50*fish3_sel_f(j-1)));
        Nspr(3,j)=Nspr(3,j-1)*mfexp(-1.*(natmort+fratio*mF40*fish1_sel_f(j-1)+(1-fratio)*mF40*fish3_sel_f(j-1)));
        Nspr(4,j)=Nspr(4,j-1)*mfexp(-1.*(natmort+fratio*mF35*fish1_sel_f(j-1)+(1-fratio)*mF35*fish3_sel_f(j-1)));
       }
      }
  
      Nspr(1,nages)=Nspr(1,nages-1)*mfexp(-1.*natmort)/(1.-mfexp(-1.*natmort));

      if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
       {
        Nspr(2,nages)=Nspr(2,nages-1)*mfexp(-1.* (natmort+fratio*mF50*fish4_sel_f(nages-1)+(1-fratio)*mF50*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF50*fish4_sel_f(nages)+(1-fratio)*mF50*fish3_sel_f(nages))));
        Nspr(3,nages)=Nspr(3,nages-1)*mfexp(-1.* (natmort+fratio*mF40*fish4_sel_f(nages-1)+(1-fratio)*mF40*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF40*fish4_sel_f(nages)+(1-fratio)*mF40*fish3_sel_f(nages))));
        Nspr(4,nages)=Nspr(4,nages-1)*mfexp(-1.* (natmort+fratio*mF35*fish4_sel_f(nages-1)+(1-fratio)*mF35*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF35*fish4_sel_f(nages)+(1-fratio)*mF35*fish3_sel_f(nages))));
       }
      if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
       {
        Nspr(2,nages)=Nspr(2,nages-1)*mfexp(-1.* (natmort+fratio*mF50*fish5_sel_f(nages-1)+(1-fratio)*mF50*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF50*fish5_sel_f(nages)+(1-fratio)*mF50*fish3_sel_f(nages))));
        Nspr(3,nages)=Nspr(3,nages-1)*mfexp(-1.* (natmort+fratio*mF40*fish5_sel_f(nages-1)+(1-fratio)*mF40*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF40*fish5_sel_f(nages)+(1-fratio)*mF40*fish3_sel_f(nages))));
        Nspr(4,nages)=Nspr(4,nages-1)*mfexp(-1.* (natmort+fratio*mF35*fish5_sel_f(nages-1)+(1-fratio)*mF35*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF35*fish5_sel_f(nages)+(1-fratio)*mF35*fish3_sel_f(nages))));
       } 
      if(ph_ifq<1) // no post-IFQ sel est
       {
        Nspr(2,nages)=Nspr(2,nages-1)*mfexp(-1.* (natmort+fratio*mF50*fish1_sel_f(nages-1)+(1-fratio)*mF50*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF50*fish1_sel_f(nages)+(1-fratio)*mF50*fish3_sel_f(nages))));
        Nspr(3,nages)=Nspr(3,nages-1)*mfexp(-1.* (natmort+fratio*mF40*fish1_sel_f(nages-1)+(1-fratio)*mF40*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF40*fish1_sel_f(nages)+(1-fratio)*mF40*fish3_sel_f(nages))));
        Nspr(4,nages)=Nspr(4,nages-1)*mfexp(-1.* (natmort+fratio*mF35*fish1_sel_f(nages-1)+(1-fratio)*mF35*fish3_sel_f(nages-1)))/ (1.-mfexp(-1.*(natmort+fratio*mF35*fish1_sel_f(nages)+(1-fratio)*mF35*fish3_sel_f(nages))));
       }
  
    for (j=1;j<=nages;j++)
     {
      // Kill them off till (spawn_fract)
       SB0    += Nspr(1,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*natmort);
       
      if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
       {
        SBF50  += Nspr(2,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF50*fish4_sel_f(j)+(1-fratio)*mF50*fish3_sel_f(j)));
        SBF40  += Nspr(3,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF40*fish4_sel_f(j)+(1-fratio)*mF40*fish3_sel_f(j)));
        SBF35  += Nspr(4,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF35*fish4_sel_f(j)+(1-fratio)*mF35*fish3_sel_f(j)));
       }
      if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
       {
        SBF50  += Nspr(2,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF50*fish5_sel_f(j)+(1-fratio)*mF50*fish3_sel_f(j)));
        SBF40  += Nspr(3,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF40*fish5_sel_f(j)+(1-fratio)*mF40*fish3_sel_f(j)));
        SBF35  += Nspr(4,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF35*fish5_sel_f(j)+(1-fratio)*mF35*fish3_sel_f(j)));
       }
      if(ph_ifq<1) // no post-IFQ sel est
       {
        SBF50  += Nspr(2,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF50*fish1_sel_f(j)+(1-fratio)*mF50*fish3_sel_f(j)));
        SBF40  += Nspr(3,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF40*fish1_sel_f(j)+(1-fratio)*mF40*fish3_sel_f(j)));
        SBF35  += Nspr(4,j)*weight_maturity_prod_f(endyr,j)*mfexp(-spawn_fract*(natmort+fratio*mF35*fish1_sel_f(j)+(1-fratio)*mF35*fish3_sel_f(j)));
       }
      }  // end NAGES loop
   
  sprpen    = 100.*square(SBF50/SB0-0.5);
  sprpen   += 100.*square(SBF40/SB0-0.4);
  sprpen   += 100.*square(SBF35/SB0-0.35);
  B40       = 0.5*SBF40*mean(pred_rec(1979,endyr-recage));
  
FUNCTION Calc_priors

// Calculate priors
    priors.initialize();
    if (active(log_sigr))
      priors(1)    = wt_sigr_prior*square(log((sigr/mfexp(sigrprior))))/(2.*square(cvsigrprior));
    if (active(log_q_srv1))
      priors(2)    =wt_q_priors* square(log_q_srv1-log(q_srv1prior))/(2.*square(cvq_srv1prior));
    if (active(steepness))
      priors(3)    = square(log(steepness/steep_prior))/(2.*cv_steep_prior); // not used in POP model
    if (active(logm))
      priors(4)    = wt_M_prior*square(logm-log(mprior))/(2.*square(cvmprior));
    if (active(log_q_srv2))
      priors(5)    = wt_q_priors*square(log_q_srv2-log(q_srv2prior))/(2.*square(cvq_srv2prior));
    if (active(log_q_srv3))
      priors(6)    = wt_q_priors*square(log_q_srv3-log(q_srv3prior))/(2.*square(cvq_srv3prior));
    if (active(log_q_srv4))
      priors(7)    = wt_q_priors*square(log_q_srv4-log(q_srv4prior))/(2.*square(cvq_srv4prior));
    if (active(log_q_srv5))
      priors(8)    = wt_q_priors*square(log_q_srv5-log(q_srv5prior))/(2.*square(cvq_srv5prior));
    if (active(log_q_srv6))
      priors(9)    = wt_q_priors*square(log_q_srv6-log(q_srv6prior))/(2.*square(cvq_srv6prior));
    if (active(log_q_srv7))
      priors(10)    = wt_q_priors*square(log_q_srv7-log(q_srv7prior))/(2.*square(cvq_srv7prior));
    if (active(log_q_srv8))
      priors(11)    = wt_q_priors*square(log_q_srv8-log(q_srv8prior))/(2.*square(cvq_srv8prior));
    if (active(log_q_LL_fish_recent))
      priors(12)    = wt_q_priors*square(log_q_LL_fish_recent-log(q_srv8prior))/(2.*square(cvq_srv8prior)); //use early IFQ fishery q prior for recent time block fishery q prior
    if (active(log_q_LL_srvy_recent))
      priors(13)    = wt_q_priors*square(log_q_LL_srvy_recent-log(q_srv1prior))/(2.*square(cvq_srv1prior)); // use LL srvy q prior for recent time block LL survey q prior
FUNCTION Surv_Likelihood
 // Calculate likelihood for survey biomass
  surv_like.initialize();
  for (i=1; i<=nyrs_srv1; i++)  {   ii=yrs_srv1(i);
  surv_like(1) += square((log(obs_srv1_biom(i)+0.0001)-log(pred_srv1(ii)+0.0001) ))/ (2.*square(obs_srv1_se(i)/obs_srv1_biom(i))); }// log-likelihood for survey biomass
  for (i=1; i<=nyrs_srv2; i++)  {   ii=yrs_srv2(i);
  surv_like(2) += square((log(obs_srv2_biom(i)+0.0001)-log(pred_srv2(ii)+0.0001) ))/ (2.*square(obs_srv2_se(i)/obs_srv2_biom(i))); }
  for (i=1; i<=nyrs_srv3; i++)  {   ii=yrs_srv3(i);
   surv_like(3) += square((log(obs_srv3_biom(i)+0.0001)-log(pred_srv3(ii)+0.0001) ))/ (2.*square(obs_srv3_se(i)/obs_srv3_biom(i))); }
  for (i=1; i<=nyrs_srv4; i++)  {   ii=yrs_srv4(i);
  surv_like(4) += square((log(obs_srv4_biom(i)+0.0001)-log(pred_srv4(ii)+0.0001) ))/ (2.*square(obs_srv4_se(i)/obs_srv4_biom(i))); }
  for (i=1; i<=nyrs_srv5; i++)  {   ii=yrs_srv5(i);
  surv_like(5) += square((log(obs_srv5_biom(i)+0.0001)-log(pred_srv5(ii)+0.0001) ))/ (2.*square(obs_srv5_se(i)/obs_srv5_biom(i))); }
  for (i=1; i<=nyrs_srv6; i++)  {   ii=yrs_srv6(i);
  surv_like(6) += square((log(obs_srv6_biom(i)+0.0001)-log(pred_srv6(ii)+0.0001) ))/ (2.*square(obs_srv6_se(i)/obs_srv6_biom(i))); }
  for (i=1; i<=nyrs_srv7; i++)  {   ii=yrs_srv7(i);
  surv_like(7) += square((log(obs_srv7_biom(i)+0.0001)-log(pred_srv7(ii)+0.0001) ))/ (2.*square(obs_srv7_se(i)/obs_srv7_biom(i))); }
   
 
   // likelihood for survey biomass 
  surv_like(1) *= wt_srv1 ;  
  surv_like(2) *= wt_srv2 ;  
  surv_like(3) *= wt_srv3 ;  
  surv_like(4) *= wt_srv4 ;  
  surv_like(5) *= wt_srv5 ;  
  surv_like(6) *= wt_srv6 ;  
  surv_like(7) *= wt_srv7 ;  
  surv_like(8) *= wt_srv8 ;  

FUNCTION Multinomial_Likelihood
// Calculate multinomial likelihoods for survey age, fishery size, and survey size and subtract "offset"
  age_like.initialize();

  for (i=1; i <= nyrs_fish1_age; i++)
    age_like(1) -= nsamples_fish1_age(i)*((oac_fish1(i) + 0.001) * log(eac_fish1(i) + 0.001)) ;

  for (i=1; i <= nyrs_srv1_age; i++)
    age_like(2) -= nsamples_srv1_age(i)*((oac_srv1(i) + 0.001) * log(eac_srv1(i) + 0.001)) ;

  for (i=1; i <= nyrs_fish1_size; i++) age_like(3) -= nsamples_fish1_size(i)*((osc_fish1_f(i) + 0.001) * log(esc_fish1_f(i) + 0.001)) ;
  for (i=1; i <= nyrs_fish1_size; i++) age_like(4) -= nsamples_fish1_size(i)*((osc_fish1_m(i) + 0.001) * log(esc_fish1_m(i) + 0.001)) ;
  for (i=1; i <= nyrs_fish2_size; i++) age_like(5) -= nsamples_fish2_size(i)*((osc_fish2(i) + 0.001) * log(esc_fish2(i) + 0.001)) ;
  for (i=1; i <= nyrs_fish3_size; i++) age_like(6) -= nsamples_fish3_size(i)*((osc_fish3_f(i) + 0.001) * log(esc_fish3_f(i) + 0.001)) ;
  for (i=1; i <= nyrs_fish3_size; i++) age_like(7) -= nsamples_fish3_size(i)*((osc_fish3_m(i) + 0.001) * log(esc_fish3_m(i) + 0.001)) ;
  for (i=1; i <= nyrs_fish4_size; i++) age_like(8) -= nsamples_fish4_size(i)*((osc_fish4(i) + 0.001) * log(esc_fish4(i) + 0.001)) ;

  for (i=1; i <= nyrs_srv1_size; i++)  age_like(9) -= nsamples_srv1_size(i)*((osc_srv1_f(i) + 0.001) * log(esc_srv1_f(i) + 0.001)) ;
  for (i=1; i <= nyrs_srv1_size; i++)  age_like(10) -= nsamples_srv1_size(i)*((osc_srv1_m(i) + 0.001) * log(esc_srv1_m(i) + 0.001)) ;
  for (i=1; i <= nyrs_srv2_size; i++)  age_like(11) -= nsamples_srv2_size(i)*((osc_srv2_f(i) + 0.001) * log(esc_srv2_f(i) + 0.001)) ;
  for (i=1; i <= nyrs_srv2_size; i++)  age_like(12) -= nsamples_srv2_size(i)*((osc_srv2_m(i) + 0.001) * log(esc_srv2_m(i) + 0.001)) ;
  for (i=1; i <= nyrs_srv7_size; i++)  age_like(13) -= nsamples_srv7_size(i)*((osc_srv7_f(i) + 0.001) * log(esc_srv7_f(i) + 0.001)) ;
  for (i=1; i <= nyrs_srv7_size; i++)  age_like(14) -= nsamples_srv7_size(i)*((osc_srv7_m(i) + 0.001) * log(esc_srv7_m(i) + 0.001)) ;
  for (i=1; i <= nyrs_srv7_age; i++)  age_like(15) -= nsamples_srv7_age(i)*((oac_srv7(i) + 0.001) * log(eac_srv7(i) + 0.001)) ;
  for (i=1; i <= nyrs_srv2_age; i++)   age_like(16) -= nsamples_srv2_age(i)*((oac_srv2(i) + 0.001) * log(eac_srv2(i) + 0.001)) ;
  
  age_like   -= offset;                       // Subract offsets

 if(data_reweight_switch==0) // no reweighting being done
  {
   age_like(1) *= wt_fish1_age;    //1               // Multiple each likelihood by their weights from .ctl file
   age_like(2) *= wt_srv1_age;  //1
   age_like(3) *= wt_fish1_size;  //1
   age_like(4) *= wt_fish1_size;  //1
   age_like(5) *= wt_fish2_size; //
   age_like(6) *= wt_fish3_size;  //1
   age_like(7) *= wt_fish3_size;  //1
   age_like(8) *= wt_fish4_size; //0
   age_like(9) *= wt_srv1_size;  //1
   age_like(10) *= wt_srv1_size;  //1
   age_like(11) *= wt_srv2_size;  //1
   age_like(12) *= wt_srv2_size;  //1
   age_like(13) *= wt_srv7_size;  //1
   age_like(14) *= wt_srv7_size;  //1
   age_like(15) *= wt_srv7_age;  //1
   age_like(16) *= wt_srv1_age;  //1
 }

 if(data_reweight_switch==1) // data reweighting being done
  {
   age_like(1) *= wt_fish1_age_iter;    //1               // Multiple each likelihood by their weights from .ctl file
   age_like(2) *= wt_srv1_age_iter;  //1
   age_like(3) *= wt_fish1_size_female_iter;  //1
   age_like(4) *= wt_fish1_size_male_iter;  //1
   age_like(5) *= wt_fish2_size; //
   age_like(6) *= wt_fish3_size_female_iter;  //1
   age_like(7) *= wt_fish3_size_male_iter;  //1
   age_like(8) *= wt_fish4_size; //0
   age_like(9) *= wt_srv1_size_female_iter;  //1
   age_like(10) *= wt_srv1_size_male_iter;  //1
   age_like(11) *= wt_srv2_size_female_iter;  //1
   age_like(12) *= wt_srv2_size_male_iter;  //1
   age_like(13) *= wt_srv7_size_female_iter;  //1
   age_like(14) *= wt_srv7_size_male_iter;  //1
   age_like(15) *= wt_srv7_age;  //1
   age_like(16) *= wt_srv2_age_iter;  //1
 }
FUNCTION Sel_Like

// #####################################################################################################################
// ########## THESE PENALTIES ARE NOT INCLUDED IN OBJ FXN APPARENTLY
// ########################################################################################################################


// Calculate penalty function for selectivity
  sel_like.initialize();
  if (active(log_srv1_sel_coffs_f) ) {

// Differences in selectivity between adjacent ages
  sel_like(1)   +=wt_sel_reg_fish1 * norm2(first_difference(first_difference(log_fish1_sel_f)));  // Constrains selectivities to be smooth
  sel_like(1)   +=wt_sel_reg_fish1 * norm2(first_difference(first_difference(log_fish1_sel_m)));  // Constrains selectivities to be smooth
  sel_like(2)   +=wt_sel_reg_fish2 * norm2(first_difference(first_difference(log_fish2_sel)));  // Constrains selectivities to be smooth
  sel_like(3)   +=wt_sel_reg_fish3 * norm2(first_difference(first_difference(log_fish3_sel_f)));  // Constrains selectivities to be smooth
  sel_like(4)   +=wt_sel_reg_fish4 * norm2(first_difference(first_difference(log_fish4_sel_f)));  // Constrains selectivities to be smooth
  sel_like(4)   +=wt_sel_reg_fish4 * norm2(first_difference(first_difference(log_fish4_sel_m)));  // Constrains selectivities to be smooth

  sel_like(5) +=wt_sel_reg_srv1 * norm2(first_difference(first_difference(log_srv1_sel_f)));  // Constrains selectivities to be smooth
  sel_like(6) +=wt_sel_reg_srv2 * norm2(first_difference(first_difference(log_srv2_sel_f)));  // Constrains selectivities to be smooth
  sel_like(5) +=wt_sel_reg_srv1 * norm2(first_difference(first_difference(log_srv1_sel_m)));  // Constrains selectivities to be smooth
  sel_like(6) +=wt_sel_reg_srv2 * norm2(first_difference(first_difference(log_srv2_sel_m)));  // Constrains selectivities to be smooth
  sel_like(7) +=wt_sel_reg_srv2 * norm2(first_difference(first_difference(log_srv7_sel_f)));  // Constrains selectivities to be smooth
  sel_like(8) +=wt_sel_reg_srv2 * norm2(first_difference(first_difference(log_srv7_sel_m)));  // Constrains selectivities to be smooth

// Differences in selectivity between adjacent ages when selectivity for first age is greater
//   Affects the degree of dome-shape

  for (j=1;j<nages;j++)
    if (log_fish1_sel_f(j)>log_fish1_sel_f(j+1))
      sel_like(9) += wt_sel_dome_fish1 *square(log_fish1_sel_f(j)-log_fish1_sel_f(j+1));  //Prevents dome-shapedness
    if (log_fish1_sel_m(j)>log_fish1_sel_m(j+1))
      sel_like(9) += wt_sel_dome_fish1 *square(log_fish1_sel_m(j)-log_fish1_sel_m(j+1));  //Prevents dome-shapedness
  for (j=1;j<nages;j++)
    if (log_fish2_sel(j)>log_fish2_sel(j+1))
      sel_like(10) += wt_sel_dome_fish2 *square(log_fish2_sel(j)-log_fish2_sel(j+1));  //Prevents dome-shapedness
  for (j=1;j<nages;j++)
    if (log_fish3_sel_f(j)>log_fish3_sel_f(j+1))
      sel_like(11) += wt_sel_dome_fish3 *square(log_fish3_sel_f(j)-log_fish3_sel_f(j+1));  //Prevents dome-shapedness
//  for (j=1;j<nages;j++)
//    if (log_fish4_sel(j)>log_fish4_sel(j+1))
//      sel_like(12) += wt_sel_dome_fish4 *square(log_fish4_sel(j)-log_fish4_sel(j+1));  //Prevents dome-shapedness


    for (j=1;j<nages;j++)
      if (log_srv1_sel_f(j)>log_srv1_sel_f(j+1))
        sel_like(13) +=wt_sel_dome_srv1 *square(log_srv1_sel_f(j)-log_srv1_sel_f(j+1));
      if (log_srv1_sel_m(j)>log_srv1_sel_m(j+1))
        sel_like(13) +=wt_sel_dome_srv1 *square(log_srv1_sel_m(j)-log_srv1_sel_m(j+1));
    for (j=1;j<nages;j++)
      if (log_srv2_sel_f(j)>log_srv2_sel_f(j+1))
        sel_like(14) +=wt_sel_dome_srv2 *square(log_srv2_sel_f(j)-log_srv2_sel_f(j+1));
      if (log_srv2_sel_m(j)>log_srv2_sel_m(j+1))
        sel_like(14) +=wt_sel_dome_srv2 *square(log_srv2_sel_m(j)-log_srv2_sel_m(j+1));
    for (j=1;j<nages;j++) {
      if (log_srv7_sel_f(j)>log_srv7_sel_f(j+1))
        sel_like(15) +=wt_sel_dome_srv2 *square(log_srv7_sel_f(j)-log_srv7_sel_f(j+1));
      if (log_srv7_sel_m(j)>log_srv7_sel_m(j+1))
        sel_like(15) +=wt_sel_dome_srv2 *square(log_srv7_sel_m(j)-log_srv7_sel_m(j+1)); }
//    for (j=1;j<nages;j++)
//      if (log_srv4_sel(j)>log_srv4_sel(j+1))
 //       sel_like(16) +=wt_sel_dome_srv4 *square(log_srv4_sel(j)-log_srv4_sel(j+1));
 }
 
FUNCTION double round(double r)   // what does this do???????????????????????????????????????

    return double((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5));

FUNCTION Get_Population_Projection

 //  Abundance at start of first projection year
 // stdev of recvar
 
   int k;
   if(mceval_phase())  // use MCMC outputs to determine proj recruitment
    {
     // random_number_generator r(1000);
     stdev_rec              = sqrt(norm2(value(log_rec_dev(1979,endyr-recage))-mean(value(log_rec_dev(1979,endyr-recage))))/(size_count(value(log_rec_dev(1979,endyr-recage)))-1));
     
     k=round(value(stdev_rec)*10000);

     N_proj_f(endyr+1,1)    = mfexp(value(log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)))/2;
     N_proj_m(endyr+1,1)    = mfexp(value(log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)))/2;
    }
   else // use recent average recruitment for projections
    {
     N_proj_f(endyr+1,1)    = mfexp(value(log(mean(pred_rec(1979,endyr-recage)))))/2;
     N_proj_m(endyr+1,1)    = mfexp(value(log(mean(pred_rec(1979,endyr-recage)))))/2;
    }
   for (j=1; j<nages-1;j++)
    {
     k=k+j;
     N_proj_f(endyr+1,j+1)  = natage_f(endyr,j)*value(S_f(endyr,j));
     N_proj_m(endyr+1,j+1)  = natage_m(endyr,j)*value(S_m(endyr,j));
    }
    
   N_proj_f(endyr+1,nages)  = value(natage_f(endyr,nages-1))*value(S_f(endyr,nages-1))+ value(natage_f(endyr,nages))*value(S_f(endyr,nages));
   N_proj_m(endyr+1,nages)  = value(natage_m(endyr,nages-1))*value(S_m(endyr,nages-1))+ value(natage_m(endyr,nages))*value(S_m(endyr,nages));

   spawn_biom_proj(endyr+1) = elem_prod(N_proj_f(endyr+1),pow(mfexp(-yieldratio*FABC_tot_proj_f-value(natmort)),spawn_fract)) * weight_maturity_prod_f(endyr);
   tot_biom_proj(endyr+1)   = N_proj_f(endyr+1)*weight_f(endyr)+N_proj_m(endyr+1)*weight_m(endyr);
   
  for (i=endyr+1;i<=endyr+15;i++)
   {
    //  F ABC 
    if (spawn_biom_proj(i)/B40 > 1.)
     {
      FABC_proj             = value(F40);
      FOFL_proj             = value(F35);
     }
    else
     {
      FABC_proj             = value(F40) * (spawn_biom_proj(i)/value(B40) - 0.05)/(1 - 0.05); 
      FOFL_proj             = value(F35) * (spawn_biom_proj(i)/value(B40) - 0.05)/(1 - 0.05);
     }

    for (j=1;j<=nages;j++)
     {
      if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
       {
        FABC_tot_proj_f(j)  = fish4_sel_f(j)* FABC_proj * fratio + fish3_sel_f(j)* FABC_proj * (1-fratio);
        FABC_tot_proj_m(j)  = fish4_sel_m(j)* FABC_proj * fratio + fish3_sel_m(j)* FABC_proj * (1-fratio);
        FOFL_tot_proj_f(j)  = fish4_sel_f(j)* FOFL_proj * fratio + fish3_sel_f(j)* FOFL_proj * (1-fratio);
        FOFL_tot_proj_m(j)  = fish4_sel_m(j)* FOFL_proj * fratio + fish3_sel_m(j)* FOFL_proj * (1-fratio);
       }
      if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
       {
        FABC_tot_proj_f(j)  = fish5_sel_f(j)* FABC_proj * fratio + fish3_sel_f(j)* FABC_proj * (1-fratio);
        FABC_tot_proj_m(j)  = fish5_sel_m(j)* FABC_proj * fratio + fish3_sel_m(j)* FABC_proj * (1-fratio);
        FOFL_tot_proj_f(j)  = fish5_sel_f(j)* FOFL_proj * fratio + fish3_sel_f(j)* FOFL_proj * (1-fratio);
        FOFL_tot_proj_m(j)  = fish5_sel_m(j)* FOFL_proj * fratio + fish3_sel_m(j)* FOFL_proj * (1-fratio);
       }
      if(ph_ifq<1) // no post-IFQ sel est
       {
        FABC_tot_proj_f(j)  = fish1_sel_f(j)* FABC_proj * fratio + fish3_sel_f(j)* FABC_proj * (1-fratio);
        FABC_tot_proj_m(j)  = fish1_sel_m(j)* FABC_proj * fratio + fish3_sel_m(j)* FABC_proj * (1-fratio);
        FOFL_tot_proj_f(j)  = fish1_sel_f(j)* FOFL_proj * fratio + fish3_sel_f(j)* FOFL_proj * (1-fratio);
        FOFL_tot_proj_m(j)  = fish1_sel_m(j)* FOFL_proj * fratio + fish3_sel_m(j)* FOFL_proj * (1-fratio);
       }

        Z_proj_f(j)         = FABC_tot_proj_f(j)+ natmort;
        Z_proj_m(j)         = FABC_tot_proj_m(j)+ natmort;
        ZOFL_proj_f(j)      = FOFL_tot_proj_f(j)+ value(natmort);
        ZOFL_proj_m(j)      = FOFL_tot_proj_m(j)+ value(natmort);
        S_proj_f(j)         = mfexp(-1.0* Z_proj_f(j));
        S_proj_m(j)         = mfexp(-1.0* Z_proj_m(j));
     }

   //  Catch 
    for (j=1;j<=nages;j++)
     { 
      catage_proj_f(i,j)     = yieldratio*N_proj_f(i,j)* FABC_tot_proj_f(j)/Z_proj_f(j)*(1.-mfexp(-Z_proj_f(j)));
      catage_proj_m(i,j)     = yieldratio*N_proj_m(i,j)* FABC_tot_proj_m(j)/Z_proj_m(j)*(1.-mfexp(-Z_proj_m(j)));
      catage_proj_OFL_f(i,j) = yieldratio*N_proj_f(i,j)* FOFL_tot_proj_f(j)/ZOFL_proj_f(j)*(1.-mfexp(-ZOFL_proj_f(j)));
      catage_proj_OFL_m(i,j) = yieldratio*N_proj_m(i,j)* FOFL_tot_proj_m(j)/ZOFL_proj_m(j)*(1.-mfexp(-ZOFL_proj_m(j)));
     }
     
      pred_catch_proj(i)     = (catage_proj_f(i)*weight_f(endyr)+catage_proj_m(i)*weight_m(endyr))/yieldratio;
      pred_catch_proj_OFL(i) = (catage_proj_OFL_f(i)*weight_f(endyr)+catage_proj_OFL_m(i)*weight_m(endyr))/yieldratio;

//  Next year's abundance
    if (i < endyr+15)
    {
     if(mceval_phase())
      {
       stdev_rec             = sqrt(norm2(value(log_rec_dev(1979,endyr-recage))-mean(value(log_rec_dev(1979,endyr-recage))))/(size_count(value(log_rec_dev(1979,endyr-recage)))-1));
       
       k=round(value(spawn_biom(endyr)*10000))+i;
       k=k+i;
       
       N_proj_f(i+1,1)       = mfexp(value(log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)))/2;
       N_proj_m(i+1,1)       = mfexp(value(log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)))/2;
      }
     else
      {
       N_proj_f(i+1,1)       = mfexp(value(log(mean(pred_rec(1979,endyr-recage)))))/2;
       N_proj_m(i+1,1)       = mfexp(value(log(mean(pred_rec(1979,endyr-recage)))))/2;
      }

     for (j=1; j<nages-1;j++)
      {
       N_proj_f(i+1,j+1)     = N_proj_f(i,j)  * mfexp(-yieldratio*FABC_tot_proj_f(j)-value(natmort));;
       N_proj_m(i+1,j+1)     = N_proj_m(i,j)  * mfexp(-yieldratio*FABC_tot_proj_m(j)-value(natmort));
      }
      
       N_proj_f(i+1,nages)   = N_proj_f(i,nages-1)* mfexp(-yieldratio*FABC_tot_proj_f(nages-1)-value(natmort))+ N_proj_f(i,nages)   * mfexp(-yieldratio*FABC_tot_proj_f(nages)-value(natmort));
       N_proj_m(i+1,nages)   = N_proj_m(i,nages-1)* mfexp(-yieldratio*FABC_tot_proj_m(nages-1)-value(natmort))+ N_proj_m(i,nages)   * mfexp(-yieldratio*FABC_tot_proj_m(nages)-value(natmort));

       spawn_biom_proj(i+1)  = elem_prod(N_proj_f(i+1),pow(mfexp(-yieldratio*FABC_tot_proj_f-value(natmort)),spawn_fract)) * weight_maturity_prod_f(endyr);  // Right way
       tot_biom_proj(i+1)    = N_proj_f(i+1)*weight_f(endyr)+N_proj_m(i+1)*weight_m(endyr);

    }  // end if statement for not in terminal year
  }    // end year for statement

     if (spawn_biom_proj(endyr+1)/B40 > 1.)
      {
       FABC    = value(F40);
       FOFL    = value(F35); 
       FABC2   = value(F40);
       FOFL2   = value(F35);
      }
     else
      {
       FABC    = value(F40) * (spawn_biom_proj(endyr+1)/value(B40) - 0.05)/(1 - 0.05); 
       FOFL    = value(F35) * (spawn_biom_proj(endyr+1)/value(B40) - 0.05)/(1 - 0.05);  
       FABC2   = value(F40) * (spawn_biom_proj(endyr+2)/value(B40) - 0.05)/(1 - 0.05); 
       FOFL2   = value(F35) * (spawn_biom_proj(endyr+2)/value(B40) - 0.05)/(1 - 0.05);
      }
      
       OFL=pred_catch_proj_OFL(endyr+1);
       ABC=pred_catch_proj(endyr+1);

FUNCTION Evaluate_Objective_Function 
  obj_fun.initialize();
  ssqcatch.initialize();
  rec_like.initialize();
  F_mort_regularity.initialize();
  M_mort_regularity.initialize();
  avg_sel_penalty.initialize();
  rec_like_bias_adj.initialize();
  Surv_Likelihood();                                  // Likelihood function for survey biomass
  ssqcatch  +=  wt_ssqcatch_fish1 *norm2(log(obs_catch_fish1+0.01)-log(pred_catch_fish1+0.01));
  ssqcatch  +=  wt_ssqcatch_fish3 *norm2(log(obs_catch_fish3+0.8)-log(pred_catch_fish3+0.8));



    for(y=(styr-nages+2);y<=endyr_rec_est;y++)   /// implement Methot and Taylor (2011) bias ramp adjustment for recruit penalty term
    {
     if(sigma_R_early_switch==1 && y<sigma_R_early_end)
      {
       rec_like_bias_adj(y)= square(log_rec_dev(y)/sigma_R_early) + b(y)*log(sigma_R_early);
      }
     else
      {
       rec_like_bias_adj(y)= square(log_rec_dev(y)/sigr) + b(y)*log(sigr);
      }
    }


  switch (SrType)
  {
    case 3:
    {
      if (rec_like_type==2)
      if(bias_ramp==1 || sigma_R_early_switch==1)
            rec_like      = wt_rec_var * 0.5*sum(rec_like_bias_adj);
       else
            rec_like      = wt_rec_var * 0.5*(norm2(log_rec_dev/sigr) + log(sigr));
     //     rec_like      = wt_rec_var*(norm2(log_rec_dev-sigr*sigr/2.)/(2.*square(sigr)) + (size_count(log_rec_dev))*log(sigr));
      else
        rec_like = wt_rec_var*(norm2(log_rec_dev));
        break;
    }
    default: 
    {
        dvar_vector stmp(styr_rec,endyr);
        for (i=styr_rec;i<=endyr;i++)
        stmp(i) = Sp_Biom(i-recage);
        srm_rec   = SRecruit(stmp);
        dvar_vector   chi(styr_rec_est,endyr_rec_est);
        chi         = log(elem_div(sam_rec(styr_rec_est,endyr_rec_est) , srm_rec(styr_rec_est,endyr_rec_est)));
        dvariable SSQRec = norm2( chi + sigrsq/2.) ;
        rec_like    = .5*SSQRec/sigrsq + nrecs_est*log(sigr); 
        rec_like   += .5*norm2( log_rec_dev(styr_rec,styr_rec_est) )/sigrsq + (styr_rec_est-styr_rec)*log(sigr) ; 
      if (endyr>endyr_rec_est)
        rec_like += .5*norm2( log_rec_dev(endyr_rec_est,endyr  ) )/sigrsq + (endyr-endyr_rec_est)  *log(sigr) ; 
      break;
    }
  }

      F_mort_regularity  = wt_fmort_reg * norm2(log_F_devs_fish1);// Penalty function for fishing mortality deviations
      F_mort_regularity += wt_fmort_reg * norm2(log_F_devs_fish3);// Penalty function for fishing mortality deviations

  if(active(log_F_devs_fish1))                          // Penalty function for fishing mortality deviations
    obj_fun         += F_mort_regularity;
    obj_fun         += sum(priors);               //Add priors
    
  if (current_phase()==ph_Fdev)
    obj_fun         += 10*(norm2(log_F_devs_fish1)+norm2(log_F_devs_fish3));         //(was-0.3) Penalty early on to scale population...
                       //1000*(square(log(mean(Fmort_fish1)/0.1))+square(log(mean(Fmort_fish3)/0.1)));         //Penalty early on to scale population...

  if (active(log_mF50) && last_phase()) 
    obj_fun         += sprpen;                    // To solve for the F40 etc.
    


  if(active(log_M_devs))                          // Penalty function for yearly natural mortality deviations
   {
     M_mort_regularity  = wt_M_reg * norm2(log_M_devs);
     obj_fun           += M_mort_regularity;
   }
  if(active(log_M_devs_age))                          // Penalty function for yearly natural mortality deviations
   {
     M_mort_regularity  = wt_M_reg * norm2(log_M_devs_age);
     obj_fun           += M_mort_regularity;
   }

      Multinomial_Likelihood();                           // Multinomial likelihood

// Sum objective function
  obj_fun           += ssqcatch ;
  obj_fun           += sum(surv_like);
  obj_fun           += sum(age_like);
  Like = obj_fun;                     // Put here to capture the data likelihood
  obj_fun           += rec_like;
  


    cout<<"monitoring SSB "<<sum(elem_prod(natage_f(endyr),weight_maturity_prod_f(endyr)))<<endl;
 
REPORT_SECTION
// Beginning of all outputting
// Goes to routine that automatically creates input file for projection model
  if (last_phase()) {
    write_projout();
    write_newproj();
    write_sarareport();
    write_HQreport();
    }
// Output file (tem.rep) which is loaded into tem.xls to display output
  
  report<<"****Executive mary Material*****"<<endl;
  report<<"     Model name"     <<endl;
  report<<model_name<<endl;
  report<<"     .dat file"     <<endl;
  report<<data_file<<endl;
  report<<"     Number parameters estimated"     <<endl;
  report<<initial_params::nvarcalc()<<endl;
  report<<"     TotalBiomass for "<<endyr+1<<endl;
  report<<tot_biom_proj(endyr+1)<<endl;
  report<<"     TotalBiomass for "<<endyr+2     <<endl;
  report<<tot_biom_proj(endyr+2)<<endl;
  report<<"     Female_Spawning Biomass for "<<endyr+1     <<endl;
  report<<spawn_biom_proj(endyr+1)<<endl;
  report<<"     Female_Spawning_Biomass for "<<endyr+2     <<endl;
  report<<spawn_biom_proj(endyr+2)<<endl;
  report<<"     B_100"     <<endl;
  report<<0.5*SB0*mean(pred_rec(1979,endyr-recage))<<endl;
  report<<"     B_40"     <<endl;
  report<<B40<<endl;
  report<<"     B_35"     <<endl;
  report<<0.5*SBF35*mean(pred_rec(1979,endyr-recage))<<endl;
  report<<"     Mean_Recruitment"     <<endl;
  report<<mean(pred_rec(1979,endyr-recage))<<endl;
  report<<"     F_40"     <<endl;
  report<<F40<<endl;
  report<<"     F_35"     <<endl;
  report<<F35<<endl;
  report<<"     F_ABC for "<<endyr+1     <<endl;
  report<<FABC<<endl;
  report<<"     F_ABC for "<<endyr+2     <<endl;
  report<<FABC2<<endl;
  report<<"     ABC for "<<endyr+1     <<endl;
  report<<pred_catch_proj(endyr+1)<<endl;
  report<<"     ABC for "<<endyr+2     <<endl;
  report<<pred_catch_proj(endyr+2)<<endl;
  report<<"     F_OFL for "<<endyr+1     <<endl;
  report<<FOFL<<endl;
  report<<"     F_OFL for "<<endyr+2     <<endl;
  report<<FOFL2<<endl;
  report<<"     OFL for "<<endyr+1     <<endl;
  report<<OFL<<endl; 
  report<<"     OFL for "<<endyr+2     <<endl;
  report<<pred_catch_proj_OFL(endyr+2)<<endl; 
  report<<"     Total likelihood"     <<endl;
  report<<obj_fun<<endl;
  report<<"     Data likelihood"     <<endl;
  report<<Like<<endl<<endl;


  report << "Iterative_Weights_Likelihoods  " << endl;

  report << wt_fish1_age_iter <<" "<<age_like(1)  <<" " ; report << "Fishery_Age_Composition_Likelihood" << endl;
  report << wt_srv1_age_iter <<" "<<age_like(2)  <<" " ; report << "Survey_Age_Composition_Likelihood_DomesticLL" << endl;
  report << wt_srv2_age_iter <<" "<<age_like(16)  <<" " ; report << "Survey_Age_Composition_Likelihood_CooperativeLL" << endl;
  report << wt_fish1_size_female_iter<<" "<<age_like(3)  <<" " ; report << "Fishery_Size_Composition_Likelihood_Fixed_Female" << endl;
  report << wt_fish1_size_male_iter<<" "<<age_like(4)  <<" " ; report << "Fishery_Size_Composition_Likelihood_Fixed_Male" << endl;
  report << wt_fish3_size_female_iter<<" "<<age_like(6) <<" " ; report << "Fishery_Size_Composition_Likelihood_Trawl_Female" << endl;
  report << wt_fish3_size_male_iter<<" "<<age_like(7)  <<" " ; report << "Fishery_Size_Composition_Likelihood_Trawl_Male" << endl;
  report << wt_srv1_size_female_iter<<" "<<age_like(9) <<" " ; report << "Survey_Size_Composition_Likelihood_Domestic LL_Female" << endl;
  report << wt_srv1_size_male_iter<<" "<<age_like(10)  <<" " ; report << "Survey_Size_Composition_Likelihood_Domestic LL_Male" << endl;
  report << wt_srv2_size_female_iter<<" "<<age_like(11)  <<" " ; report << "Survey_Size_Composition_Likelihood_Cooperative_LL_Female" << endl;
  report << wt_srv2_size_male_iter<<" "<<age_like(12)  <<" " ; report << "Survey_Size_Composition_Likelihood_Cooperative_LL_Male" << endl;
  report << wt_srv7_size_female_iter<<" "<<age_like(13)<<" " ; report << "Survey_Size_Composition_Likelihood_GOATrawl_Female" << endl;
  report << wt_srv7_size_male_iter<<" "<<age_like(14)  <<" " ; report << "Survey_Size_Composition_Likelihood_GOATrawl_Male" << endl;


  report<<" ************   Some more parameter estimates and their SDs ************"<<endl;
 
  if(last_phase()) {
    // add standard deviation data types    
  report<<"   q_domestic   "<<endl;
  report<<q_srv1<<" "<<q_srv1.sd<<endl;
  report<<"   q_cooperative  "<<endl;
   report<<q_srv2<<" "<<q_srv2.sd<<endl;
  report<<"   q_trawl  "<<endl;
   report<<q_srv7<<" "<<q_srv7.sd<<endl;
 /*
  report<<natmort<<" "<<nattymort.sd<<endl;
  report<<"  sigr   "<<endl;  
  report<<sigr<<" "<<cigar.sd<<endl;  
  report<<"   log_mean_rec"<<endl;
  report<<log_mean_rec<<" "<<LMR.sd<<endl;
 */
  report<<"   F_40"<<endl;
  report<<F40<<" "<<F40.sd<<endl;
  report<<"    tot_biom"<<endl;
  report<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  report<<"   spawn_biom"<<endl;
  report<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  report<<"    B40"<<endl;
  report<<B40<<" "<<B40.sd<<endl;
  report<<"   ABC"<<endl;
  report<<pred_catch_proj(endyr+1)<<" "<<pred_catch_proj.sd(endyr+1)<<endl<<endl;
 
 }

  report << "Parameter Phases"<<endl;
  report << "ph_mean_rec ph_recdev ph_Rzero ph_steepness ph_sigr" <<endl;
  report <<  ph_mean_rec << " "<< ph_recdev <<" "<<  ph_Rzero<<" "<<ph_steepness<<" "<<ph_sigr  <<endl;
  report << "ph_avg_F ph_Fdev ph_F50 "<<endl;
  report <<  ph_avg_F << " "<< ph_Fdev <<" "<<  ph_F50 <<endl;
  report << "ph_m ph_Mdevs ph_Mdevs_age ph_mdelta"<<endl;
  report <<  ph_m << " "<< ph_Mdevs <<" "<<  ph_Mdevs_age<<" "<<ph_mdelta  <<endl;
  report << "ph_fish_sel ph_fish2_sel ph_ifq ph_ifq_block2 ph_fish_sel_delt ph_fish_sel_delt_alt"<<endl;
  report <<  ph_fish_sel << " "<< ph_fish2_sel << " "<< ph_ifq <<" "<<  ph_ifq_block2 <<" "<<ph_fish_sel_delt<<" "<<ph_fish_sel_delt_alt<< endl;
  report << "ph_srv1_sel ph_srv2_sel ph_LL_block2 ph_srv_sel_delt ph_srv_sel_delt_alt"<<endl;
  report <<  ph_srv1_sel << " "<< ph_srv2_sel<<" "<<ph_LL_block2<<" "<< ph_srv_sel_delt<<" "<<ph_srv_sel_delt_alt <<endl;
  report << "ph_q_srv1 ph_q_srv2 ph_q_srv3 ph_q_srv4 ph_q_srv5" <<endl;
  report <<  ph_q_srv1 << " "<< ph_q_srv2 <<" "<<  ph_q_srv3<<" "<<ph_q_srv4<<" "<<ph_q_srv5  <<endl;
  report << "ph_q_srv6 ph_q_srv7 ph_q_srv8 ph_srv2_q2 " <<endl;
  report <<  ph_q_srv6 << " "<< ph_q_srv7 <<" "<<  ph_q_srv8<<" "<<ph_srv2_q2  <<endl<<endl<<endl;

  report<<model_name<<endl;
  report<<data_file<<endl;
  report<<"Num_parameters_Estimated "<<initial_params::nvarcalc()<<endl;
  report << "Year "<< yy <<endl;
  report << "Pred_Catch Fixed Gear "<< pred_catch_fish1<<endl<<"Pred catch trawl "<<pred_catch_fish3 <<endl;
  report << "Obs_Catch Fixed Gear "<< obs_catch_fish1<<endl<<" Obs_Catch Trawl "<<obs_catch_fish3 <<endl;

  report << "Survival Female"<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<S_f(i) <<endl; report<<endl;
  report << "Survival Male"<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<S_m(i) <<endl; report<<endl;

  report << "Natural Mortality"<<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<M(i) <<endl; report<<endl;
  
  report << "Numbers Females "<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<natage_f(i) <<endl; report<<endl;
  report << "Numbers Males"<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<natage_m(i) <<endl; report<<endl;

  report << "Numbers at Length Females "<<len_bin_labels <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<num_len_f(i) <<endl; report<<endl;
  report << "Numbers at Length  Males"<<len_bin_labels <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<num_len_m(i) <<endl; report<<endl;

  report << "Age "<<aa <<endl;
  report << "Fishery_sel1 Females"<<fish1_sel_f  <<endl;
  report << "Fishery_sel1 Males"<<fish1_sel_m  <<endl;
  report << "Fishery_sel2 "<<fish2_sel  <<endl;
  report << "Fishery_sel3_f "<<fish3_sel_f  <<endl;
  report << "Fishery_sel3_m "<<fish3_sel_m  <<endl;
  report << "Fishery_sel4_f "<<fish4_sel_f  <<endl;
  report << "Fishery_sel4_m "<<fish4_sel_m  <<endl;
  report << "Fishery_sel5_f "<<fish5_sel_f  <<endl;
  report << "Fishery_sel5_m "<<fish5_sel_m  <<endl;
  report << "Survey_sel1 Female"<<srv1_sel_f  <<endl<<endl;
  report << "Survey_sel1 male"<<srv1_sel_m  <<endl<<endl;
  report << "Survey_sel2 Female"<<srv2_sel_f  <<endl<<endl;
  report << "Survey_sel2 male"<<srv2_sel_m  <<endl<<endl;
  report << "Survey_sel7 Female"<<srv7_sel_f  <<endl<<endl;
  report << "Survey_sel7 male"<<srv7_sel_m  <<endl<<endl;
  report << "Survey_sel10 Female"<<srv10_sel_f  <<endl<<endl;
  report << "Survey_sel10 male"<<srv10_sel_m  <<endl<<endl;
  
  sdnr_fish1_age = 0;
  sdnr_fish1_size = 0;
  sdnr_fish3_size = 0;
  sdnr_srv1_age = 0;
  sdnr_srv2_age = 0;
  sdnr_srv1_size = 0;
  sdnr_srv2_size = 0;
  sdnr_srv7_size = 0;


  report << "Obs_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) report << yrs_fish1_age(i)<<" "<<oac_fish1(i)<<endl; report<<endl;
  report << "Pred_P_fish1_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) report << yrs_fish1_age(i)<<" "<<eac_fish1(i) <<endl; report<<endl;

  report << "Obs_P_fish1_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++)  report << yrs_fish1_size(i)<<" "<<osc_fish1_f(i)<<endl; report<<endl; 
  report << "Pred_P_fish1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) report << yrs_fish1_size(i)<<" "<<esc_fish1_f(i) <<endl; report<<endl;

    report << "Obs_P_fish1_size Male"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) report << yrs_fish1_size(i)<<" "<<osc_fish1_m(i)<<endl; report<<endl;
  report << "Pred_P_fish1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) report << yrs_fish1_size(i)<<" "<<esc_fish1_m(i) <<endl; report<<endl;

   report << "Obs_P_fish3_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++)  report << yrs_fish3_size(i)<<" "<<osc_fish3_m(i)<<endl; report<<endl; 
  report << "Pred_P_fish3_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) report << yrs_fish3_size(i)<<" "<<esc_fish3_m(i) <<endl; report<<endl;

  report << "Obs_P_fish3_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++)   report << yrs_fish3_size(i)<<" "<<osc_fish3_f(i)<<endl; report<<endl; 
  report << "Pred_P_fish3_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) report << yrs_fish3_size(i)<<" "<<esc_fish3_f(i) <<endl; report<<endl;

  report << "Obs_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++)  report << yrs_srv1_age(i)<<" "<<oac_srv1(i)<<endl; report<<endl; 
  report << "Pred_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; report<<endl;

     report << "Obs_P_srv2_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv2_age;i++)   report << yrs_srv2_age(i)<<" "<<oac_srv2(i)<<endl; report<<endl; 
  report << "Pred_P_srv2_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv2_age;i++) report << yrs_srv2_age(i)<<" "<<eac_srv2(i) <<endl; report<<endl;

  report << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<osc_srv1_f(i)<<endl; report<<endl; 
  report << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<esc_srv1_f(i) <<endl; report<<endl;

    report << "Obs_P_srv1_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++)  report << yrs_srv1_size(i)<<" "<<osc_srv1_m(i)<<endl; report<<endl; 
  report << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<esc_srv1_m(i) <<endl; report<<endl;

  report << "Obs_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++)  report << yrs_srv2_size(i)<<" "<<osc_srv2_f(i)<<endl; report<<endl; 
  report << "Pred_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<esc_srv2_f(i) <<endl; report<<endl;

    report << "Obs_P_srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++)  report << yrs_srv2_size(i)<<" "<<osc_srv2_m(i)<<endl; report<<endl; 
  report << "Pred_P_srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<esc_srv2_m(i) <<endl; report<<endl;
  
    report << "Obs_P_srv7_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++)  report << yrs_srv7_size(i)<<" "<<osc_srv7_f(i)<<endl; report<<endl; 
  report << "Pred_P_srv7_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) report << yrs_srv7_size(i)<<" "<<esc_srv7_f(i) <<endl; report<<endl;

    report << "Obs_P_srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++)  report << yrs_srv7_size(i)<<" "<<osc_srv7_m(i)<<endl; report<<endl;
  report << "Pred_P_srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) report << yrs_srv7_size(i)<<" "<<esc_srv7_m(i) <<endl; report<<endl;


  report << "Obs_full_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) {
      sdnr_fish1_age +=sdnr(eac_fish1(i),oac_fish1(i),wt_fish1_age*double(nsamples_fish1_age(i)))/nyrs_fish1_age;
      report << yrs_fish1_age(i)<<" "<<oac_fish1(i) 
      <<" eff_N "<<(1-eac_fish1(i))*eac_fish1(i)/norm2(oac_fish1(i)-eac_fish1(i))  <<" N "<<nsamples_fish1_age(i)
      <<" SDNR "<< sdnr(eac_fish1(i),oac_fish1(i),wt_fish1_age*double(nsamples_fish1_age(i)))<<endl; report<<endl; }
  report << "Pred_full_P_fish1_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) report << yrs_fish1_age(i)<<" "<<eac_fish1(i) <<endl; report<<endl;

  report << "Obs_full_P_fish1_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++)  {
   sdnr_fish1_size += sdnr(esc_fish1_f(i),osc_fish1_f(i),wt_fish1_size*double(nsamples_fish1_size(i)))/nyrs_fish1_size/2;
  report << yrs_fish1_size(i)<<" "<<osc_fish1_f(i) 
      <<" eff_N "<<(1-esc_fish1_f(i))*esc_fish1_f(i)/norm2(osc_fish1_f(i)-esc_fish1_f(i))  <<" N "<<nsamples_fish1_size(i)
      <<" SDNR "<< sdnr(esc_fish1_f(i),osc_fish1_f(i),wt_fish1_size*double(nsamples_fish1_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_fish1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) report << yrs_fish1_size(i)<<" "<<esc_fish1_f(i) <<endl; report<<endl;

    report << "Obs_full_P_fish1_size Male"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) {
    sdnr_fish1_size += sdnr(esc_fish1_m(i),osc_fish1_m(i),wt_fish1_size*double(nsamples_fish1_size(i)))/nyrs_fish1_size/2;
  report << yrs_fish1_size(i)<<" "<<osc_fish1_m(i) 
      <<" eff_N "<<(1-esc_fish1_m(i))*esc_fish1_m(i)/norm2(osc_fish1_m(i)-esc_fish1_m(i))  <<" N "<<nsamples_fish1_size(i)
      <<" SDNR "<< sdnr(esc_fish1_m(i),osc_fish1_m(i),wt_fish1_size*double(nsamples_fish1_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_fish1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) report << yrs_fish1_size(i)<<" "<<esc_fish1_m(i) <<endl; report<<endl;

   report << "Obs_full_P_fish3_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) {
    sdnr_fish3_size+= sdnr(esc_fish3_m(i),osc_fish3_m(i),wt_fish3_size*double(nsamples_fish3_size(i)))/nyrs_fish3_size/2;
  report << yrs_fish3_size(i)<<" "<<osc_fish3_m(i) 
      <<" eff_N "<<(1-esc_fish3_m(i))*esc_fish3_m(i)/norm2(osc_fish3_m(i)-esc_fish3_m(i))  <<" N "<<nsamples_fish3_size(i)
      <<" SDNR "<< sdnr(esc_fish3_m(i),osc_fish3_m(i),wt_fish3_size*double(nsamples_fish3_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_fish3_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) report << yrs_fish3_size(i)<<" "<<esc_fish3_m(i) <<endl; report<<endl;

  report << "Obs_full_P_fish3_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++)  {
     sdnr_fish3_size+= sdnr(esc_fish3_f(i),osc_fish3_f(i),wt_fish3_size*double(nsamples_fish3_size(i)))/nyrs_fish3_size/2;
  report << yrs_fish3_size(i)<<" "<<osc_fish3_f(i) 
      <<" eff_N "<<(1-esc_fish3_f(i))*esc_fish3_f(i)/norm2(osc_fish3_f(i)-esc_fish3_f(i))  <<" N "<<nsamples_fish3_size(i)
      <<" SDNR "<< sdnr(esc_fish3_f(i),osc_fish3_f(i),wt_fish3_size*double(nsamples_fish3_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_fish3_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) report << yrs_fish3_size(i)<<" "<<esc_fish3_f(i) <<endl; report<<endl;

  report << "Obs_full_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) {
   sdnr_srv1_age+=sdnr(eac_srv1(i),oac_srv1(i),wt_srv1_age*double(nsamples_srv1_age(i)))/nyrs_srv1_age;
   report << yrs_srv1_age(i)<<" "<<oac_srv1(i) 
      <<" eff_N "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i)) <<" N "<<nsamples_srv1_age(i)
      <<" SDNR "<< sdnr(eac_srv1(i),oac_srv1(i),wt_srv1_age*double(nsamples_srv1_age(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; report<<endl;

     report << "Obs_full_P_srv2_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv2_age;i++) {
    sdnr_srv2_age+=sdnr(eac_srv2(i),oac_srv2(i),wt_srv1_age*double(nsamples_srv2_age(i)))/nyrs_srv2_age;
    report << yrs_srv2_age(i)<<" "<<oac_srv2(i) 
      <<" eff_N "<<(1-eac_srv2(i))*eac_srv2(i)/norm2(oac_srv2(i)-eac_srv2(i)) <<" N "<<nsamples_srv2_age(i)
      <<" SDNR "<< sdnr(eac_srv2(i),oac_srv2(i),wt_srv1_age*double(nsamples_srv2_age(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv2_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv2_age;i++) report << yrs_srv2_age(i)<<" "<<eac_srv2(i) <<endl; report<<endl;

  report << "Obs_full_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) {
     sdnr_srv1_size+=sdnr(esc_srv1_f(i),osc_srv1_f(i),wt_srv1_size*double(nsamples_srv1_size(i)))/nyrs_srv1_size/2;
    report << yrs_srv1_size(i)<<" "<<osc_srv1_f(i) 
      <<" eff_N "<<(1-esc_srv1_f(i))*esc_srv1_f(i)/norm2(osc_srv1_f(i)-esc_srv1_f(i)) <<" N "<<nsamples_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1_f(i),osc_srv1_f(i),wt_srv1_size*double(nsamples_srv1_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<esc_srv1_f(i) <<endl; report<<endl;

    report << "Obs_full_P_srv1_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) {
    sdnr_srv1_size+=sdnr(esc_srv1_m(i),osc_srv1_m(i),wt_srv1_size*double(nsamples_srv1_size(i)))/nyrs_srv1_size/2;
     report << yrs_srv1_size(i)<<" "<<osc_srv1_m(i) 
      <<" eff_N "<<(1-esc_srv1_m(i))*esc_srv1_m(i)/norm2(osc_srv1_m(i)-esc_srv1_m(i)) <<" N "<<nsamples_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1_m(i),osc_srv1_m(i),wt_srv1_size*double(nsamples_srv1_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<esc_srv1_m(i) <<endl; report<<endl;

  report << "Obs_full_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) {
    sdnr_srv2_size+=sdnr(esc_srv2_f(i),osc_srv2_f(i),wt_srv2_size*double(nsamples_srv2_size(i)))/nyrs_srv2_size/2;
     report << yrs_srv2_size(i)<<" "<<osc_srv2_f(i) 
      <<" eff_N "<<(1-esc_srv2_f(i))*esc_srv2_f(i)/norm2(osc_srv2_f(i)-esc_srv2_f(i)) <<" N "<<nsamples_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2_f(i),osc_srv2_f(i),wt_srv2_size*double(nsamples_srv2_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<esc_srv2_f(i) <<endl; report<<endl;

    report << "Obs_full_P_srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) {
        sdnr_srv2_size+=sdnr(esc_srv2_m(i),osc_srv2_m(i),wt_srv2_size*double(nsamples_srv2_size(i)))/nyrs_srv2_size/2;
 report << yrs_srv2_size(i)<<" "<<osc_srv2_m(i) 
      <<" eff_N "<<(1-esc_srv2_m(i))*esc_srv2_m(i)/norm2(osc_srv2_m(i)-esc_srv2_m(i)) <<" N "<<nsamples_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2_m(i),osc_srv2_m(i),wt_srv2_size*double(nsamples_srv2_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<esc_srv2_m(i) <<endl; report<<endl;
  
    report << "Obs_full_P_srv7_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) {
    sdnr_srv7_size+= sdnr(esc_srv7_f(i),osc_srv7_f(i),wt_srv7_size*double(nsamples_srv7_size(i)))/nyrs_srv7_size/2;
     report << yrs_srv7_size(i)<<" "<<osc_srv7_f(i) 
      <<" eff_N "<<(1-esc_srv7_f(i))*esc_srv7_f(i)/norm2(osc_srv7_f(i)-esc_srv7_f(i)) <<" N "<<nsamples_srv7_size(i)
      <<" SDNR "<< sdnr(esc_srv7_f(i),osc_srv7_f(i),wt_srv7_size*double(nsamples_srv7_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv7_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) report << yrs_srv7_size(i)<<" "<<esc_srv7_f(i) <<endl; report<<endl;

    report << "Obs_full_P_srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) {
        sdnr_srv7_size+= sdnr(esc_srv7_m(i),osc_srv7_m(i),wt_srv7_size*double(nsamples_srv7_size(i)))/nyrs_srv7_size/2;
  report << yrs_srv7_size(i)<<" "<<osc_srv7_m(i) 
      <<" eff_N "<<(1-esc_srv7_m(i))*esc_srv7_m(i)/norm2(osc_srv7_m(i)-esc_srv7_m(i)) <<" N "<<nsamples_srv7_size(i)
      <<" SDNR "<< sdnr(esc_srv7_m(i),osc_srv7_m(i),wt_srv7_size*double(nsamples_srv7_size(i)))<<endl; report<<endl; }
  report << "Pred_full_P_srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) report << yrs_srv7_size(i)<<" "<<esc_srv7_m(i) <<endl; report<<endl;


  report << "Survey Biomass " <<endl;
  report << "Year:     " << yrs_srv1  <<endl;
  report << "Predicted:   " << pred_srv1  <<endl;
  report << "Observed:   " << obs_srv1_biom  <<endl<< endl;

  report << "Weight "<<aa<< endl;
  report << " "<< weight_f(endyr) << endl;

  report << "Weight Females "<<aa<< endl;
  report << " "<< weight_f(endyr) << endl;

  report << "Weight Males "<<aa<< endl;
  report << " "<< weight_m(endyr) << endl;

  report << "Fully_selected_F "<<aa <<endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     report << " "<< Fmort_fish1*max(fish4_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     report << " "<< Fmort_fish1*max(fish5_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    } 
   if(ph_ifq<1) // no post-IFQ sel est
    {
     report << " "<< Fmort_fish1*max(fish1_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }

  report << "Year " << yy<< endl;
  report << "SpBiom "<< spawn_biom <<endl;
  report << "Tot_biom "<< tot_biom   <<endl;
    // report << tot_biom.sd <<endl;

  report << "Fishery Selectivity " <<endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     report << fish4_sel_f / max(fish4_sel_f) <<endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     report << fish5_sel_f / max(fish5_sel_f) <<endl;
    }
   if(ph_ifq<1) // no post-IFQ sel est
    {
     report << fish1_sel_f / max(fish1_sel_f) <<endl;
    }
    
  report << "F35 F40 F50 "<<endl;
  report <<  F35 << " "<< mF40 <<" "<<  F50 <<endl;
  report << "SSB projection: "<< spawn_biom_proj<<endl<<"Catch projection: "<<pred_catch_proj<<endl;
  report <<  "B40: "<< B40<< endl;

  report << "Wts_n_Likelihoods  " << endl;
  report << wt_ssqcatch_fish1 <<" "<<ssqcatch <<" " ; report << "SSQ_Catch_Likelihood" << endl;
  report << wt_srv3     <<" "<<surv_like(3)  <<" " ; report << "Domestic_Survey_Abundance_Index_Likelihood" << endl;
  report << wt_srv4     <<" "<<surv_like(4)  <<" " ; report << "Cooperative_Survey_Abundance_Index_Likelihood" << endl;
  report << wt_fish1_age <<" "<<age_like(1)  <<" " ; report << "Fishery_Age_Composition_Likelihood" << endl;
  report << wt_srv1_age <<" "<<age_like(2)  <<" " ; report << "Survey_Age_Composition_Likelihood_DomesticLL" << endl;
  report << wt_srv1_age <<" "<<age_like(16)  <<" " ; report << "Survey_Age_Composition_Likelihood_CooperativeLL" << endl;
  report << wt_fish1_size<<" "<<age_like(3)+age_like(4)  <<" " ; report << "Fishery_Size_Composition_Likelihood_Fixed" << endl;
  report << wt_fish3_size<<" "<<age_like(6)+age_like(7)  <<" " ; report << "Fishery_Size_Composition_Likelihood_Trawl" << endl;
  report << wt_srv1_size<<" "<<age_like(9)+age_like(10)  <<" " ; report << "Survey_Size_Composition_Likelihood_Domestic LL" << endl;
  report << wt_srv2_size<<" "<<age_like(11)+age_like(12)  <<" " ; report << "Survey_Size_Composition_Likelihood_Cooperative_LL" << endl;
  report << wt_srv7_size<<" "<<age_like(13)+age_like(14)  <<" " ; report << "Survey_Size_Composition_Likelihood_GOATrawl" << endl;
  report << wt_rec_var <<" "<<rec_like     <<" " ; report << "Recruitment_Deviations_Likelihood" << endl;

  report << wt_sel_reg_fish1 <<" "<<sel_like(1)      <<" " ; report << "Fish_sel_Regularity_Penalty  "<<endl  ;
  report << wt_sel_reg_srv1 <<" "<<sel_like(2)      <<" " ; report << "Surv_sel_Regularity_Penalty  "<<endl  ;
  report << wt_sel_dome_fish1<<" "<<sel_like(3)      <<" " ; report << "Fish_Sel_Domeshapedness_Penalty "<<endl  ;
  report << wt_sel_dome_srv1<<" "<<sel_like(4)      <<" " ; report << "Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  report << "0"   <<" "<<avg_sel_penalty  <<" " ; report << "Average_Selectivity_Condition_Penalty "<<endl  ;
  report << wt_fmort_reg     <<" "<<F_mort_regularity<<" " ; report << "Fishing_Mortality_Regularity_Penalty" << endl;
  report << wt_M_reg     <<" "<<M_mort_regularity<<" " ; report << "Natural_Mortality_Regularity_Penalty" << endl;
 // report << wt_ssqcatch     <<" "<<F_dev_penalty<<" " ; report << "Fishing_Mortality_Deviations_Penalty" << endl;
  report << " "<<priors(1)  <<" " ; report << "priors sigr"     <<endl;
  report << " "<<priors(2)  <<" " ; report << "priors q_1" <<endl;
   report << " "<<priors(5)  <<" " ; report << "priors q_2" <<endl;
  report << " "<<priors(6)  <<" " ; report << "priors q_3" <<endl;
  report << " "<<priors(7)  <<" " ; report << "priors q_4" <<endl;
  report << " "<<priors(8)  <<" " ; report << "priors q_5" <<endl;
  report << " "<<priors(9)  <<" " ; report << "priors q_6" <<endl;
  report << " "<<priors(10)  <<" " ; report << "priors q_7" <<endl;
  report << " "<<priors(10)  <<" " ; report << "priors q_8" <<endl;
  report << " "<<priors(4)  <<" " ; report << "priors M"<<endl;
  report << " "<<obj_fun    <<" " ; report << "obj_fun"         <<endl;
  report << " "<<Like       <<" " ; report << "data likelihood" <<endl;//(2*square(sigr))+ size_count(log_rec_dev)*log(sigr)<<endl;
 if(last_phase()) { 
  report <<" SDNR1 "<< wt_srv1*std_dev(elem_div((pred_srv1(yrs_srv1)-obs_srv1_biom),obs_srv1_se))<<endl;
  report <<" SDNR2 "<< wt_srv2*std_dev(elem_div((pred_srv2(yrs_srv2)-obs_srv2_biom),obs_srv2_se))<<endl;
  report <<" SDNR3 "<< wt_srv3*std_dev(elem_div((pred_srv3(yrs_srv3)-obs_srv3_biom),obs_srv3_se))<<endl;
  report <<" SDNR4 "<< wt_srv4*std_dev(elem_div((pred_srv4(yrs_srv4)-obs_srv4_biom),obs_srv4_se))<<endl;
  report <<" SDNR5 "<< wt_srv5*std_dev(elem_div((pred_srv5(yrs_srv5)-obs_srv5_biom),obs_srv5_se))<<endl;
  report <<" SDNR6 "<< wt_srv6*std_dev(elem_div((pred_srv6(yrs_srv6)-obs_srv6_biom),obs_srv6_se))<<endl;
  report <<" SDNR7 "<< wt_srv7*std_dev(elem_div((pred_srv7(yrs_srv7)-obs_srv7_biom),obs_srv7_se))<<endl;
//  report <<" SDNR8 "<< wt_srv8*std_dev(elem_div((pred_srv8(yrs_srv8)-obs_srv8_biom),obs_srv8_se))<<endl;
 

    }
  report << "SigmaR: "<<sigr<< " Nat_Mort_Base: "<<natmort<<" Male delta M "<<mdelta<<" Spawning Per Recruit "<< " "<<SBF40<<" "<<SB0<<" Virgin SPR "<<endl;
  report << "Stock-recruitment, type: "<<SrType<<" 1=Ricker, 2=B-Holt, 3=Mean"<<endl;
  report << "Year SSB SR_Pred R_Est "<<endl;
  report << "Weighted likelihods for comps broken down"<<endl;
  report << "Fish1 Age, Survey 1 Age, Fish 1 Size Female, Fish 1 Size Male, Fish 2 Size, Fish3 Size Female, F3 Size Male, Srv1 Size F, Srv1 Size M, Srv2 Size F, Srv2 Size M, Srv7 Size F, Srv7 Size M"<<endl;
  report <<age_like<<" age _like"<<endl<<surv_like<<" surv like "<<endl<<sel_like<<" sel_like"<<endl;
  report << "Unweighted likelihods for comps broken down"<<endl;
  report << "Fish1 Age, Survey 1 Age, Fish 1 Size Female, Fish 1 Size Male, Fish 2 Size, Fish3 Size Female, F3 Size Male, Srv1 Size F, Srv1 Size M, Srv2 Size F, Srv2 Size M, Srv7 Size F, Srv7 Size M"<<endl;
  report <<age_like(1) / (0.00001+wt_fish1_age)<< " "<<(age_like(2)+0.0001) / (0.00001+wt_srv1_age)<< " "<<(age_like(3)+0.0001) / (0.00001+wt_fish1_size)<< " "<<(age_like(4)+0.0001) / (0.00001+wt_fish1_size)<< " "<<(age_like(5)+0.0001) / (0.00001+wt_fish2_size)<< " "<<(age_like(6)+0.0001) / (0.00001+wt_fish3_size)<< " "<<(age_like(7)+0.0001) / (0.00001+wt_fish3_size)<< " "<<(age_like(8)+0.0001) / (0.00001+wt_fish4_size)<< " "<<(age_like(9)+0.0001) / (0.00001+wt_srv1_size)<< " "<<(age_like(10)+0.0001) / (0.00001+wt_srv1_size)<< " "<<(age_like(11)+0.0001) / (0.00001+wt_srv2_size)<< " "<<(age_like(12)+0.0001) / (0.00001+wt_srv2_size)<< " "<<(age_like(13)+0.0001) / (0.00001+wt_srv7_size)<< " "<<(age_like(14)+0.0001) / (0.00001+wt_srv7_size)<<" "<<(age_like(15)+0.0001) / (0.00001+wt_srv7_age)<< " "<<(age_like(16)+0.0001) / (0.00001+wt_srv1_age) <<" unw_age)_like"<<endl;
  report <<surv_like(1) / (0.00001+wt_srv1)<< " "<<surv_like(2) / (0.00001+wt_srv2)<< " "<<surv_like(3) / (0.00001+wt_srv3)<< " "<<surv_like(4) / (0.00001+wt_srv4)<< " "<<surv_like(5) / (0.00001+wt_srv5)<< " "<<surv_like(6) / (0.00001+wt_srv6)<< " "<<(surv_like(7)+0.0001) / (0.00001+wt_srv7)<< " "<<(surv_like(8)+0.0001) / (0.00001+wt_srv8)<< " "<<"unw_surv_like "<<endl<<sel_like<<" sel_like"<<endl;
  report << "Survey Biomass 2" <<endl;
  report << "Year:     " << yrs_srv2  <<endl;
  report << "Predicted:   " << pred_srv2  <<endl;
  report << "Observed:   " << obs_srv2_biom  <<endl<< endl;
  report << "Survey Biomass 3" <<endl;
  report << "Year:     " << yrs_srv3  <<endl;
  report << "Predicted:   " << pred_srv3  <<endl;
  report << "Observed:   " << obs_srv3_biom  <<endl<< endl;

  report << "Survey Biomass 4" <<endl;
  report << "Year:     " << yrs_srv4  <<endl;
  report << "Predicted:   " << pred_srv4  <<endl;
  report << "Observed:   " << obs_srv4_biom  <<endl<< endl;

  report << "Survey Biomass 5" <<endl;
  report << "Year:     " << yrs_srv5  <<endl;
  report << "Predicted:   " << pred_srv5  <<endl;
  report << "Observed:   " << obs_srv5_biom  <<endl<< endl;

  report << "Survey Biomass 6" <<endl;
  report << "Year:     " << yrs_srv6  <<endl;
  report << "Predicted:   " << pred_srv6  <<endl;
  report << "Observed:   " << obs_srv6_biom  <<endl<< endl;
  report << "Survey Biomass 7" <<endl;
  report << "Year:     " << yrs_srv7  <<endl;
  report << "Predicted:   " << pred_srv7  <<endl;
  report << "Observed:   " << obs_srv7_biom  <<endl<< endl;
  report << "qs"<<endl<<q_srv1<<endl<<q_srv2<<endl<<q_srv3<<endl<<q_srv4<<endl<<q_srv5<<endl<<q_srv6<<endl<<q_srv7<<endl<<q_srv8<<endl;
  report << "Year SSB SRR Recr "<<endl;
  for (i=styr;i<=endyr;i++)
  report<< i <<" "<<Sp_Biom(i-recage)<<" "<<srm_rec(i)<<" "<<sam_rec(i)<<endl;
  report<< "Age/Length residuals"<<endl;
  report << "fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) report << yrs_fish1_age(i)<<" "<<oac_fish1(i)-eac_fish1(i)<<endl; 

  report << "fish1_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) report << yrs_fish1_size(i)<<" "<<osc_fish1_f(i)-esc_fish1_f(i)<<endl; 

  report << "fish1_size Male"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) report << yrs_fish1_size(i)<<" "<<osc_fish1_m(i)-esc_fish1_m(i)<<endl; 

  report << "fish3_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) report << yrs_fish3_size(i)<<" "<<osc_fish3_f(i)-esc_fish3_f(i)<<endl;
   
  report << "fish3_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) report << yrs_fish3_size(i)<<" "<<osc_fish3_m(i)-esc_fish3_m(i)<<endl;

  report << "srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<oac_srv1(i)-eac_srv1(i)<<endl; 

  report << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<osc_srv1_f(i)-esc_srv1_f(i)<<endl;

  report << "srv1_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<osc_srv1_m(i)-esc_srv1_m(i)<<endl;
  
  report << "srv2_size Females"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<osc_srv2_f(i)-esc_srv2_f(i)<<endl;

  report << "srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<osc_srv2_m(i)-esc_srv2_m(i)<<endl;

    report << "srv7_size Females"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) report << yrs_srv7_size(i)<<" "<<osc_srv7_f(i)-esc_srv7_f(i)<<endl;

  report << "srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) report << yrs_srv7_size(i)<<" "<<osc_srv7_m(i)-esc_srv7_m(i)<<endl;

  report << "N_proj_f "<<endl<<N_proj_f<<endl<<"N_proj_m "<<N_proj_m<<endl<<" spawn_bio next year"<<endl<<spawn_biom_proj(endyr+1)<<endl
  <<"ABC full"<<ABC3<<endl;
  
  report <<" spawn_bio projected"<<endl<<spawn_biom_proj<<endl;
  report << "Specified catch projection"<<(catage_proj_f*weight_f(endyr)+catage_proj_m*weight_m(endyr))<<endl<<"ABC projection: "<<pred_catch_proj<<endl;
  // New SDNR stuff
  report<<" average SDNRs for compositional data"<<endl;
  report<<"Fishery Ages "<<sdnr_fish1_age<<endl;
  report<<"Fishery 1 sizes "<<sdnr_fish1_size<<endl;
  report<<"Fishery 3 sizes "<<sdnr_fish3_size<<endl;
  report<<"Survey 1 ages "<<sdnr_srv1_age<<endl;
  report<<"Survey 2 ages "<<sdnr_srv2_age<<endl;
  report<<"Survey 1 sizes "<<sdnr_srv1_size<<endl;
  report<<"Survey 2 sizes "<<sdnr_srv2_size<<endl;
  report<<"Survey 7 sizes "<<sdnr_srv7_size<<endl;
  write_sarareport();
    save_gradients(gradients);
  # include "sable-r-report.cxx"

FUNCTION write_sarareport
   ofstream sarareport("SABLE_SARA.dat");

 sarareport << "SABLE        # stock  " << endl;
 sarareport << "AK       # region     (AI AK BOG BSAI EBS GOA SEO WCWYK)" << endl;
 sarareport << endyr << "       # ASSESS_YEAR - year assessment is presented to the SSC" << endl;
 sarareport << "3a         # TIER  (1a 1b 2a 2b 3a 3b 4 5 6) " << endl;
 sarareport << "none       # TIER2  if mixed (none 1a 1b 2a 2b 3a 3b 4 5 6)" << endl;
 sarareport << "Update       # UPDATE (new benchmark full partial)" << endl;
  sarareport << "!!!!REPLACE WITH LOW 2.5% SSB FROM MCMC" << "     # Minimum B  Lower 95% confidence interval for spawning biomass (kt) in assessment year" << endl;
 sarareport << "!!!!REPLACE WITH high 97.5% SSB FROM MCMC" << "     # Maximum B  Upper 95% confidence interval for spawning biomass (kt) in assessment year" << endl;
 sarareport << "!!!!Replace with B35 from SUmmary Table" << "     # BMSY (kt)  is equilibrium spawning biomass at MSY (Tiers 1-2) or 7/8 x B40% (Tier 3)" << endl;
 sarareport << "Custom AMAK       # MODEL - Required only if NMFS toolbox software used; optional otherwise " << endl;
 sarareport << "NA         # VERSION - Required only if NMFS toolbox software used; optional otherwise" << endl;
 sarareport << "2          # number of sexes  if 1 sex=ALL elseif 2 sex=(FEMALE, MALE) " << endl;
 sarareport << "2          # number of fisheries" << endl;
 sarareport << "1000000          # multiplier for recruitment, N at age, and survey number (1,1000,1000000)" << endl;
 sarareport << "2          # recruitment age used by model or size" << endl;
 sarareport << "2          # age+ or mmCW+ used for biomass estimate" << endl;
 sarareport << "\"Single age\"        # Fishing mortality type such as \"Single age\" or \"exploitation rate\"" << endl;
 sarareport << "\"Age model\"         # Fishing mortality source such as \"Model\" or \"(total catch (t))/(survey biomass (t))\"" << endl;
 sarareport << "\"Age of maximum F\"  # Fishing mortality range such as \"Age of maximum F\"" << endl; 
 sarareport << "#FISHERYDESC -list of fisheries (ALL TWL LGL POT FIX FOR DOM TWLJAN LGLMAY POTAUG ...)" << endl; 
 sarareport << "FIX TWL" << endl; 
 sarareport <<"#FISHERYYEAR - list years used in the model " << endl;
 sarareport << yy << endl; 
 sarareport<<"#AGE - list of ages used in the model"<<endl;
 sarareport << aa << endl; 
 sarareport <<"#RECRUITMENT - Number of recruits by year " << endl;
 sarareport  << pred_rec << endl;     
 sarareport <<"#SPAWNBIOMASS - Spawning biomass by year in kilo tons " << endl;
 sarareport  << spawn_biom << endl;  
 sarareport <<"#TOTALBIOMASS - Total biomass by year in kilo tons " << endl;
 sarareport  << tot_biom << endl;
 sarareport <<"#TOTFSHRYMORT - Fishing mortality rate by year " << endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     sarareport  << Fmort_fish1*max(fish4_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     sarareport  << Fmort_fish1*max(fish5_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    } 
   if(ph_ifq<1) // no post-IFQ sel est
    {
     sarareport  << Fmort_fish1*max(fish1_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }
 sarareport <<"#TOTALCATCH - Total catch by year in kilo tons " << endl;
 sarareport  << pred_catch_fish1 + pred_catch_fish3 << endl;
 sarareport <<"#MATURITY - Maturity ratio by age (females only)" << endl;  
 sarareport  << maturity(endyr) << endl; 
 sarareport <<"#SPAWNWT - Average spawning weight (in kg) by age"<< endl; 
 sarareport << weight_maturity_prod_f(endyr) << endl;
 sarareport <<"#NATMORT - Natural mortality rate for females then males"<< endl; 
  for (i=1;  i<=nages;  i++) 
   sarareport  << M(endyr,i) <<"  ";
   sarareport<< endl;   
 for (i=1;  i<=nages;  i++) 
   sarareport  << M(endyr,i) <<"  ";
   sarareport<< endl;   
    sarareport << "#N_AT_AGE - Estimated numbers of female (first) then male (second) fish at age " << endl;
  for (i=styr; i<=endyr;i++)
   sarareport <<natage_f(i)<< "  ";
   sarareport<<endl;
  for (i=styr; i<=endyr;i++)
   sarareport <<natage_m(i)<< "  ";
   sarareport<<endl;
 sarareport <<"#FSHRY_WT_KG - Fishery weight at age (in kg) females (first) males (second), only one fishery"<< endl;   
   sarareport << weight_f(endyr)  << endl;
   sarareport << weight_m(endyr) <<endl;
     sarareport << weight_f(endyr)  << endl;
   sarareport << weight_m(endyr) <<endl;
   sarareport<<endl;       
 sarareport << "#SELECTIVITY - Estimated fixed gear fishery selectivity, females first, males next" << endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     sarareport << fish4_sel_f << endl;
     sarareport << fish4_sel_m << endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     sarareport << fish5_sel_f << endl;
     sarareport << fish5_sel_m << endl;
    } 
   if(ph_ifq<1) // no post-IFQ sel est
    {
     sarareport << fish1_sel_f << endl;
     sarareport << fish1_sel_m << endl;
    }
 sarareport << "#SELECTIVITY - Estimated trawl gear fishery selectivity, females first, males next" << endl;
 sarareport << fish3_sel_f << endl;
 sarareport << fish3_sel_m << endl;
  sarareport << "#SURVEYDESC"<<endl;
 sarareport<<"GOA_trawl_survey"<<endl;
 sarareport<<"#SURVEYMULT"<<endl;
 sarareport<<"1000"<<endl;
 sarareport << "#GOA_trawl_survey - Gulf of Alaska survey biomass (Year, Obs_biomass) " << endl;
 sarareport << yrs_srv7 << endl;
 sarareport<< obs_srv7_biom << endl;
   sarareport << "#SURVEYDESC"<<endl;
   sarareport<<"AK_LL_survey_cooperative"<<endl;
 sarareport<<"#SURVEYMULT"<<endl;
 sarareport<<"1000"<<endl;
 sarareport << "#AK Longline survey - cooperative with Japanese relative population numbers " << endl;
 sarareport << yrs_srv4 << endl;
 sarareport<< obs_srv4_biom << endl;
   sarareport << "#SURVEYDESC"<<endl;
 sarareport<<"AK_LL_survey"<<endl;
 sarareport<<"#SURVEYMULT"<<endl;
 sarareport<<"1000"<<endl;
 sarareport << "#AK Longline survey - relative population numbers " << endl;
 sarareport << yrs_srv3 << endl;
 sarareport<< obs_srv3_biom << endl;
   sarareport << "#SURVEYDESC"<<endl;
  sarareport<<"AK_LL_fishery CPUE"<<endl;
 sarareport<<"#SURVEYMULT"<<endl;
 sarareport<<"1000"<<endl;
 sarareport << "#AK Longline fishery CPUE - relative population weight " << endl;
 sarareport << yrs_srv5 << endl;
 sarareport<< obs_srv5_biom << endl;
 sarareport<<"#STOCKNOTES"<<endl;
 sarareport<<"\"SAFE report indicates that this stock was not subjected to overfishing in "<<endyr-1<< " and is neither overfished nor approaching a condition of being overfished in "<<endyr<<".\""<<endl;
 sarareport<<"\"A full assessment was undertaken in "<<endyr<<" using  model 23.5.\""<<endl;


FUNCTION write_HQreport
   ofstream HQreport("SABLEAK_HQ.dat");
 HQreport<<"#ASSESSMENT_SUMMARY -------------------------------------------------------------------------------------------------"<<endl;
 HQreport << "#STOCK  " << endl;
 HQreport << "\"SABLE\"" << endl;
 HQreport << "#STOCK_NAME"<<endl;
 HQreport << "\"Sablefish - Eastern Bering Sea / Aleutian Islands / Gulf of Alaska\"" << endl;
  HQreport << "#REGION" <<endl;
 HQreport << "\"AK\"" << endl;
  HQreport << "#ASMT_TYPE"<<endl;
 HQreport << "\"Operational\"" << endl;
  HQreport << "#ASMT_YEAR"<<endl;
 HQreport << endyr << endl;
  HQreport <<"#ASMT_MONTH"<<endl;
  HQreport << "\"Dec\"" << endl;
   HQreport <<"#TIER"<<endl;
 HQreport << "\"3a\"" << endl;
 HQreport <<"#NUM_SEXES"<<endl;
 HQreport << 2  << endl;
  HQreport <<"#NUM_FISHERIES"<<endl;
 HQreport <<2 << endl;
  HQreport <<"#REC_MULT"<<endl;
 HQreport << 1000000 << endl;
  HQreport <<"#RECAGE"<<endl;
 HQreport << 2  << endl;
  HQreport <<"#COMPLEX"<<endl;
  HQreport << "\"NA\"" << endl;
   HQreport <<"#LAST_DATA_YEAR"<<endl;
   HQreport <<  endyr  <<endl;
  HQreport <<"#ASMT_MODEL_CATEGORY"<<endl;
   HQreport <<"\"6 - Statistical Catch-at-Age\""<<endl;
   HQreport <<"#ASMT_MODEL"<<endl;
   HQreport <<"\"Custom SCAA: Custom Statistical Catch-at-Age\""<<endl;
   HQreport <<"#MODEL_VERSION"<<endl;
   HQreport <<"\"23.5\""<<endl;
   HQreport <<"#ENSEMBLE"<<endl;
   HQreport <<"\"NA\""<<endl;
   HQreport <<"#LEAD_LAB"<<endl;
   HQreport <<"\"AFSC\""<<endl;
   HQreport <<"#POC_EMAIL"<<endl;
   HQreport <<"\"daniel.goethel@noaa.gov\""<<endl;
   HQreport <<"#REVIEW_RESULT"<<endl;
   HQreport <<"\"Full acceptance\""<<endl;
   HQreport <<"#CATCH_INPUT_DATA"<<endl;
   HQreport <<5<<endl;
   HQreport <<"#ABUNDANCE_INPUT_DATA"<<endl;
   HQreport <<4<<endl;
   HQreport <<"#BIOLOGICAL_INPUT_DATA"<<endl;
   HQreport <<4<<endl;
   HQreport <<"#SIZEAGE_COMP_INPUT_DATA"<<endl;
   HQreport <<4<<endl;
   HQreport <<"#ECOSYSTEM_LINKAGE"<<endl;
   HQreport <<2<<endl;

 HQreport<<"#FISHING_MORTALITY_ESTIMATES ----------------------------------------------------------------------------------------"<<endl;
   HQreport <<"#F_YEAR"<<endl;
   HQreport <<endyr-1<<endl;
   HQreport <<"#F_BASIS"<<endl;
   HQreport<<"\"F for Fully-Selected Fish\""<<endl;
   HQreport <<"#F_UNIT"<<endl;
   HQreport<<"\"Fully-Selected F\""<<endl;
   HQreport <<"#BEST_F_ESTIMATE"<<endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     HQreport  << Fmort_fish1(endyr-1)*max(fish4_sel_f)+Fmort_fish3(endyr-1)*max(fish3_sel_f) <<endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     HQreport  << Fmort_fish1(endyr-1)*max(fish5_sel_f)+Fmort_fish3(endyr-1)*max(fish3_sel_f) <<endl;
    } 
   if(ph_ifq<1) // no post-IFQ sel est
    {
     HQreport  << Fmort_fish1(endyr-1)*max(fish1_sel_f)+Fmort_fish3(endyr-1)*max(fish3_sel_f) <<endl;
    }
   HQreport <<"#F_LIMIT"<<endl;
   HQreport <<"!!!!!REPLACE WITH F SIS VALUE FROM SPREADSHEET"<<endl;
   HQreport <<"#F_LIMIT_BASIS"<<endl;
   HQreport<<"\"F from 2023 asmt corresponding to specified 2022 OFL\""<<endl;
   HQreport <<"#F_MSY"<<endl;
   HQreport<<FOFL<<endl;
   HQreport <<"#F_MSY_BASIS"<<endl;
   HQreport<<"\"F35% as proxy\""<<endl;

 HQreport<<"#BIOMASS_ESTIMATES --------------------------------------------------------------------------------------------------"<<endl;

   HQreport <<"#B_YEAR"<<endl;
   HQreport <<endyr<<endl;
   HQreport <<"#B_BASIS"<<endl;
   HQreport <<"\"Mature Female Biomass\""<<endl;
   HQreport <<"#B_UNIT"<<endl;
   HQreport <<"\"Metric tons\""<<endl;
   HQreport <<"#BEST_B_ESTIMATE"<<endl;
   HQreport <<spawn_biom(endyr)*1000<<endl;
   HQreport <<"#LOWER_B_ESTIMATE"<<endl;
  HQreport << "!!!!REPLACE WITH LOW 2.5% SSB FROM MCMC" <<endl;
   HQreport <<"#UPPER_B_ESTIMATE"<<endl;
 HQreport << "!!!!!REPLACE WITH high 97.5% SSB FROM MCMC" <<endl;
    HQreport <<"#ESTIMATE_METHOD"<<endl;
   HQreport <<"\"Credible\""<<endl;
    HQreport <<"#INTERVAL_SIZE"<<endl;
   HQreport <<95<<endl;
    HQreport <<"#B_MSY"<<endl;
   HQreport << "!!!!Replace with B_35 from Summary Table"<<endl;
    HQreport <<"#B_MSY_BASIS"<<endl;
   HQreport <<"\"B35%\""<<endl;

 HQreport<<"#TIME_SERIES_ESTIMATES ----------------------------------------------------------------------------------------------"<<endl;

 HQreport <<"#FISHERYYEAR" << endl;
 HQreport << yy << endl;
  HQreport<<"#AGE"<<endl;
 HQreport << aa << endl;
 HQreport <<"#RECRUITMENT" << endl;
 HQreport  << pred_rec << endl;     
 HQreport <<"#SPAWNBIOMASS" << endl;
 HQreport  << spawn_biom*1000 << endl;  
 HQreport <<"#TOTALBIOMASS" << endl;
 HQreport  << tot_biom*1000 << endl;
 HQreport <<"#TOTFSHRYMORT" << endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     HQreport  << Fmort_fish1*max(fish4_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     HQreport  << Fmort_fish1*max(fish5_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    } 
   if(ph_ifq<1) // no post-IFQ sel est
    {
     HQreport  << Fmort_fish1*max(fish1_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }
 HQreport <<"#TOTALCATCH" << endl;
 HQreport  << (pred_catch_fish1 + pred_catch_fish3)*1000 << endl;

  HQreport << "#SURVEYDESC"<<endl;
 HQreport<<"GOA_trawl_survey; AK Cooperative Longline survey; AK US Longline survey; AK Longline fishery CPUE"<<endl;

 HQreport<<"#STOCKNOTES"<<endl;
 HQreport<<"\"SAFE report indicates that this stock was not subjected to overfishing in "<<endyr-1<< " and is neither overfished nor approaching a condition of being overfished in "<<endyr<<".\""<<endl;
 HQreport<<"\"A full assessment was undertaken in "<<endyr<<" using  model 23.5.\""<<endl;

 HQreport<<"#SURVEY_ESTIMATES [OPTIONAL] ------------------------------------------------------------------------------------------"<<endl;

 HQreport << "#GOA_trawl_survey" << endl;
 HQreport << yrs_srv7 << endl;
 HQreport<< obs_srv7_biom << endl;
 HQreport<<"#AK Cooperative Longline survey"<<endl;
 HQreport << yrs_srv4 << endl;
 HQreport<< obs_srv4_biom << endl;
 HQreport<<"#AK US Longline survey"<<endl;
 HQreport << yrs_srv3 << endl;
 HQreport<< obs_srv3_biom << endl;
  HQreport<<"#AK Longline fishery CPUE"<<endl;
 HQreport << yrs_srv5 << endl;
 HQreport<< obs_srv5_biom << endl;

FUNCTION double sdnr(const dvar_vector& pred,const dvector& obs,double m)
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred)+0.000001;
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs+0.000001-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;

FUNCTION write_projout
 ofstream projout("projold.dat");
// Function to write out data file for projection model....
 projout <<"#_Random_number_seed"<<endl;
 projout <<"1234"<<endl;
 projout <<"#_Number_of_fisheries"<<endl;
 projout <<"1"<<endl;
 projout <<"#_Number_of_projection_years"<<endl;
 projout <<"14"<<endl;
 projout <<"#_Number_of_simulations"<<endl;
 projout <<"1000"<<endl;
 projout <<"#_Begin_year_of_projection" <<endl;
 projout <<endyr<<endl;
 projout <<"#_Number_of_ages"<<endl;
 projout <<nages<<endl;
 for (j=1;j<=nages;j++) natmortv(j) =M(endyr,j); 
 projout <<"#_Natural_Mortality" << aa << endl;
 projout <<natmortv<<endl;
 projout <<"#_2001_TAC_or_best_estimate_of_catch_2001"<<endl;
 projout <<obs_catch_fish1(endyr)<<endl;
 projout <<"#_F_ratio(must_sum_to_one_only_one_fishery)"<<endl;
 projout <<"1"<<endl;
 projout <<"#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)"<<endl;
 projout << mean(Fmort_fish1(endyr-4,endyr))<<endl;
 projout <<"#_Author_F_as_fraction_F_40%"<<endl;
 projout <<"1"<<endl;
 projout <<"#_Spawn_month"<<endl;
 projout << spawn_fract*12+1<<endl;
 projout <<"#_Wt_at_age_spawners"<<aa<<endl<<weight_f(endyr)<< endl;
 projout <<"#_Wt_at_age_fishery" <<aa<<endl<<weight_f(endyr)<< endl;
 projout <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2"<<aa<<endl<<elem_prod(elem_div(natage_f(endyr),(natage_f(endyr)+natage_m(endyr))),maturity(endyr))<< endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     sel_rep_proj_f=fish4_sel_f / max(fish4_sel_f);
     sel_rep_proj_m=fish4_sel_m / max(fish4_sel_m);
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     sel_rep_proj_f=fish5_sel_f / max(fish5_sel_f);
     sel_rep_proj_m=fish5_sel_m / max(fish5_sel_m);
    }
   if(ph_ifq<1) // no post-IFQ sel est
    {
     sel_rep_proj_f=fish1_sel_f / max(fish1_sel_f);
     sel_rep_proj_m=fish1_sel_m / max(fish1_sel_m);
    }
 projout <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<sel_rep_proj_f<< endl;
 projout <<"#_Numbers_at_age_end_year"<<aa<<endl<<natage_f(endyr)+natage_m(endyr)<< endl;
 projout <<"#_N_recruitment_years"<<endl<<endyr-1979-1<< endl;
 projout <<"#_Recruitment_start_at_1977_yearclass=1979_for_age_2_recruits"<<yy(1979,endyr-2)<<endl<<pred_rec(1979,endyr -2)<< endl;
 projout.close();

FUNCTION write_newproj
 ofstream newproj("proj.dat");
// Function to write out data file for new Ianelli 2005 projection model....
 newproj <<"#Species name here:"<<endl;
 newproj <<"GOA_SABLE"<<endl;
 newproj <<"#SSL Species?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Constant buffer of Dorn?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Number of fisheries?"<<endl;
 newproj <<"2"<<endl;
 newproj <<"#Number of sexes?"<<endl;
 newproj <<"2"<<endl;
 newproj <<"#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)"<<endl;
 newproj << mean(Fmort_fish1(endyr-4,endyr)+Fmort_fish3(endyr-4,endyr))<<endl;
 newproj <<"#_Author_F_as_fraction_F_40%"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#ABC SPR" <<endl;
 newproj <<"0.4"<<endl;
 newproj <<"#MSY SPR" <<endl;
 newproj <<"0.35"<<endl;
 newproj <<"#_Spawn_month"<<endl;
 newproj << spawn_fract*12+1<<endl;
 newproj <<"#_Number_of_ages"<<endl;
 newproj <<nages<<endl;
 newproj <<"#_F_ratio(must_sum_to_one_only_one_fishery)"<<endl;
 newproj <<fratio<<" "<<1-fratio<<endl;
 for (j=1;j<=nages;j++) natmortv(j) = M(endyr,j); 
 newproj <<"#_Natural_Mortality" << aa << endl;
 newproj <<natmortv<<endl;
 newproj <<"#_Natural_Mortality" << aa << endl;
 newproj <<natmortv<<endl;
 newproj <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2"<<aa<<endl<<maturity(endyr)<< endl;
 newproj <<"#_Wt_at_age_spawners"<<aa<<endl<<weight_f(endyr)<< endl;
 newproj <<"#_Wt_at_age_fishery1 female" <<aa<<endl<<weight_f(endyr)<< endl;
 newproj <<"#_Wt_at_age_fishery2 female" <<aa<<endl<<weight_f(endyr)<< endl;
 newproj <<"#_Wt_at_age_spawners1 male"<<aa<<endl<<weight_m(endyr)<< endl;
 newproj <<"#_Wt_at_age_fishery2 male" <<aa<<endl<<weight_m(endyr)<< endl;
 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<sel_rep_proj_f<< endl<<fish3_sel_f/max(fish3_sel_f)<<endl;
 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<sel_rep_proj_m<< endl<<fish3_sel_m/max(fish3_sel_m)<<endl;
 newproj <<"#_Numbers_at_age_end_year"<<aa<<endl<<natage_f(endyr)<<endl<<natage_m(endyr)<< endl;
 newproj <<"#_N_recruitment_years"<<endl<<endyr-1979-1<< endl;
 newproj <<"#_Recruitment_start_at_1977_yearclass=1979_for_age_2_recruits"<<yy(1979,endyr-2)<<endl<<pred_rec(1979,endyr-2)<< endl;
 newproj <<"#_Spawners per recruitment (starting at 1977)"<<endl<<spawn_biom(1977,endyr-5)<< endl;
 newproj.close();

FUNCTION write_fullrep
     ofstream fullrep("sable.rep");

  fullrep<<"****Executive mary Material*****"<<endl;
  fullrep<<"     Model name"     <<endl;
  fullrep<<model_name<<endl;
  fullrep<<"     .dat file"     <<endl;
  fullrep<<data_file<<endl;
  fullrep<<"     Number parameters estimated"     <<endl;
  fullrep<<initial_params::nvarcalc()<<endl;
  fullrep<<"     TotalBiomass for "<<endyr+1<<endl;
  fullrep<<tot_biom_proj(endyr+1)<<endl;
  fullrep<<"     TotalBiomass for "<<endyr+2     <<endl;
  fullrep<<tot_biom_proj(endyr+2)<<endl;
  fullrep<<"     Female_Spawning Biomass for "<<endyr+1     <<endl;
  fullrep<<spawn_biom_proj(endyr+1)<<endl;
  fullrep<<"     Female_Spawning_Biomass for "<<endyr+2     <<endl;
  fullrep<<spawn_biom_proj(endyr+2)<<endl;
  fullrep<<"     B_100"     <<endl;
  fullrep<<SB0*mean(pred_rec(1979,endyr-recage))/2<<endl;
  fullrep<<"     B_40"     <<endl;
  fullrep<<B40<<endl;
  fullrep<<"     B_35"     <<endl;
  fullrep<<SBF35*mean(pred_rec(1979,endyr-recage))/2<<endl;
  fullrep<<"     Mean_Recruitment"     <<endl;
  fullrep<<mean(pred_rec(1979,endyr-recage))<<endl;
  fullrep<<"     F_40"     <<endl;
  fullrep<<F40<<endl;
  fullrep<<"     F_35"     <<endl;
  fullrep<<F35<<endl;
  fullrep<<"     F_ABC for "<<endyr+1     <<endl;
  fullrep<<FABC<<endl;
  fullrep<<"     F_ABC for "<<endyr+2     <<endl;
  fullrep<<FABC2<<endl;
  fullrep<<"     ABC for "<<endyr+1     <<endl;
  fullrep<<pred_catch_proj(endyr+1)<<endl;
  fullrep<<"     ABC for "<<endyr+2     <<endl;
  fullrep<<pred_catch_proj(endyr+2)<<endl;
  fullrep<<"     F_OFL for "<<endyr+1     <<endl;
  fullrep<<FOFL<<endl;
  fullrep<<"     F_OFL for "<<endyr+2     <<endl;
  fullrep<<FOFL2<<endl;
  fullrep<<"     OFL for "<<endyr+1     <<endl;
  fullrep<<OFL<<endl; 
  fullrep<<"     OFL for "<<endyr+2     <<endl;
  fullrep<<pred_catch_proj_OFL(endyr+2)<<endl; 
  fullrep<<"     Total likelihood"     <<endl;
  fullrep<<obj_fun<<endl;
  fullrep<<"     Data likelihood"     <<endl;
  fullrep<<Like<<endl<<endl;
  
  fullrep<<" ************   Some more parameter estimates and their SDs ************"<<endl;
 
  if(last_phase()) {
    // add standard deviation data types    
  fullrep<<"   q_domestic   "<<endl;
  fullrep<<q_srv1<<" "<<q_srv1.sd<<endl;
  fullrep<<"   q_cooperative  "<<endl;
   fullrep<<q_srv2<<" "<<q_srv2.sd<<endl;
  fullrep<<"   q_trawl  "<<endl;
   fullrep<<q_srv7<<" "<<q_srv7.sd<<endl;
 /*
  fullrep<<natmort<<" "<<nattymort.sd<<endl;
  fullrep<<"  sigr   "<<endl;  
  fullrep<<sigr<<" "<<cigar.sd<<endl;  
  fullrep<<"   log_mean_rec"<<endl;
  fullrep<<log_mean_rec<<" "<<LMR.sd<<endl;
 */
  fullrep<<"   F_40"<<endl;
  fullrep<<F40<<" "<<F40.sd<<endl;
  fullrep<<"    tot_biom"<<endl;
  fullrep<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  fullrep<<"   spawn_biom"<<endl;
  fullrep<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  fullrep<<"    B40"<<endl;
  fullrep<<B40<<" "<<B40.sd<<endl;
  fullrep<<"   ABC"<<endl;
  fullrep<<pred_catch_proj(endyr+1)<<" "<<pred_catch_proj.sd(endyr+1)<<endl<<endl;
 
 }

  fullrep << "Parameter Phases"<<endl;
  fullrep << "ph_mean_rec ph_recdev ph_Rzero ph_steepness ph_sigr" <<endl;
  fullrep <<  ph_mean_rec << " "<< ph_recdev <<" "<<  ph_Rzero<<" "<<ph_steepness<<" "<<ph_sigr  <<endl;
  fullrep << "ph_avg_F ph_Fdev ph_F50 "<<endl;
  fullrep <<  ph_avg_F << " "<< ph_Fdev <<" "<<  ph_F50 <<endl;
  fullrep << "ph_m ph_Mdevs ph_Mdevs_age ph_mdelta"<<endl;
  fullrep <<  ph_m << " "<< ph_Mdevs <<" "<<  ph_Mdevs_age<<" "<<ph_mdelta  <<endl;
  fullrep << "ph_fish_sel ph_fish2_sel ph_ifq ph_ifq_block2 ph_fish_sel_delt ph_fish_sel_delt_alt"<<endl;
  fullrep <<  ph_fish_sel << " "<< ph_fish2_sel << " "<< ph_ifq <<" "<<  ph_ifq_block2 <<" "<<ph_fish_sel_delt<<" "<<ph_fish_sel_delt_alt<< endl;
  fullrep << "ph_srv1_sel ph_srv2_sel ph_LL_block2 ph_srv_sel_delt ph_srv_sel_delt_alt"<<endl;
  fullrep <<  ph_srv1_sel << " "<< ph_srv2_sel<<" "<<ph_LL_block2<<" "<< ph_srv_sel_delt<<" "<<ph_srv_sel_delt_alt <<endl;
  fullrep << "ph_q_srv1 ph_q_srv2 ph_q_srv3 ph_q_srv4 ph_q_srv5" <<endl;
  fullrep <<  ph_q_srv1 << " "<< ph_q_srv2 <<" "<<  ph_q_srv3<<" "<<ph_q_srv4<<" "<<ph_q_srv5  <<endl;
  fullrep << "ph_q_srv6 ph_q_srv7 ph_q_srv8 ph_srv2_q2 " <<endl;
  fullrep <<  ph_q_srv6 << " "<< ph_q_srv7 <<" "<<  ph_q_srv8<<" "<<ph_srv2_q2  <<endl<<endl<<endl;


  fullrep<<model_name<<endl;
  fullrep<<data_file<<endl;
  fullrep<<"Num_parameters_Estimated "<<initial_params::nvarcalc()<<endl;
  fullrep << "Year "<< yy <<endl;
  fullrep << "Pred_Catch Fixed Gear "<< pred_catch_fish1<<endl<<"Pred catch trawl "<<pred_catch_fish3 <<endl;
  fullrep << "Obs_Catch Fixed Gear "<< obs_catch_fish1<<endl<<" Obs_Catch Trawl "<<obs_catch_fish3 <<endl;

  fullrep << "Survival Female"<<aa <<endl;
  for (i=styr;i<=endyr;i++) fullrep << i<<" "<<S_f(i) <<endl; fullrep<<endl;
  fullrep << "Survival Male"<<aa <<endl;
  for (i=styr;i<=endyr;i++) fullrep << i<<" "<<S_m(i) <<endl; fullrep<<endl;

  fullrep << "Natural Mortality"<<aa <<endl;
  for (i=styr;i<=endyr;i++) fullrep << i<<" "<<M(i) <<endl; fullrep<<endl;
  
  fullrep << "Numbers Females "<<aa <<endl;
  for (i=styr;i<=endyr;i++) fullrep << i<<" "<<natage_f(i) <<endl; fullrep<<endl;
  fullrep << "Numbers Males"<<aa <<endl;
  for (i=styr;i<=endyr;i++) fullrep << i<<" "<<natage_m(i) <<endl; fullrep<<endl;

  fullrep << "Numbers at Length Females "<<len_bin_labels <<endl;
  for (i=styr;i<=endyr;i++) fullrep << i<<" "<<num_len_f(i) <<endl; fullrep<<endl;
  fullrep << "Numbers at Length  Males"<<len_bin_labels <<endl;
  for (i=styr;i<=endyr;i++) fullrep << i<<" "<<num_len_m(i) <<endl; fullrep<<endl;

  fullrep << "Age "<<aa <<endl;
  fullrep << "Fishery_sel1 Females"<<fish1_sel_f  <<endl;
  fullrep << "Fishery_sel1 Males"<<fish1_sel_m  <<endl;
  fullrep << "Fishery_sel2 "<<fish2_sel  <<endl;
  fullrep << "Fishery_sel3_f "<<fish3_sel_f  <<endl;
  fullrep << "Fishery_sel3_m "<<fish3_sel_m  <<endl;
  fullrep << "Fishery_sel4_f "<<fish4_sel_f  <<endl;
  fullrep << "Fishery_sel4_m "<<fish4_sel_m  <<endl;
  fullrep << "Fishery_sel5_f "<<fish5_sel_f  <<endl;
  fullrep << "Fishery_sel5_m "<<fish5_sel_m  <<endl;
  fullrep << "Survey_sel1 Female"<<srv1_sel_f  <<endl<<endl;
  fullrep << "Survey_sel1 male"<<srv1_sel_m  <<endl<<endl;
  fullrep << "Survey_sel2 Female"<<srv2_sel_f  <<endl<<endl;
  fullrep << "Survey_sel2 male"<<srv2_sel_m  <<endl<<endl;
  fullrep << "Survey_sel7 Female"<<srv7_sel_f  <<endl<<endl;
  fullrep << "Survey_sel7 male"<<srv7_sel_m  <<endl<<endl;
  fullrep << "Survey_sel10 Female"<<srv10_sel_f  <<endl<<endl;
  fullrep << "Survey_sel10 male"<<srv10_sel_m  <<endl<<endl;

  sdnr_fish1_age = 0;
  sdnr_fish1_size = 0;
  sdnr_fish3_size = 0;
  sdnr_srv1_age = 0;
  sdnr_srv2_age = 0;
  sdnr_srv1_size = 0;
  sdnr_srv2_size = 0;
  sdnr_srv7_size = 0;


  fullrep << "Obs_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) {
      sdnr_fish1_age +=sdnr(eac_fish1(i),oac_fish1(i),wt_fish1_age*double(nsamples_fish1_age(i)))/nyrs_fish1_age;
      fullrep << yrs_fish1_age(i)<<" "<<oac_fish1(i) 
      <<" eff_N "<<(1-eac_fish1(i))*eac_fish1(i)/norm2(oac_fish1(i)-eac_fish1(i))  <<" N "<<nsamples_fish1_age(i)
      <<" SDNR "<< sdnr(eac_fish1(i),oac_fish1(i),wt_fish1_age*double(nsamples_fish1_age(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_fish1_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) fullrep << yrs_fish1_age(i)<<" "<<eac_fish1(i) <<endl; fullrep<<endl;

  fullrep << "Obs_P_fish1_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++)  {
   sdnr_fish1_size += sdnr(esc_fish1_f(i),osc_fish1_f(i),wt_fish1_size*double(nsamples_fish1_size(i)))/nyrs_fish1_size/2;
  fullrep << yrs_fish1_size(i)<<" "<<osc_fish1_f(i) 
      <<" eff_N "<<(1-esc_fish1_f(i))*esc_fish1_f(i)/norm2(osc_fish1_f(i)-esc_fish1_f(i))  <<" N "<<nsamples_fish1_size(i)
      <<" SDNR "<< sdnr(esc_fish1_f(i),osc_fish1_f(i),wt_fish1_size*double(nsamples_fish1_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_fish1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) fullrep << yrs_fish1_size(i)<<" "<<esc_fish1_f(i) <<endl; fullrep<<endl;

    fullrep << "Obs_P_fish1_size Male"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) {
    sdnr_fish1_size += sdnr(esc_fish1_m(i),osc_fish1_m(i),wt_fish1_size*double(nsamples_fish1_size(i)))/nyrs_fish1_size/2;
  fullrep << yrs_fish1_size(i)<<" "<<osc_fish1_m(i) 
      <<" eff_N "<<(1-esc_fish1_m(i))*esc_fish1_m(i)/norm2(osc_fish1_m(i)-esc_fish1_m(i))  <<" N "<<nsamples_fish1_size(i)
      <<" SDNR "<< sdnr(esc_fish1_m(i),osc_fish1_m(i),wt_fish1_size*double(nsamples_fish1_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_fish1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) fullrep << yrs_fish1_size(i)<<" "<<esc_fish1_m(i) <<endl; fullrep<<endl;

   fullrep << "Obs_P_fish3_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) {
    sdnr_fish3_size+= sdnr(esc_fish3_m(i),osc_fish3_m(i),wt_fish3_size*double(nsamples_fish3_size(i)))/nyrs_fish3_size/2;
  fullrep << yrs_fish3_size(i)<<" "<<osc_fish3_m(i) 
      <<" eff_N "<<(1-esc_fish3_m(i))*esc_fish3_m(i)/norm2(osc_fish3_m(i)-esc_fish3_m(i))  <<" N "<<nsamples_fish3_size(i)
      <<" SDNR "<< sdnr(esc_fish3_m(i),osc_fish3_m(i),wt_fish3_size*double(nsamples_fish3_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_fish3_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) fullrep << yrs_fish3_size(i)<<" "<<esc_fish3_m(i) <<endl; fullrep<<endl;

  fullrep << "Obs_P_fish3_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++)  {
     sdnr_fish3_size+= sdnr(esc_fish3_f(i),osc_fish3_f(i),wt_fish3_size*double(nsamples_fish3_size(i)))/nyrs_fish3_size/2;
  fullrep << yrs_fish3_size(i)<<" "<<osc_fish3_f(i) 
      <<" eff_N "<<(1-esc_fish3_f(i))*esc_fish3_f(i)/norm2(osc_fish3_f(i)-esc_fish3_f(i))  <<" N "<<nsamples_fish3_size(i)
      <<" SDNR "<< sdnr(esc_fish3_f(i),osc_fish3_f(i),wt_fish3_size*double(nsamples_fish3_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_fish3_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) fullrep << yrs_fish3_size(i)<<" "<<esc_fish3_f(i) <<endl; fullrep<<endl;

  fullrep << "Obs_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) {
   sdnr_srv1_age+=sdnr(eac_srv1(i),oac_srv1(i),wt_srv1_age*double(nsamples_srv1_age(i)))/nyrs_srv1_age;
   fullrep << yrs_srv1_age(i)<<" "<<oac_srv1(i) 
      <<" eff_N "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i)) <<" N "<<nsamples_srv1_age(i)
      <<" SDNR "<< sdnr(eac_srv1(i),oac_srv1(i),wt_srv1_age*double(nsamples_srv1_age(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) fullrep << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; fullrep<<endl;

     fullrep << "Obs_P_srv2_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv2_age;i++) {
    sdnr_srv2_age+=sdnr(eac_srv2(i),oac_srv2(i),wt_srv1_age*double(nsamples_srv2_age(i)))/nyrs_srv2_age;
    fullrep << yrs_srv2_age(i)<<" "<<oac_srv2(i) 
      <<" eff_N "<<(1-eac_srv2(i))*eac_srv2(i)/norm2(oac_srv2(i)-eac_srv2(i)) <<" N "<<nsamples_srv2_age(i)
      <<" SDNR "<< sdnr(eac_srv2(i),oac_srv2(i),wt_srv1_age*double(nsamples_srv2_age(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv2_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv2_age;i++) fullrep << yrs_srv2_age(i)<<" "<<eac_srv2(i) <<endl; fullrep<<endl;

  fullrep << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) {
     sdnr_srv1_size+=sdnr(esc_srv1_f(i),osc_srv1_f(i),wt_srv1_size*double(nsamples_srv1_size(i)))/nyrs_srv1_size/2;
    fullrep << yrs_srv1_size(i)<<" "<<osc_srv1_f(i) 
      <<" eff_N "<<(1-esc_srv1_f(i))*esc_srv1_f(i)/norm2(osc_srv1_f(i)-esc_srv1_f(i)) <<" N "<<nsamples_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1_f(i),osc_srv1_f(i),wt_srv1_size*double(nsamples_srv1_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) fullrep << yrs_srv1_size(i)<<" "<<esc_srv1_f(i) <<endl; fullrep<<endl;

    fullrep << "Obs_P_srv1_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) {
    sdnr_srv1_size+=sdnr(esc_srv1_m(i),osc_srv1_m(i),wt_srv1_size*double(nsamples_srv1_size(i)))/nyrs_srv1_size/2;
     fullrep << yrs_srv1_size(i)<<" "<<osc_srv1_m(i) 
      <<" eff_N "<<(1-esc_srv1_m(i))*esc_srv1_m(i)/norm2(osc_srv1_m(i)-esc_srv1_m(i)) <<" N "<<nsamples_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1_m(i),osc_srv1_m(i),wt_srv1_size*double(nsamples_srv1_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) fullrep << yrs_srv1_size(i)<<" "<<esc_srv1_m(i) <<endl; fullrep<<endl;

  fullrep << "Obs_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) {
    sdnr_srv2_size+=sdnr(esc_srv2_f(i),osc_srv2_f(i),wt_srv2_size*double(nsamples_srv2_size(i)))/nyrs_srv2_size/2;
     fullrep << yrs_srv2_size(i)<<" "<<osc_srv2_f(i) 
      <<" eff_N "<<(1-esc_srv2_f(i))*esc_srv2_f(i)/norm2(osc_srv2_f(i)-esc_srv2_f(i)) <<" N "<<nsamples_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2_f(i),osc_srv2_f(i),wt_srv2_size*double(nsamples_srv2_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) fullrep << yrs_srv2_size(i)<<" "<<esc_srv2_f(i) <<endl; fullrep<<endl;

    fullrep << "Obs_P_srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) {
        sdnr_srv2_size+=sdnr(esc_srv2_m(i),osc_srv2_m(i),wt_srv2_size*double(nsamples_srv2_size(i)))/nyrs_srv2_size/2;
 fullrep << yrs_srv2_size(i)<<" "<<osc_srv2_m(i) 
      <<" eff_N "<<(1-esc_srv2_m(i))*esc_srv2_m(i)/norm2(osc_srv2_m(i)-esc_srv2_m(i)) <<" N "<<nsamples_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2_m(i),osc_srv2_m(i),wt_srv2_size*double(nsamples_srv2_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) fullrep << yrs_srv2_size(i)<<" "<<esc_srv2_m(i) <<endl; fullrep<<endl;
  
    fullrep << "Obs_P_srv7_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) {
    sdnr_srv7_size+= sdnr(esc_srv7_f(i),osc_srv7_f(i),wt_srv7_size*double(nsamples_srv7_size(i)))/nyrs_srv7_size/2;
     fullrep << yrs_srv7_size(i)<<" "<<osc_srv7_f(i) 
      <<" eff_N "<<(1-esc_srv7_f(i))*esc_srv7_f(i)/norm2(osc_srv7_f(i)-esc_srv7_f(i)) <<" N "<<nsamples_srv7_size(i)
      <<" SDNR "<< sdnr(esc_srv7_f(i),osc_srv7_f(i),wt_srv7_size*double(nsamples_srv7_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv7_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) fullrep << yrs_srv7_size(i)<<" "<<esc_srv7_f(i) <<endl; fullrep<<endl;

    fullrep << "Obs_P_srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) {
        sdnr_srv7_size+= sdnr(esc_srv7_m(i),osc_srv7_m(i),wt_srv7_size*double(nsamples_srv7_size(i)))/nyrs_srv7_size/2;
  fullrep << yrs_srv7_size(i)<<" "<<osc_srv7_m(i) 
      <<" eff_N "<<(1-esc_srv7_m(i))*esc_srv7_m(i)/norm2(osc_srv7_m(i)-esc_srv7_m(i)) <<" N "<<nsamples_srv7_size(i)
      <<" SDNR "<< sdnr(esc_srv7_m(i),osc_srv7_m(i),wt_srv7_size*double(nsamples_srv7_size(i)))<<endl; fullrep<<endl; }
  fullrep << "Pred_P_srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) fullrep << yrs_srv7_size(i)<<" "<<esc_srv7_m(i) <<endl; fullrep<<endl;

  fullrep << "Survey Biomass " <<endl;
  fullrep << "Year:     " << yrs_srv1  <<endl;
  fullrep << "Predicted:   " << pred_srv1  <<endl;
  fullrep << "Observed:   " << obs_srv1_biom  <<endl<< endl;

  fullrep << "Weight "<<aa<< endl;
  fullrep << " "<< weight_f << endl;

  fullrep << "Weight Females "<<aa<< endl;
  fullrep << " "<< weight_f << endl;

  fullrep << "Weight Males "<<aa<< endl;
  fullrep << " "<< weight_m << endl;

  fullrep << "Fully_selected_F "<<aa <<endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     fullrep << " "<< Fmort_fish1*max(fish4_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     fullrep << " "<< Fmort_fish1*max(fish5_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    } 
   if(ph_ifq<1) // no post-IFQ sel est
    {
     fullrep << " "<< Fmort_fish1*max(fish1_sel_f)+Fmort_fish3*max(fish3_sel_f) <<endl;
    }
  fullrep << "Year " << yy<< endl;
  fullrep << "SpBiom "<< spawn_biom <<endl;
  fullrep << "Tot_biom "<< tot_biom   <<endl;
    // fullrep << tot_biom.sd <<endl;

  fullrep << "Fishery Selectivity " <<endl;
   if((ph_ifq>0) && (ph_ifq_block2<1))  // post-IFQ sel does not  includes a recent time block
    {
     fullrep << fish4_sel_f / max(fish4_sel_f) <<endl;
    }
   if((ph_ifq>0) && (ph_ifq_block2>0)) // post-IFQ sel includes a recent time block
    {
     fullrep << fish5_sel_f / max(fish5_sel_f) <<endl;
    }
   if(ph_ifq<1) // no post-IFQ sel est
    {
     fullrep << fish1_sel_f / max(fish1_sel_f) <<endl;
    }
  fullrep << "F35 F40 F50 "<<endl;
  fullrep <<  F35 << " "<< mF40 <<" "<<  F50 <<endl <<endl <<endl <<endl <<endl <<endl;
  fullrep << "SSB projection: "<< spawn_biom_proj<<endl<<"Catch projection: "<<pred_catch_proj<<endl;
  fullrep <<  "B40: "<< B40<< endl;

  fullrep << "Wts_n_Likelihoods  " << endl;
  fullrep << wt_ssqcatch_fish1 <<" "<<ssqcatch <<" " ; fullrep << "SSQ_Catch_Likelihood" << endl;
  fullrep << wt_srv3     <<" "<<surv_like(3)  <<" " ; fullrep << "Domestic_Survey_Abundance_Index_Likelihood" << endl;
  fullrep << wt_srv4     <<" "<<surv_like(4)  <<" " ; fullrep << "Cooperative_Survey_Abundance_Index_Likelihood" << endl;
  fullrep << wt_fish1_age <<" "<<age_like(1)  <<" " ; fullrep << "Fishery_Age_Composition_Likelihood" << endl;
  fullrep << wt_srv1_age <<" "<<age_like(2)  <<" " ; fullrep << "Survey_Age_Composition_Likelihood_DomesticLL" << endl;
  fullrep << wt_srv1_age <<" "<<age_like(16)  <<" " ; fullrep << "Survey_Age_Composition_Likelihood_CooperativeLL" << endl;
  fullrep << wt_fish1_size<<" "<<age_like(3)+age_like(4)  <<" " ; fullrep << "Fishery_Size_Composition_Likelihood_Fixed" << endl;
  fullrep << wt_fish3_size<<" "<<age_like(6)+age_like(7)  <<" " ; fullrep << "Fishery_Size_Composition_Likelihood_Trawl" << endl;
  fullrep << wt_srv1_size<<" "<<age_like(9)+age_like(10)  <<" " ; fullrep << "Survey_Size_Composition_Likelihood_Domestic LL" << endl;
  fullrep << wt_srv2_size<<" "<<age_like(11)+age_like(12)  <<" " ; fullrep << "Survey_Size_Composition_Likelihood_Cooperative_LL" << endl;
  fullrep << wt_srv7_size<<" "<<age_like(13)+age_like(14)  <<" " ; fullrep << "Survey_Size_Composition_Likelihood_GOATrawl" << endl;
  fullrep << wt_rec_var <<" "<<rec_like     <<" " ; fullrep << "Recruitment_Deviations_Likelihood" << endl;

  fullrep << wt_sel_reg_fish1 <<" "<<sel_like(1)      <<" " ; fullrep << "Fish_sel_Regularity_Penalty  "<<endl  ;
  fullrep << wt_sel_reg_srv1 <<" "<<sel_like(2)      <<" " ; fullrep << "Surv_sel_Regularity_Penalty  "<<endl  ;
  fullrep << wt_sel_dome_fish1<<" "<<sel_like(3)      <<" " ; fullrep << "Fish_Sel_Domeshapedness_Penalty "<<endl  ;
  fullrep << wt_sel_dome_srv1<<" "<<sel_like(4)      <<" " ; fullrep << "Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  fullrep << "0"   <<" "<<avg_sel_penalty  <<" " ; fullrep << "Average_Selectivity_Condition_Penalty "<<endl  ;
  fullrep << wt_fmort_reg     <<" "<<F_mort_regularity<<" " ; fullrep << "Fishing_Mortality_Regularity_Penalty" << endl;
  fullrep << wt_M_reg     <<" "<<M_mort_regularity<<" " ; fullrep << "Natural_Mortality_Regularity_Penalty" << endl;
 // fullrep << wt_ssqcatch     <<" "<<F_dev_penalty<<" " ; fullrep << "Fishing_Mortality_Deviations_Penalty" << endl;
  fullrep << " "<<priors(1)  <<" " ; fullrep << "priors sigr"     <<endl;
  fullrep << " "<<priors(2)  <<" " ; fullrep << "priors q_1" <<endl;
   fullrep << " "<<priors(5)  <<" " ; fullrep << "priors q_2" <<endl;
  fullrep << " "<<priors(6)  <<" " ; fullrep << "priors q_3" <<endl;
  fullrep << " "<<priors(7)  <<" " ; fullrep << "priors q_4" <<endl;
  fullrep << " "<<priors(8)  <<" " ; fullrep << "priors q_5" <<endl;
  fullrep << " "<<priors(9)  <<" " ; fullrep << "priors q_6" <<endl;
  fullrep << " "<<priors(10)  <<" " ; fullrep << "priors q_7" <<endl;
  fullrep << " "<<priors(10)  <<" " ; fullrep << "priors q_8" <<endl;
  fullrep << " "<<priors(4)  <<" " ; fullrep << "priors M"<<endl;
  fullrep << " "<<obj_fun    <<" " ; fullrep << "obj_fun"         <<endl;
  fullrep << " "<<Like       <<" " ; fullrep << "data likelihood" <<endl;//(2*square(sigr))+ size_count(log_rec_dev)*log(sigr)<<endl;
 if(last_phase()) { 
  fullrep <<" SDNR1 "<< wt_srv1*std_dev(elem_div((pred_srv1(yrs_srv1)-obs_srv1_biom),obs_srv1_se))<<endl;
  fullrep <<" SDNR2 "<< wt_srv2*std_dev(elem_div((pred_srv2(yrs_srv2)-obs_srv2_biom),obs_srv2_se))<<endl;
  fullrep <<" SDNR3 "<< wt_srv3*std_dev(elem_div((pred_srv3(yrs_srv3)-obs_srv3_biom),obs_srv3_se))<<endl;
  fullrep <<" SDNR4 "<< wt_srv4*std_dev(elem_div((pred_srv4(yrs_srv4)-obs_srv4_biom),obs_srv4_se))<<endl;
  fullrep <<" SDNR5 "<< wt_srv5*std_dev(elem_div((pred_srv5(yrs_srv5)-obs_srv5_biom),obs_srv5_se))<<endl;
  fullrep <<" SDNR6 "<< wt_srv6*std_dev(elem_div((pred_srv6(yrs_srv6)-obs_srv6_biom),obs_srv6_se))<<endl;
  fullrep <<" SDNR7 "<< wt_srv7*std_dev(elem_div((pred_srv7(yrs_srv7)-obs_srv7_biom),obs_srv7_se))<<endl;
//  fullrep <<" SDNR8 "<< wt_srv8*std_dev(elem_div((pred_srv8(yrs_srv8)-obs_srv8_biom),obs_srv8_se))<<endl;
 

    }
  fullrep << "SigmaR: "<<sigr<< " Nat_Mort_Base: "<<natmort<<" Male delta M "<<mdelta<<" Spawning Per Recruit "<< " "<<SBF40<<" "<<SB0<<" Virgin SPR "<<endl;
  fullrep << "Stock-recruitment, type: "<<SrType<<" 1=Ricker, 2=B-Holt, 3=Mean"<<endl;
  fullrep << "Year SSB SR_Pred R_Est "<<endl;
  fullrep << "Weighted likelihods for comps broken down"<<endl;
  fullrep << "Fish1 Age, Survey 1 Age, Fish 1 Size Female, Fish 1 Size Male, Fish 2 Size, Fish3 Size Female, F3 Size Male, Srv1 Size F, Srv1 Size M, Srv2 Size F, Srv2 Size M, Srv7 Size F, Srv7 Size M"<<endl;
  fullrep <<age_like<<" age _like"<<endl<<surv_like<<" surv like "<<endl<<sel_like<<" sel_like"<<endl;
  fullrep << "Unweighted likelihods for comps broken down"<<endl;
  fullrep << "Fish1 Age, Survey 1 Age, Fish 1 Size Female, Fish 1 Size Male, Fish 2 Size, Fish3 Size Female, F3 Size Male, Srv1 Size F, Srv1 Size M, Srv2 Size F, Srv2 Size M, Srv7 Size F, Srv7 Size M"<<endl;
  fullrep <<age_like(1) / (0.00001+wt_fish1_age)<< " "<<(age_like(2)+0.0001) / (0.00001+wt_srv1_age)<< " "<<(age_like(3)+0.0001) / (0.00001+wt_fish1_size)<< " "<<(age_like(4)+0.0001) / (0.00001+wt_fish1_size)<< " "<<(age_like(5)+0.0001) / (0.00001+wt_fish2_size)<< " "<<(age_like(6)+0.0001) / (0.00001+wt_fish3_size)<< " "<<(age_like(7)+0.0001) / (0.00001+wt_fish3_size)<< " "<<(age_like(8)+0.0001) / (0.00001+wt_fish4_size)<< " "<<(age_like(9)+0.0001) / (0.00001+wt_srv1_size)<< " "<<(age_like(10)+0.0001) / (0.00001+wt_srv1_size)<< " "<<(age_like(11)+0.0001) / (0.00001+wt_srv2_size)<< " "<<(age_like(12)+0.0001) / (0.00001+wt_srv2_size)<< " "<<(age_like(13)+0.0001) / (0.00001+wt_srv7_size)<< " "<<(age_like(14)+0.0001) / (0.00001+wt_srv7_size)<<" "<<(age_like(15)+0.0001) / (0.00001+wt_srv7_age)<< " "<<(age_like(16)+0.0001) / (0.00001+wt_srv1_age) <<" unw_age)_like"<<endl;
  fullrep <<surv_like(1) / (0.00001+wt_srv1)<< " "<<surv_like(2) / (0.00001+wt_srv2)<< " "<<surv_like(3) / (0.00001+wt_srv3)<< " "<<surv_like(4) / (0.00001+wt_srv4)<< " "<<surv_like(5) / (0.00001+wt_srv5)<< " "<<surv_like(6) / (0.00001+wt_srv6)<< " "<<(surv_like(7)+0.0001) / (0.00001+wt_srv7)<< " "<<(surv_like(8)+0.0001) / (0.00001+wt_srv8)<< " "<<"unw_surv_like "<<endl<<sel_like<<" sel_like"<<endl;
  fullrep << "Survey Biomass 2" <<endl;
  fullrep << "Year:     " << yrs_srv2  <<endl;
  fullrep << "Predicted:   " << pred_srv2  <<endl;
  fullrep << "Observed:   " << obs_srv2_biom  <<endl<< endl;
  fullrep << "Survey Biomass 3" <<endl;
  fullrep << "Year:     " << yrs_srv3  <<endl;
  fullrep << "Predicted:   " << pred_srv3  <<endl;
  fullrep << "Observed:   " << obs_srv3_biom  <<endl<< endl;

  fullrep << "Survey Biomass 4" <<endl;
  fullrep << "Year:     " << yrs_srv4  <<endl;
  fullrep << "Predicted:   " << pred_srv4  <<endl;
  fullrep << "Observed:   " << obs_srv4_biom  <<endl<< endl;

  fullrep << "Survey Biomass 5" <<endl;
  fullrep << "Year:     " << yrs_srv5  <<endl;
  fullrep << "Predicted:   " << pred_srv5  <<endl;
  fullrep << "Observed:   " << obs_srv5_biom  <<endl<< endl;

  fullrep << "Survey Biomass 6" <<endl;
  fullrep << "Year:     " << yrs_srv6  <<endl;
  fullrep << "Predicted:   " << pred_srv6  <<endl;
  fullrep << "Observed:   " << obs_srv6_biom  <<endl<< endl;
  fullrep << "Survey Biomass 7" <<endl;
  fullrep << "Year:     " << yrs_srv7  <<endl;
  fullrep << "Predicted:   " << pred_srv7  <<endl;
  fullrep << "Observed:   " << obs_srv7_biom  <<endl<< endl;
  fullrep << "qs"<<endl<<q_srv1<<endl<<q_srv2<<endl<<q_srv3<<endl<<q_srv4<<endl<<q_srv5<<endl<<q_srv6<<endl<<q_srv7<<endl<<q_srv8<<endl;
  fullrep << "Year SSB SRR Recr "<<endl;
  for (i=styr;i<=endyr;i++)
  fullrep<< i <<" "<<Sp_Biom(i-recage)<<" "<<srm_rec(i)<<" "<<sam_rec(i)<<endl;
  fullrep<< "Age/Length residuals"<<endl;
  fullrep << "fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish1_age;i++) fullrep << yrs_fish1_age(i)<<" "<<oac_fish1(i)-eac_fish1(i)<<endl; 

  fullrep << "fish1_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) fullrep << yrs_fish1_size(i)<<" "<<osc_fish1_f(i)-esc_fish1_f(i)<<endl; 

  fullrep << "fish1_size Male"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish1_size;i++) fullrep << yrs_fish1_size(i)<<" "<<osc_fish1_m(i)-esc_fish1_m(i)<<endl; 

  fullrep << "fish3_size Female"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) fullrep << yrs_fish3_size(i)<<" "<<osc_fish3_f(i)-esc_fish3_f(i)<<endl;
   
  fullrep << "fish3_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish3_size;i++) fullrep << yrs_fish3_size(i)<<" "<<osc_fish3_m(i)-esc_fish3_m(i)<<endl;

  fullrep << "srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) fullrep << yrs_srv1_age(i)<<" "<<oac_srv1(i)-eac_srv1(i)<<endl; 

  fullrep << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) fullrep << yrs_srv1_size(i)<<" "<<osc_srv1_f(i)-esc_srv1_f(i)<<endl;

  fullrep << "srv1_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) fullrep << yrs_srv1_size(i)<<" "<<osc_srv1_m(i)-esc_srv1_m(i)<<endl;
  
  fullrep << "srv2_size Females"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) fullrep << yrs_srv2_size(i)<<" "<<osc_srv2_f(i)-esc_srv2_f(i)<<endl;

  fullrep << "srv2_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) fullrep << yrs_srv2_size(i)<<" "<<osc_srv2_m(i)-esc_srv2_m(i)<<endl;

    fullrep << "srv7_size Females"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) fullrep << yrs_srv7_size(i)<<" "<<osc_srv7_f(i)-esc_srv7_f(i)<<endl;

  fullrep << "srv7_size Males"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv7_size;i++) fullrep << yrs_srv7_size(i)<<" "<<osc_srv7_m(i)-esc_srv7_m(i)<<endl;

  fullrep << "N_proj_f "<<endl<<N_proj_f<<endl<<"N_proj_m "<<N_proj_m<<endl<<" spawn_bio next year"<<endl<<spawn_biom_proj(endyr+1)<<endl
  <<"ABC full"<<ABC3<<endl;
  
  fullrep <<" spawn_bio projected"<<endl<<spawn_biom_proj<<endl;
  fullrep << "Specified catch projection"<<(catage_proj_f*weight_f(endyr)+catage_proj_m*weight_m(endyr))<<endl<<"ABC projection: "<<pred_catch_proj<<endl;
  // New SDNR stuff
  fullrep<<" average SDNRs for compositional data"<<endl;
  fullrep<<"Fishery Ages "<<sdnr_fish1_age<<endl;
  fullrep<<"Fishery 1 sizes "<<sdnr_fish1_size<<endl;
  fullrep<<"Fishery 3 sizes "<<sdnr_fish3_size<<endl;
  fullrep<<"Survey 1 ages "<<sdnr_srv1_age<<endl;
  fullrep<<"Survey 2 ages "<<sdnr_srv2_age<<endl;
  fullrep<<"Survey 1 sizes "<<sdnr_srv1_size<<endl;
  fullrep<<"Survey 2 sizes "<<sdnr_srv2_size<<endl;
  fullrep<<"Survey 7 sizes "<<sdnr_srv7_size<<endl;

    fullrep<<"weight_f "<<weight_f<<endl;
    fullrep<<"weight_m "<<weight_m<<endl;
    fullrep<<"size_age_f "<<size_age_f<<endl;
    fullrep<<"size_age_m "<<size_age_m<<endl;
    fullrep<<"maturity "<<maturity<<endl;

RUNTIME_SECTION
  convergence_criteria 1.e-7  
  maximum_function_evaluations 20000

TOP_OF_MAIN_SECTION
 gradient_structure::set_MAX_NVAR_OFFSET(1000);
   gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
   gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
   gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  arrmblsize=3900000;

FINAL_SECTION
 write_fullrep();
