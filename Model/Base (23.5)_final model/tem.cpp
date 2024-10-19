#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"
  #include "statsLib.h"
  #include "qfclib.h"
  #include "mhp-s-funcs.cpp"                // Include S-compatible output functions (needs preceding)
  #include <contrib.h>
  #define EOUT(var) cout <<#var<<" "<<var<<endl;
  adstring model_name;
  adstring data_file;
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <tem.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("evalout.prj");;
 ad_comm::change_datafile_name("tem.ctl");    //Read in phases, penalties and priors from "tem.ctl"
 *(ad_comm::global_datafile) >>  model_name;  //define model name from .ctl file
 *(ad_comm::global_datafile) >>  data_file;   //define .dat file from .ctl file (tem_2020_na_wh.dat...replace 2020 with current year)
  data_reweight_switch.allocate("data_reweight_switch");
  growth_blocks.allocate("growth_blocks");
  growth_cutoffs.allocate(1,growth_blocks,"growth_cutoffs");
  weight_blocks.allocate("weight_blocks");
  weight_cutoffs.allocate(1,weight_blocks,"weight_cutoffs");
  maturity_blocks.allocate("maturity_blocks");
  maturity_cutoffs.allocate(1,maturity_blocks,"maturity_cutoffs");
  SrType.allocate("SrType");
  styr_rec_est.allocate("styr_rec_est");
  endyr_rec_est.allocate("endyr_rec_est");
 nrecs_est = endyr_rec_est-styr_rec_est+1;
  rec_like_type.allocate("rec_like_type");
  bias_ramp.allocate("bias_ramp");
  bmax.allocate("bmax");
  b_a50.allocate("b_a50");
  b_year_st.allocate("b_year_st");
  b_year_end.allocate("b_year_end");
  sigma_R_early_switch.allocate("sigma_R_early_switch");
  sigma_R_early.allocate("sigma_R_early");
  sigma_R_early_end.allocate("sigma_R_early_end");
  ph_mean_rec.allocate("ph_mean_rec");
  ph_recdev.allocate("ph_recdev");
  ph_steepness.allocate("ph_steepness");
  ph_Rzero.allocate("ph_Rzero");
  ph_sigr.allocate("ph_sigr");
  ph_m.allocate("ph_m");
  ph_mdelta.allocate("ph_mdelta");
  ph_Mdevs.allocate("ph_Mdevs");
  ph_Mdevs_age.allocate("ph_Mdevs_age");
  ph_avg_F.allocate("ph_avg_F");
  ph_Fdev.allocate("ph_Fdev");
  ph_F50.allocate("ph_F50");
  ph_fish_sel.allocate("ph_fish_sel");
  ph_fish2_sel.allocate("ph_fish2_sel");
  ph_fish_sel_delt.allocate("ph_fish_sel_delt");
  ph_fish_sel_delt_alt.allocate("ph_fish_sel_delt_alt");
  ph_srv1_sel.allocate("ph_srv1_sel");
  ph_srv2_sel.allocate("ph_srv2_sel");
  ph_srv_sel_delt.allocate("ph_srv_sel_delt");
  ph_srv_sel_delt_alt.allocate("ph_srv_sel_delt_alt");
  ph_ifq.allocate("ph_ifq");
  ph_ifq_block2.allocate("ph_ifq_block2");
  ph_LL_block2.allocate("ph_LL_block2");
  yr_sel_chg_fish.allocate("yr_sel_chg_fish");
  yr_sel_chg_srv1.allocate("yr_sel_chg_srv1");
  ph_q_srv1.allocate("ph_q_srv1");
  ph_q_srv2.allocate("ph_q_srv2");
  ph_q_srv3.allocate("ph_q_srv3");
  ph_q_srv4.allocate("ph_q_srv4");
  ph_q_srv5.allocate("ph_q_srv5");
  ph_q_srv6.allocate("ph_q_srv6");
  ph_q_srv7.allocate("ph_q_srv7");
  ph_q_srv8.allocate("ph_q_srv8");
  ph_srv2_q2.allocate("ph_srv2_q2");
  ph_q_LL_srv_rec.allocate("ph_q_LL_srv_rec");
  ph_q_IFQ_rec.allocate("ph_q_IFQ_rec");
  fsh1_sel_opt.allocate("fsh1_sel_opt");
  fsh2_sel_opt.allocate("fsh2_sel_opt");
  fsh3_sel_opt.allocate("fsh3_sel_opt");
  fsh4_sel_opt.allocate("fsh4_sel_opt");
  fsh5_sel_opt.allocate("fsh5_sel_opt");
  srv1_sel_opt.allocate("srv1_sel_opt");
  srv2_sel_opt.allocate("srv2_sel_opt");
  srv7_sel_opt.allocate("srv7_sel_opt");
  srv10_sel_opt.allocate("srv10_sel_opt");
  n_fish_sel_ages.allocate("n_fish_sel_ages");
  n_srv1_sel_ages.allocate("n_srv1_sel_ages");
  mprior.allocate("mprior");
 log_mprior = log(mprior);
  cvmprior.allocate("cvmprior");
  steep_prior.allocate("steep_prior");
  cv_steep_prior.allocate("cv_steep_prior");
  sigrprior.allocate("sigrprior");
  cvsigrprior.allocate("cvsigrprior");
  q_srv1prior.allocate("q_srv1prior");
  cvq_srv1prior.allocate("cvq_srv1prior");
 log_q1prior = log(q_srv1prior);
  q_srv2prior.allocate("q_srv2prior");
  cvq_srv2prior.allocate("cvq_srv2prior");
 log_q2prior = log(q_srv2prior);
  q_srv3prior.allocate("q_srv3prior");
  cvq_srv3prior.allocate("cvq_srv3prior");
 log_q3prior = log(q_srv3prior);
  q_srv4prior.allocate("q_srv4prior");
  cvq_srv4prior.allocate("cvq_srv4prior");
 log_q4prior = log(q_srv4prior);
  q_srv5prior.allocate("q_srv5prior");
  cvq_srv5prior.allocate("cvq_srv5prior");
 log_q5prior = log(q_srv5prior);
  q_srv6prior.allocate("q_srv6prior");
  cvq_srv6prior.allocate("cvq_srv6prior");
 log_q6prior = log(q_srv6prior);
  q_srv7prior.allocate("q_srv7prior");
  cvq_srv7prior.allocate("cvq_srv7prior");
 log_q7prior = log(q_srv7prior);
  q_srv8prior.allocate("q_srv8prior");
  cvq_srv8prior.allocate("cvq_srv8prior");
 log_q8prior = log(q_srv8prior);
  yr_catchwt.allocate("yr_catchwt");
  wt_ssqcatch_fish1.allocate("wt_ssqcatch_fish1");
  wt_ssqcatch_fish3.allocate("wt_ssqcatch_fish3");
  wt_srv1.allocate("wt_srv1");
  wt_srv2.allocate("wt_srv2");
  wt_srv3.allocate("wt_srv3");
  wt_srv4.allocate("wt_srv4");
  wt_srv5.allocate("wt_srv5");
  wt_srv6.allocate("wt_srv6");
  wt_srv7.allocate("wt_srv7");
  wt_srv8.allocate("wt_srv8");
  wt_fish1_age.allocate("wt_fish1_age");
  wt_srv1_age.allocate("wt_srv1_age");
  wt_fish1_size.allocate("wt_fish1_size");
  wt_srv1_size.allocate("wt_srv1_size");
  wt_fish2_size.allocate("wt_fish2_size");
  wt_srv2_size.allocate("wt_srv2_size");
  wt_fish3_size.allocate("wt_fish3_size");
  wt_srv7_size.allocate("wt_srv7_size");
  wt_fish4_size.allocate("wt_fish4_size");
  wt_srv7_age.allocate("wt_srv7_age");
  wt_srv_extra_size.allocate("wt_srv_extra_size");
  wt_srv5_size.allocate("wt_srv5_size");
  wt_fish6_size.allocate("wt_fish6_size");
  wt_srv6_size.allocate("wt_srv6_size");
  wt_rec_var.allocate("wt_rec_var");
  wt_fish1_age_iter.allocate("wt_fish1_age_iter");
  wt_srv1_age_iter.allocate("wt_srv1_age_iter");
  wt_srv2_age_iter.allocate("wt_srv2_age_iter");
  wt_fish1_size_male_iter.allocate("wt_fish1_size_male_iter");
  wt_fish1_size_female_iter.allocate("wt_fish1_size_female_iter");
  wt_srv1_size_male_iter.allocate("wt_srv1_size_male_iter");
  wt_srv1_size_female_iter.allocate("wt_srv1_size_female_iter");
  wt_srv2_size_male_iter.allocate("wt_srv2_size_male_iter");
  wt_srv2_size_female_iter.allocate("wt_srv2_size_female_iter");
  wt_fish3_size_male_iter.allocate("wt_fish3_size_male_iter");
  wt_fish3_size_female_iter.allocate("wt_fish3_size_female_iter");
  wt_srv7_size_male_iter.allocate("wt_srv7_size_male_iter");
  wt_srv7_size_female_iter.allocate("wt_srv7_size_female_iter");
  wt_sel_reg_fish1.allocate("wt_sel_reg_fish1");
  wt_sel_reg_fish2.allocate("wt_sel_reg_fish2");
  wt_sel_reg_fish3.allocate("wt_sel_reg_fish3");
  wt_sel_reg_fish4.allocate("wt_sel_reg_fish4");
  wt_sel_dome_fish1.allocate("wt_sel_dome_fish1");
  wt_sel_dome_fish2.allocate("wt_sel_dome_fish2");
  wt_sel_dome_fish3.allocate("wt_sel_dome_fish3");
  wt_sel_dome_fish4.allocate("wt_sel_dome_fish4");
  wt_sel_reg_srv1.allocate("wt_sel_reg_srv1");
  wt_sel_reg_srv2.allocate("wt_sel_reg_srv2");
  wt_sel_dome_srv1.allocate("wt_sel_dome_srv1");
  wt_sel_dome_srv2.allocate("wt_sel_dome_srv2");
  wt_fmort_reg.allocate("wt_fmort_reg");
  wt_M_reg.allocate("wt_M_reg");
  wt_q_priors.allocate("wt_q_priors");
  wt_M_prior.allocate("wt_M_prior");
  wt_sigr_prior.allocate("wt_sigr_prior");
  hist_hal_prop.allocate("hist_hal_prop");
  yieldratio.allocate("yieldratio");
 cout<<yieldratio<<endl;
 cout<<"done with ctl"<<endl;
 ad_comm::change_datafile_name(data_file);    // Read data from the data file (i.e., tem_2020_na_wh.dat)
  styr.allocate("styr");
  endyr.allocate("endyr");
  recage.allocate("recage");
  nages.allocate("nages");
  nlenbins.allocate("nlenbins");
  len_bin_labels.allocate(1,nlenbins,"len_bin_labels");
  nyrs = endyr - styr + 1;
  styr_rec = (styr - nages) + 1;     // First year of recruitment
 styr_sp  = styr_rec - recage ;     // First year of spawning biomass  
 endyr_sp = endyr   - recage - 1;   // endyr year of (main) spawning biomass
 endyr_rec= endyr_rec_est;  // 
  yy.allocate(styr,endyr);
 yy.fill_seqadd(styr,1) ;
  aa.allocate(1,nages);
 aa.fill_seqadd(recage,1) ;
  spawn_fract.allocate("spawn_fract");
 spawn_fract = (spawn_fract - 1) / 12;
  p_mature1.allocate(1,nages,"p_mature1");
  p_mature2.allocate(1,nages,"p_mature2");
  p_mature3.allocate(1,nages,"p_mature3");
  wt_m1.allocate(1,nages,"wt_m1");
  wt_f1.allocate(1,nages,"wt_f1");
  wt_m2.allocate(1,nages,"wt_m2");
  wt_f2.allocate(1,nages,"wt_f2");
  wt_all.allocate(1,nages,"wt_all");
  wt_old.allocate(1,nages,"wt_old");
  prop_m.allocate(styr,endyr,"prop_m");
  prop_m2.allocate(styr,endyr,"prop_m2");
  obs_catch_fish1.allocate(styr,endyr,"obs_catch_fish1");
  obs_catch_fish3.allocate(styr,endyr,"obs_catch_fish3");
  nyrs_srv1.allocate("nyrs_srv1");
  yrs_srv1.allocate(1,nyrs_srv1,"yrs_srv1");
  obs_srv1_biom.allocate(1,nyrs_srv1,"obs_srv1_biom");
  obs_srv1_se.allocate(1,nyrs_srv1,"obs_srv1_se");
  obs_srv1_lci.allocate(1,nyrs_srv1,"obs_srv1_lci");
  obs_srv1_uci.allocate(1,nyrs_srv1,"obs_srv1_uci");
  nyrs_srv2.allocate("nyrs_srv2");
  yrs_srv2.allocate(1,nyrs_srv2,"yrs_srv2");
  obs_srv2_biom.allocate(1,nyrs_srv2,"obs_srv2_biom");
  obs_srv2_se.allocate(1,nyrs_srv2,"obs_srv2_se");
  obs_srv2_lci.allocate(1,nyrs_srv2,"obs_srv2_lci");
  obs_srv2_uci.allocate(1,nyrs_srv2,"obs_srv2_uci");
 cout<<"Number of srv1 rpn is:"<<nyrs_srv2<<endl;
  nyrs_srv3.allocate("nyrs_srv3");
  yrs_srv3.allocate(1,nyrs_srv3,"yrs_srv3");
  obs_srv3_biom.allocate(1,nyrs_srv3,"obs_srv3_biom");
  obs_srv3_se.allocate(1,nyrs_srv3,"obs_srv3_se");
  obs_srv3_lci.allocate(1,nyrs_srv3,"obs_srv3_lci");
  obs_srv3_uci.allocate(1,nyrs_srv3,"obs_srv3_uci");
  nyrs_srv4.allocate("nyrs_srv4");
  yrs_srv4.allocate(1,nyrs_srv4,"yrs_srv4");
  obs_srv4_biom.allocate(1,nyrs_srv4,"obs_srv4_biom");
  obs_srv4_se.allocate(1,nyrs_srv4,"obs_srv4_se");
  obs_srv4_lci.allocate(1,nyrs_srv4,"obs_srv4_lci");
  obs_srv4_uci.allocate(1,nyrs_srv4,"obs_srv4_uci");
  nyrs_srv5.allocate("nyrs_srv5");
  yrs_srv5.allocate(1,nyrs_srv5,"yrs_srv5");
  obs_srv5_biom.allocate(1,nyrs_srv5,"obs_srv5_biom");
  obs_srv5_se.allocate(1,nyrs_srv5,"obs_srv5_se");
  obs_srv5_lci.allocate(1,nyrs_srv5,"obs_srv5_lci");
  obs_srv5_uci.allocate(1,nyrs_srv5,"obs_srv5_uci");
  nyrs_srv6.allocate("nyrs_srv6");
  yrs_srv6.allocate(1,nyrs_srv6,"yrs_srv6");
  obs_srv6_biom.allocate(1,nyrs_srv6,"obs_srv6_biom");
  obs_srv6_se.allocate(1,nyrs_srv6,"obs_srv6_se");
  obs_srv6_lci.allocate(1,nyrs_srv6,"obs_srv6_lci");
  obs_srv6_uci.allocate(1,nyrs_srv6,"obs_srv6_uci");
  nyrs_srv7.allocate("nyrs_srv7");
  yrs_srv7.allocate(1,nyrs_srv7,"yrs_srv7");
  obs_srv7_biom.allocate(1,nyrs_srv7,"obs_srv7_biom");
  obs_srv7_se.allocate(1,nyrs_srv7,"obs_srv7_se");
  obs_srv7_lci.allocate(1,nyrs_srv7,"obs_srv7_lci");
  obs_srv7_uci.allocate(1,nyrs_srv7,"obs_srv7_uci");
 cout<<"Number of trawl years:"<<nyrs_srv7<<endl;
  nyrs_fish1_age.allocate("nyrs_fish1_age");
  yrs_fish1_age.allocate(1,nyrs_fish1_age,"yrs_fish1_age");
  nsamples_fish1_age.allocate(1,nyrs_fish1_age,"nsamples_fish1_age");
  oac_fish1.allocate(1,nyrs_fish1_age,1,nages,"oac_fish1");
  nyrs_srv1_age.allocate("nyrs_srv1_age");
  yrs_srv1_age.allocate(1,nyrs_srv1_age,"yrs_srv1_age");
  nsamples_srv1_age.allocate(1,nyrs_srv1_age,"nsamples_srv1_age");
  oac_srv1.allocate(1,nyrs_srv1_age,1,nages,"oac_srv1");
  nyrs_srv2_age.allocate("nyrs_srv2_age");
  yrs_srv2_age.allocate(1,nyrs_srv2_age,"yrs_srv2_age");
  nsamples_srv2_age.allocate(1,nyrs_srv2_age,"nsamples_srv2_age");
  oac_srv2.allocate(1,nyrs_srv2_age,1,nages,"oac_srv2");
 cout<<"Number of ages srv2:"<<nyrs_srv2_age<<endl;
  nyrs_fish1_size.allocate("nyrs_fish1_size");
  yrs_fish1_size.allocate(1,nyrs_fish1_size,"yrs_fish1_size");
  nsamples_fish1_size.allocate(1,nyrs_fish1_size,"nsamples_fish1_size");
  osc_fish1_m.allocate(1,nyrs_fish1_size,1,nlenbins,"osc_fish1_m");
  osc_fish1_f.allocate(1,nyrs_fish1_size,1,nlenbins,"osc_fish1_f");
  nyrs_fish2_size.allocate("nyrs_fish2_size");
  yrs_fish2_size.allocate(1,nyrs_fish2_size,"yrs_fish2_size");
  nsamples_fish2_size.allocate(1,nyrs_fish2_size,"nsamples_fish2_size");
  osc_fish2.allocate(1,nyrs_fish2_size,1,nlenbins,"osc_fish2");
 cout<<"Number of fish 2 sizes:"<< nyrs_fish2_size <<endl;
  nyrs_fish3_size.allocate("nyrs_fish3_size");
  yrs_fish3_size.allocate(1,nyrs_fish3_size,"yrs_fish3_size");
  nsamples_fish3_size.allocate(1,nyrs_fish3_size,"nsamples_fish3_size");
  osc_fish3_m.allocate(1,nyrs_fish3_size,1,nlenbins,"osc_fish3_m");
  osc_fish3_f.allocate(1,nyrs_fish3_size,1,nlenbins,"osc_fish3_f");
 cout<<"Number of fish 3 sizes:"<< nyrs_fish3_size <<endl;
 cout<<"Years fish 4 sizes:"<< yrs_fish3_size <<endl;
  nyrs_fish4_size.allocate("nyrs_fish4_size");
  yrs_fish4_size.allocate(1,nyrs_fish4_size,"yrs_fish4_size");
  nsamples_fish4_size.allocate(1,nyrs_fish4_size,"nsamples_fish4_size");
  osc_fish4.allocate(1,nyrs_fish4_size,1,nlenbins,"osc_fish4");
 cout<<"Number of fish 4 sizes:"<< nyrs_fish4_size <<endl;
  nyrs_srv1_size.allocate("nyrs_srv1_size");
  yrs_srv1_size.allocate(1,nyrs_srv1_size,"yrs_srv1_size");
  nsamples_srv1_size.allocate(1,nyrs_srv1_size,"nsamples_srv1_size");
  osc_srv1_m.allocate(1,nyrs_srv1_size,1,nlenbins,"osc_srv1_m");
  osc_srv1_f.allocate(1,nyrs_srv1_size,1,nlenbins,"osc_srv1_f");
 cout<<"Number of srv1 size is:"<<nyrs_srv1_size<<endl;
  nyrs_srv2_size.allocate("nyrs_srv2_size");
  yrs_srv2_size.allocate(1,nyrs_srv2_size,"yrs_srv2_size");
  nsamples_srv2_size.allocate(1,nyrs_srv2_size,"nsamples_srv2_size");
  osc_srv2_m.allocate(1,nyrs_srv2_size,1,nlenbins,"osc_srv2_m");
  osc_srv2_f.allocate(1,nyrs_srv2_size,1,nlenbins,"osc_srv2_f");
  nyrs_srv7_size.allocate("nyrs_srv7_size");
  yrs_srv7_size.allocate(1,nyrs_srv7_size,"yrs_srv7_size");
  nsamples_srv7_size.allocate(1,nyrs_srv7_size,"nsamples_srv7_size");
  osc_srv7_m.allocate(1,nyrs_srv7_size,1,nlenbins,"osc_srv7_m");
  osc_srv7_f.allocate(1,nyrs_srv7_size,1,nlenbins,"osc_srv7_f");
  nyrs_srv7_age.allocate("nyrs_srv7_age");
  yrs_srv7_age.allocate(1,nyrs_srv7_age,"yrs_srv7_age");
  nsamples_srv7_age.allocate(1,nyrs_srv7_age,"nsamples_srv7_age");
  oac_srv7.allocate(1,nyrs_srv7_age,1,nages,"oac_srv7");
 cout<<"Number of survey age years is "<<nyrs_srv7_age<<endl;
  sizeage1_m.allocate(1,nages,1,nlenbins,"sizeage1_m");
  sizeage1_f.allocate(1,nages,1,nlenbins,"sizeage1_f");
  sizeage2_m.allocate(1,nages,1,nlenbins,"sizeage2_m");
  sizeage2_f.allocate(1,nages,1,nlenbins,"sizeage2_f");
  sizeage3_m.allocate(1,nages,1,nlenbins,"sizeage3_m");
  sizeage3_f.allocate(1,nages,1,nlenbins,"sizeage3_f");
  sizeage_all.allocate(1,nages,1,nlenbins,"sizeage_all");
  ageage.allocate(1,nages,1,nages,"ageage");
  eof.allocate("eof");
 cout<<"The universal answer is "<<eof<<endl;
  offset.allocate(1,16);
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
}

void model_parameters::initializationfunction(void)
{
  log_mean_rec.set_initial_value(2.0);
  log_sigr.set_initial_value(sigrprior);
  logm.set_initial_value(-2.0);
  log_avg_F_fish1.set_initial_value(-1.0);
  log_avg_F_fish3.set_initial_value(-1.0);
  log_mF50.set_initial_value(-2.0);
  log_mF40.set_initial_value(-2.0);
  log_mF35.set_initial_value(-2.0);
  log_q_srv1.set_initial_value(0.5);
  log_q_srv2.set_initial_value(0.5);
  log_q_srv3.set_initial_value(0.5);
  log_q_srv4.set_initial_value(0.5);
  log_q_srv5.set_initial_value(0.5);
  log_q_srv6.set_initial_value(0.5);
  log_q_srv7.set_initial_value(0.5);
  log_q_srv8.set_initial_value(-1.0);
  log_q_LL_fish_recent.set_initial_value(-1.0);
  log_a50_fish1_f.set_initial_value(1.5);
  log_a50_fish1_m.set_initial_value(1.5);
  log_a50_fish2.set_initial_value(1.5);
  log_a50_fish3_m.set_initial_value(1.0);
  log_a50_fish3_f.set_initial_value(1.0);
  log_a50_fish4_m.set_initial_value(1.5);
  log_a50_fish4_f.set_initial_value(1.5);
  log_a50_fish5_m.set_initial_value(1.5);
  log_a50_fish5_f.set_initial_value(1.5);
  log_a50_srv1_f.set_initial_value(1.5);
  log_a50_srv1_m.set_initial_value(1.5);
  log_a50_srv2_f.set_initial_value(1.5);
  log_a50_srv2_m.set_initial_value(1.5);
  log_a50_srv7_f.set_initial_value(1.0);
  log_a50_srv7_m.set_initial_value(1.0);
  log_a50_srv10_f.set_initial_value(1.5);
  log_a50_srv10_m.set_initial_value(1.5);
  log_delta_fish1_f.set_initial_value(1.5);
  log_delta_fish1_m.set_initial_value(1.5);
  log_delta_fish2.set_initial_value(1.5);
  log_delta_fish3_f.set_initial_value(1.5);
  log_delta_fish3_m.set_initial_value(1.5);
  log_delta_fish4_f.set_initial_value(1.5);
  log_delta_fish4_m.set_initial_value(1.5);
  log_delta_fish5_f.set_initial_value(1.5);
  log_delta_fish5_m.set_initial_value(1.5);
  log_delta_srv1_f.set_initial_value(1.5);
  log_delta_srv1_m.set_initial_value(1.5);
  log_delta_srv2_f.set_initial_value(1.5);
  log_delta_srv2_m.set_initial_value(1.5);
  log_delta_srv10_f.set_initial_value(1.5);
  log_delta_srv10_m.set_initial_value(1.5);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  steepness.allocate(0.2001,0.999,ph_steepness,"steepness");
  log_Rzero.allocate(1,5,ph_Rzero,"log_Rzero");
  log_rec_dev.allocate(styr-nages+2,endyr_rec_est,-10,10,ph_recdev,"log_rec_dev");
  log_mean_rec.allocate(-1,10,ph_mean_rec,"log_mean_rec");
  log_sigr.allocate(-2,1,ph_sigr,"log_sigr");
  log_mF50.allocate(-5,0.5,ph_F50,"log_mF50");
  log_mF40.allocate(-5,0.5,ph_F50,"log_mF40");
  log_mF35.allocate(-5,0.5,ph_F50,"log_mF35");
  logm.allocate(-4,0.5,ph_m,"logm");
  log_M_devs.allocate(styr,endyr,-5,5,ph_Mdevs,"log_M_devs");
  log_M_devs_age.allocate(1,nages,-5,5,ph_Mdevs_age,"log_M_devs_age");
  mdelta.allocate(-0.05,0.05,ph_mdelta,"mdelta");
  natmort.allocate("natmort");
  #ifndef NO_AD_INITIALIZE
  natmort.initialize();
  #endif
  M_year.allocate(styr,endyr,"M_year");
  #ifndef NO_AD_INITIALIZE
    M_year.initialize();
  #endif
  M_age.allocate(1,nages,"M_age");
  #ifndef NO_AD_INITIALIZE
    M_age.initialize();
  #endif
  M.allocate(styr,endyr,1,nages,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  natmortv.allocate(1,nages,"natmortv");
  #ifndef NO_AD_INITIALIZE
    natmortv.initialize();
  #endif
  log_avg_F_fish1.allocate(-10,2,ph_avg_F,"log_avg_F_fish1");
  log_F_devs_fish1.allocate(styr,endyr,-5,5,ph_Fdev,"log_F_devs_fish1");
  log_avg_F_fish3.allocate(-10,2,ph_avg_F,"log_avg_F_fish3");
  log_F_devs_fish3.allocate(styr+3,endyr,-5,5,ph_Fdev,"log_F_devs_fish3");
  Fmort_fish1.allocate(styr,endyr,"Fmort_fish1");
  #ifndef NO_AD_INITIALIZE
    Fmort_fish1.initialize();
  #endif
  Fmort_fish3.allocate(styr,endyr,"Fmort_fish3");
  #ifndef NO_AD_INITIALIZE
    Fmort_fish3.initialize();
  #endif
  F_fish1_f.allocate(styr,endyr,1,nages,"F_fish1_f");
  #ifndef NO_AD_INITIALIZE
    F_fish1_f.initialize();
  #endif
  F_fish1_m.allocate(styr,endyr,1,nages,"F_fish1_m");
  #ifndef NO_AD_INITIALIZE
    F_fish1_m.initialize();
  #endif
  F_fish3_f.allocate(styr,endyr,1,nages,"F_fish3_f");
  #ifndef NO_AD_INITIALIZE
    F_fish3_f.initialize();
  #endif
  F_fish3_m.allocate(styr,endyr,1,nages,"F_fish3_m");
  #ifndef NO_AD_INITIALIZE
    F_fish3_m.initialize();
  #endif
  hist_hal_F.allocate("hist_hal_F");
  #ifndef NO_AD_INITIALIZE
  hist_hal_F.initialize();
  #endif
  log_q_srv1.allocate(-5,5,ph_q_srv1,"log_q_srv1");
  log_q_srv2.allocate(-5,5,ph_q_srv2,"log_q_srv2");
  log_q_srv3.allocate(-5,5,ph_q_srv3,"log_q_srv3");
  log_q_srv4.allocate(-5,5,ph_q_srv4,"log_q_srv4");
  log_q_srv5.allocate(-5,5,ph_q_srv5,"log_q_srv5");
  log_q_srv6.allocate(-5,5,ph_q_srv6,"log_q_srv6");
  log_q_srv7.allocate(-5,5,ph_q_srv7,"log_q_srv7");
  log_q_srv8.allocate(-15,5,ph_q_srv8,"log_q_srv8");
  log_q_srv9.allocate(-5,5,ph_srv2_q2,"log_q_srv9");
  log_q_LL_srvy_recent.allocate(-5,5,ph_q_LL_srv_rec,"log_q_LL_srvy_recent");
  log_q_LL_fish_recent.allocate(-15,5,ph_q_IFQ_rec,"log_q_LL_fish_recent");
  log_fish1_sel_coffs_f.allocate(1,n_fish_sel_ages,phase_selcoff_fsh1,"log_fish1_sel_coffs_f");
  log_a50_fish1_f.allocate(-1,4,ph_fish_sel,"log_a50_fish1_f");
  log_delta_fish1_f.allocate(-5,4,ph_fish_sel_delt,"log_delta_fish1_f");
  a50_fish1_f.allocate("a50_fish1_f");
  #ifndef NO_AD_INITIALIZE
  a50_fish1_f.initialize();
  #endif
  d50_fish1_f.allocate(phase_dlogist_fsh1,"d50_fish1_f");
  delta_fish1_f.allocate("delta_fish1_f");
  #ifndef NO_AD_INITIALIZE
  delta_fish1_f.initialize();
  #endif
  gamma_fish1_f.allocate(phase_dlogist_fsh1,"gamma_fish1_f");
  log_fish1_sel_f.allocate(1,nages,"log_fish1_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_fish1_sel_f.initialize();
  #endif
  fish1_sel_f.allocate(1,nages,"fish1_sel_f");
  #ifndef NO_AD_INITIALIZE
    fish1_sel_f.initialize();
  #endif
  log_avgfish1sel_f.allocate("log_avgfish1sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgfish1sel_f.initialize();
  #endif
  log_fish1_sel_coffs_m.allocate(1,n_fish_sel_ages,phase_selcoff_fsh1,"log_fish1_sel_coffs_m");
  log_a50_fish1_m.allocate(-1,4,ph_fish_sel,"log_a50_fish1_m");
  log_delta_fish1_m.allocate(-5,4,ph_fish_sel_delt_alt,"log_delta_fish1_m");
  a50_fish1_m.allocate("a50_fish1_m");
  #ifndef NO_AD_INITIALIZE
  a50_fish1_m.initialize();
  #endif
  d50_fish1_m.allocate(phase_dlogist_fsh1,"d50_fish1_m");
  delta_fish1_m.allocate("delta_fish1_m");
  #ifndef NO_AD_INITIALIZE
  delta_fish1_m.initialize();
  #endif
  gamma_fish1_m.allocate(phase_dlogist_fsh1,"gamma_fish1_m");
  log_fish1_sel_m.allocate(1,nages,"log_fish1_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_fish1_sel_m.initialize();
  #endif
  fish1_sel_m.allocate(1,nages,"fish1_sel_m");
  #ifndef NO_AD_INITIALIZE
    fish1_sel_m.initialize();
  #endif
  log_avgfish1sel_m.allocate("log_avgfish1sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgfish1sel_m.initialize();
  #endif
  log_fish2_sel_coffs.allocate(1,n_fish_sel_ages,phase_selcoff_fsh2,"log_fish2_sel_coffs");
  log_a50_fish2.allocate(-1,2.5,ph_fish2_sel,"log_a50_fish2");
  log_delta_fish2.allocate(-5,4,ph_fish_sel_delt_alt,"log_delta_fish2");
  a50_fish2.allocate("a50_fish2");
  #ifndef NO_AD_INITIALIZE
  a50_fish2.initialize();
  #endif
  delta_fish2.allocate("delta_fish2");
  #ifndef NO_AD_INITIALIZE
  delta_fish2.initialize();
  #endif
  log_fish2_sel.allocate(1,nages,"log_fish2_sel");
  #ifndef NO_AD_INITIALIZE
    log_fish2_sel.initialize();
  #endif
  fish2_sel.allocate(1,nages,"fish2_sel");
  #ifndef NO_AD_INITIALIZE
    fish2_sel.initialize();
  #endif
  log_avgfish2sel.allocate("log_avgfish2sel");
  #ifndef NO_AD_INITIALIZE
  log_avgfish2sel.initialize();
  #endif
  log_fish3_sel_coffs_f.allocate(1,n_fish_sel_ages,phase_selcoff_fsh3,"log_fish3_sel_coffs_f");
  log_a50_fish3_f.allocate(-1,4,ph_fish_sel,"log_a50_fish3_f");
  log_delta_fish3_f.allocate(-5,4,ph_fish_sel_delt,"log_delta_fish3_f");
  a50_fish3_f.allocate("a50_fish3_f");
  #ifndef NO_AD_INITIALIZE
  a50_fish3_f.initialize();
  #endif
  d50_fish3_f.allocate(-1,"d50_fish3_f");
  delta_fish3_f.allocate("delta_fish3_f");
  #ifndef NO_AD_INITIALIZE
  delta_fish3_f.initialize();
  #endif
  log_gamma_fish3_f.allocate(-0.00001,-0.00000000001,-4,"log_gamma_fish3_f");
  log_fish3_sel_f.allocate(1,nages,"log_fish3_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_fish3_sel_f.initialize();
  #endif
  fish3_sel_f.allocate(1,nages,"fish3_sel_f");
  #ifndef NO_AD_INITIALIZE
    fish3_sel_f.initialize();
  #endif
  log_avgfish3sel_f.allocate("log_avgfish3sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgfish3sel_f.initialize();
  #endif
  gamma_fish3_f.allocate("gamma_fish3_f");
  #ifndef NO_AD_INITIALIZE
  gamma_fish3_f.initialize();
  #endif
  log_fish3_sel_coffs_m.allocate(1,n_fish_sel_ages,phase_selcoff_fsh3,"log_fish3_sel_coffs_m");
  log_a50_fish3_m.allocate(-1,4,ph_fish_sel,"log_a50_fish3_m");
  log_delta_fish3_m.allocate(-5,4,ph_fish_sel_delt_alt,"log_delta_fish3_m");
  a50_fish3_m.allocate("a50_fish3_m");
  #ifndef NO_AD_INITIALIZE
  a50_fish3_m.initialize();
  #endif
  d50_fish3_m.allocate(-1,"d50_fish3_m");
  delta_fish3_m.allocate("delta_fish3_m");
  #ifndef NO_AD_INITIALIZE
  delta_fish3_m.initialize();
  #endif
  log_gamma_fish3_m.allocate(-0.18,-0.173,-4,"log_gamma_fish3_m");
  log_fish3_sel_m.allocate(1,nages,"log_fish3_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_fish3_sel_m.initialize();
  #endif
  fish3_sel_m.allocate(1,nages,"fish3_sel_m");
  #ifndef NO_AD_INITIALIZE
    fish3_sel_m.initialize();
  #endif
  log_avgfish3sel_m.allocate("log_avgfish3sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgfish3sel_m.initialize();
  #endif
  gamma_fish3_m.allocate("gamma_fish3_m");
  #ifndef NO_AD_INITIALIZE
  gamma_fish3_m.initialize();
  #endif
  log_fish4_sel_coffs_f.allocate(1,n_fish_sel_ages,phase_selcoff_fsh4,"log_fish4_sel_coffs_f");
  log_a50_fish4_f.allocate(-1,4,ph_ifq,"log_a50_fish4_f");
  log_delta_fish4_f.allocate(-5,4,ph_fish_sel_delt,"log_delta_fish4_f");
  a50_fish4_f.allocate("a50_fish4_f");
  #ifndef NO_AD_INITIALIZE
  a50_fish4_f.initialize();
  #endif
  d50_fish4_f.allocate(phase_dlogist_fsh4,"d50_fish4_f");
  delta_fish4_f.allocate("delta_fish4_f");
  #ifndef NO_AD_INITIALIZE
  delta_fish4_f.initialize();
  #endif
  gamma_fish4_f.allocate(phase_dlogist_fsh4,"gamma_fish4_f");
  log_fish4_sel_f.allocate(1,nages,"log_fish4_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_fish4_sel_f.initialize();
  #endif
  fish4_sel_f.allocate(1,nages,"fish4_sel_f");
  #ifndef NO_AD_INITIALIZE
    fish4_sel_f.initialize();
  #endif
  log_avgfish4sel_f.allocate("log_avgfish4sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgfish4sel_f.initialize();
  #endif
  log_fish4_sel_coffs_m.allocate(1,n_fish_sel_ages,phase_selcoff_fsh4,"log_fish4_sel_coffs_m");
  log_a50_fish4_m.allocate(-1,4,ph_ifq,"log_a50_fish4_m");
  log_delta_fish4_m.allocate(-5,4,ph_fish_sel_delt,"log_delta_fish4_m");
  a50_fish4_m.allocate("a50_fish4_m");
  #ifndef NO_AD_INITIALIZE
  a50_fish4_m.initialize();
  #endif
  d50_fish4_m.allocate(phase_dlogist_fsh4,"d50_fish4_m");
  delta_fish4_m.allocate("delta_fish4_m");
  #ifndef NO_AD_INITIALIZE
  delta_fish4_m.initialize();
  #endif
  gamma_fish4_m.allocate(-4,"gamma_fish4_m");
  log_fish4_sel_m.allocate(1,nages,"log_fish4_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_fish4_sel_m.initialize();
  #endif
  fish4_sel_m.allocate(1,nages,"fish4_sel_m");
  #ifndef NO_AD_INITIALIZE
    fish4_sel_m.initialize();
  #endif
  log_avgfish4sel_m.allocate("log_avgfish4sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgfish4sel_m.initialize();
  #endif
  log_fish5_sel_coffs_f.allocate(1,n_fish_sel_ages,phase_selcoff_fsh5,"log_fish5_sel_coffs_f");
  log_a50_fish5_f.allocate(-1,4,ph_ifq_block2,"log_a50_fish5_f");
  log_delta_fish5_f.allocate(-5,4,ph_ifq_block2,"log_delta_fish5_f");
  a50_fish5_f.allocate("a50_fish5_f");
  #ifndef NO_AD_INITIALIZE
  a50_fish5_f.initialize();
  #endif
  d50_fish5_f.allocate(phase_dlogist_fsh5,"d50_fish5_f");
  delta_fish5_f.allocate("delta_fish5_f");
  #ifndef NO_AD_INITIALIZE
  delta_fish5_f.initialize();
  #endif
  gamma_fish5_f.allocate(phase_dlogist_fsh5,"gamma_fish5_f");
  log_fish5_sel_f.allocate(1,nages,"log_fish5_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_fish5_sel_f.initialize();
  #endif
  fish5_sel_f.allocate(1,nages,"fish5_sel_f");
  #ifndef NO_AD_INITIALIZE
    fish5_sel_f.initialize();
  #endif
  log_avgfish5sel_f.allocate("log_avgfish5sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgfish5sel_f.initialize();
  #endif
  log_fish5_sel_coffs_m.allocate(1,n_fish_sel_ages,phase_selcoff_fsh5,"log_fish5_sel_coffs_m");
  log_a50_fish5_m.allocate(-1,4,ph_ifq_block2,"log_a50_fish5_m");
  log_delta_fish5_m.allocate(-5,4,ph_ifq_block2,"log_delta_fish5_m");
  a50_fish5_m.allocate("a50_fish5_m");
  #ifndef NO_AD_INITIALIZE
  a50_fish5_m.initialize();
  #endif
  d50_fish5_m.allocate(phase_dlogist_fsh5,"d50_fish5_m");
  delta_fish5_m.allocate("delta_fish5_m");
  #ifndef NO_AD_INITIALIZE
  delta_fish5_m.initialize();
  #endif
  gamma_fish5_m.allocate(-4,"gamma_fish5_m");
  log_fish5_sel_m.allocate(1,nages,"log_fish5_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_fish5_sel_m.initialize();
  #endif
  fish5_sel_m.allocate(1,nages,"fish5_sel_m");
  #ifndef NO_AD_INITIALIZE
    fish5_sel_m.initialize();
  #endif
  log_avgfish5sel_m.allocate("log_avgfish5sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgfish5sel_m.initialize();
  #endif
  sel_rep_proj_f.allocate(1,nages,"sel_rep_proj_f");
  #ifndef NO_AD_INITIALIZE
    sel_rep_proj_f.initialize();
  #endif
  sel_rep_proj_m.allocate(1,nages,"sel_rep_proj_m");
  #ifndef NO_AD_INITIALIZE
    sel_rep_proj_m.initialize();
  #endif
  log_srv1_sel_coffs_f.allocate(1,n_srv1_sel_ages,phase_selcoff_srv1,"log_srv1_sel_coffs_f");
  log_a50_srv1_f.allocate(-1,4,ph_srv1_sel,"log_a50_srv1_f");
  log_delta_srv1_f.allocate(-5,4,ph_srv_sel_delt,"log_delta_srv1_f");
  a50_srv1_f.allocate("a50_srv1_f");
  #ifndef NO_AD_INITIALIZE
  a50_srv1_f.initialize();
  #endif
  d50_srv1_f.allocate(phase_dlogist_srv1,"d50_srv1_f");
  delta_srv1_f.allocate("delta_srv1_f");
  #ifndef NO_AD_INITIALIZE
  delta_srv1_f.initialize();
  #endif
  gamma_srv1_f.allocate(phase_dlogist_srv1,"gamma_srv1_f");
  log_srv1_sel_f.allocate(1,nages,"log_srv1_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_srv1_sel_f.initialize();
  #endif
  srv1_sel_f.allocate(1,nages,"srv1_sel_f");
  #ifndef NO_AD_INITIALIZE
    srv1_sel_f.initialize();
  #endif
  log_avgsrv1sel_f.allocate("log_avgsrv1sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv1sel_f.initialize();
  #endif
  log_srv1_sel_coffs_m.allocate(1,n_srv1_sel_ages,phase_selcoff_srv1,"log_srv1_sel_coffs_m");
  log_a50_srv1_m.allocate(-1,4,ph_srv1_sel,"log_a50_srv1_m");
  log_delta_srv1_m.allocate(-5,4,ph_srv_sel_delt,"log_delta_srv1_m");
  a50_srv1_m.allocate("a50_srv1_m");
  #ifndef NO_AD_INITIALIZE
  a50_srv1_m.initialize();
  #endif
  d50_srv1_m.allocate(phase_dlogist_srv1,"d50_srv1_m");
  delta_srv1_m.allocate("delta_srv1_m");
  #ifndef NO_AD_INITIALIZE
  delta_srv1_m.initialize();
  #endif
  gamma_srv1_m.allocate(phase_dlogist_srv1,"gamma_srv1_m");
  log_srv1_sel_m.allocate(1,nages,"log_srv1_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_srv1_sel_m.initialize();
  #endif
  srv1_sel_m.allocate(1,nages,"srv1_sel_m");
  #ifndef NO_AD_INITIALIZE
    srv1_sel_m.initialize();
  #endif
  log_avgsrv1sel_m.allocate("log_avgsrv1sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv1sel_m.initialize();
  #endif
  log_srv2_sel_coffs_f.allocate(1,n_srv1_sel_ages,phase_selcoff_srv2,"log_srv2_sel_coffs_f");
  log_a50_srv2_f.allocate(-1,4,ph_srv2_sel,"log_a50_srv2_f");
  log_delta_srv2_f.allocate(-5,4,ph_srv_sel_delt_alt,"log_delta_srv2_f");
  a50_srv2_f.allocate("a50_srv2_f");
  #ifndef NO_AD_INITIALIZE
  a50_srv2_f.initialize();
  #endif
  d50_srv2_f.allocate(phase_dlogist_srv2,"d50_srv2_f");
  gamma_srv2_f.allocate(phase_dlogist_srv2,"gamma_srv2_f");
  delta_srv2_f.allocate("delta_srv2_f");
  #ifndef NO_AD_INITIALIZE
  delta_srv2_f.initialize();
  #endif
  log_srv2_sel_f.allocate(1,nages,"log_srv2_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_srv2_sel_f.initialize();
  #endif
  srv2_sel_f.allocate(1,nages,"srv2_sel_f");
  #ifndef NO_AD_INITIALIZE
    srv2_sel_f.initialize();
  #endif
  log_avgsrv2sel_f.allocate("log_avgsrv2sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv2sel_f.initialize();
  #endif
  log_srv2_sel_coffs_m.allocate(1,n_srv1_sel_ages,phase_selcoff_srv2,"log_srv2_sel_coffs_m");
  log_a50_srv2_m.allocate(-1,4,ph_srv2_sel,"log_a50_srv2_m");
  log_delta_srv2_m.allocate(-5,4,ph_srv_sel_delt_alt,"log_delta_srv2_m");
  a50_srv2_m.allocate("a50_srv2_m");
  #ifndef NO_AD_INITIALIZE
  a50_srv2_m.initialize();
  #endif
  d50_srv2_m.allocate(phase_dlogist_srv2,"d50_srv2_m");
  gamma_srv2_m.allocate(phase_dlogist_srv2,"gamma_srv2_m");
  delta_srv2_m.allocate("delta_srv2_m");
  #ifndef NO_AD_INITIALIZE
  delta_srv2_m.initialize();
  #endif
  log_srv2_sel_m.allocate(1,nages,"log_srv2_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_srv2_sel_m.initialize();
  #endif
  srv2_sel_m.allocate(1,nages,"srv2_sel_m");
  #ifndef NO_AD_INITIALIZE
    srv2_sel_m.initialize();
  #endif
  log_avgsrv2sel_m.allocate("log_avgsrv2sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv2sel_m.initialize();
  #endif
  log_srv7_sel_coffs_f.allocate(1,n_srv1_sel_ages,phase_selcoff_srv7,"log_srv7_sel_coffs_f");
  log_a50_srv7_f.allocate(-3,5,ph_srv1_sel,"log_a50_srv7_f");
  log_delta_srv7_f.allocate(-4,4,-4,"log_delta_srv7_f");
  a50_srv7_f.allocate("a50_srv7_f");
  #ifndef NO_AD_INITIALIZE
  a50_srv7_f.initialize();
  #endif
  log_d50_srv7_f.allocate(phase_dlogist_srv7,"log_d50_srv7_f");
  delta_srv7_f.allocate("delta_srv7_f");
  #ifndef NO_AD_INITIALIZE
  delta_srv7_f.initialize();
  #endif
  log_gamma_srv7_f.allocate(-10,-0.0000000000000001,-4,"log_gamma_srv7_f");
  log_srv7_sel_f.allocate(1,nages,"log_srv7_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_srv7_sel_f.initialize();
  #endif
  srv7_sel_f.allocate(1,nages,"srv7_sel_f");
  #ifndef NO_AD_INITIALIZE
    srv7_sel_f.initialize();
  #endif
  log_avgsrv7sel_f.allocate("log_avgsrv7sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv7sel_f.initialize();
  #endif
  gamma_srv7_f.allocate("gamma_srv7_f");
  #ifndef NO_AD_INITIALIZE
  gamma_srv7_f.initialize();
  #endif
  d50_srv7_f.allocate("d50_srv7_f");
  #ifndef NO_AD_INITIALIZE
  d50_srv7_f.initialize();
  #endif
  log_srv7_sel_coffs_m.allocate(1,n_srv1_sel_ages,phase_selcoff_srv7,"log_srv7_sel_coffs_m");
  log_a50_srv7_m.allocate(-3,5,ph_srv1_sel,"log_a50_srv7_m");
  log_delta_srv7_m.allocate(-4,4,-4,"log_delta_srv7_m");
  a50_srv7_m.allocate("a50_srv7_m");
  #ifndef NO_AD_INITIALIZE
  a50_srv7_m.initialize();
  #endif
  log_d50_srv7_m.allocate(phase_dlogist_srv7,"log_d50_srv7_m");
  delta_srv7_m.allocate("delta_srv7_m");
  #ifndef NO_AD_INITIALIZE
  delta_srv7_m.initialize();
  #endif
  log_gamma_srv7_m.allocate(-10,-0.0000000000000001,-4,"log_gamma_srv7_m");
  log_srv7_sel_m.allocate(1,nages,"log_srv7_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_srv7_sel_m.initialize();
  #endif
  srv7_sel_m.allocate(1,nages,"srv7_sel_m");
  #ifndef NO_AD_INITIALIZE
    srv7_sel_m.initialize();
  #endif
  log_avgsrv7sel_m.allocate("log_avgsrv7sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv7sel_m.initialize();
  #endif
  gamma_srv7_m.allocate("gamma_srv7_m");
  #ifndef NO_AD_INITIALIZE
  gamma_srv7_m.initialize();
  #endif
  d50_srv7_m.allocate("d50_srv7_m");
  #ifndef NO_AD_INITIALIZE
  d50_srv7_m.initialize();
  #endif
  log_srv10_sel_coffs_f.allocate(1,n_srv1_sel_ages,phase_selcoff_srv10,"log_srv10_sel_coffs_f");
  log_a50_srv10_f.allocate(-1,4,ph_LL_block2,"log_a50_srv10_f");
  log_delta_srv10_f.allocate(-5,4,ph_srv_sel_delt_alt,"log_delta_srv10_f");
  a50_srv10_f.allocate("a50_srv10_f");
  #ifndef NO_AD_INITIALIZE
  a50_srv10_f.initialize();
  #endif
  d50_srv10_f.allocate(phase_dlogist_srv10,"d50_srv10_f");
  delta_srv10_f.allocate("delta_srv10_f");
  #ifndef NO_AD_INITIALIZE
  delta_srv10_f.initialize();
  #endif
  gamma_srv10_f.allocate(phase_dlogist_srv10,"gamma_srv10_f");
  log_srv10_sel_f.allocate(1,nages,"log_srv10_sel_f");
  #ifndef NO_AD_INITIALIZE
    log_srv10_sel_f.initialize();
  #endif
  srv10_sel_f.allocate(1,nages,"srv10_sel_f");
  #ifndef NO_AD_INITIALIZE
    srv10_sel_f.initialize();
  #endif
  log_avgsrv10sel_f.allocate("log_avgsrv10sel_f");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv10sel_f.initialize();
  #endif
  log_srv10_sel_coffs_m.allocate(1,n_srv1_sel_ages,phase_selcoff_srv10,"log_srv10_sel_coffs_m");
  log_a50_srv10_m.allocate(-1,4,ph_LL_block2,"log_a50_srv10_m");
  log_delta_srv10_m.allocate(-5,4,ph_srv_sel_delt_alt,"log_delta_srv10_m");
  a50_srv10_m.allocate("a50_srv10_m");
  #ifndef NO_AD_INITIALIZE
  a50_srv10_m.initialize();
  #endif
  d50_srv10_m.allocate(phase_dlogist_srv10,"d50_srv10_m");
  delta_srv10_m.allocate("delta_srv10_m");
  #ifndef NO_AD_INITIALIZE
  delta_srv10_m.initialize();
  #endif
  gamma_srv10_m.allocate(phase_dlogist_srv10,"gamma_srv10_m");
  log_srv10_sel_m.allocate(1,nages,"log_srv10_sel_m");
  #ifndef NO_AD_INITIALIZE
    log_srv10_sel_m.initialize();
  #endif
  srv10_sel_m.allocate(1,nages,"srv10_sel_m");
  #ifndef NO_AD_INITIALIZE
    srv10_sel_m.initialize();
  #endif
  log_avgsrv10sel_m.allocate("log_avgsrv10sel_m");
  #ifndef NO_AD_INITIALIZE
  log_avgsrv10sel_m.initialize();
  #endif
  F50.allocate("F50");
  #ifndef NO_AD_INITIALIZE
  F50.initialize();
  #endif
  F40.allocate("F40");
  F35.allocate("F35");
  #ifndef NO_AD_INITIALIZE
  F35.initialize();
  #endif
  mF50.allocate("mF50");
  #ifndef NO_AD_INITIALIZE
  mF50.initialize();
  #endif
  mF40.allocate("mF40");
  #ifndef NO_AD_INITIALIZE
  mF40.initialize();
  #endif
  mF35.allocate("mF35");
  #ifndef NO_AD_INITIALIZE
  mF35.initialize();
  #endif
  SB0.allocate("SB0");
  #ifndef NO_AD_INITIALIZE
  SB0.initialize();
  #endif
  SBF50.allocate("SBF50");
  #ifndef NO_AD_INITIALIZE
  SBF50.initialize();
  #endif
  SBF40.allocate("SBF40");
  #ifndef NO_AD_INITIALIZE
  SBF40.initialize();
  #endif
  SBF35.allocate("SBF35");
  #ifndef NO_AD_INITIALIZE
  SBF35.initialize();
  #endif
  sprpen.allocate("sprpen");
  #ifndef NO_AD_INITIALIZE
  sprpen.initialize();
  #endif
  Nspr.allocate(1,4,1,nages,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  Z_f.allocate(styr,endyr,1,nages,"Z_f");
  #ifndef NO_AD_INITIALIZE
    Z_f.initialize();
  #endif
  Z_m.allocate(styr,endyr,1,nages,"Z_m");
  #ifndef NO_AD_INITIALIZE
    Z_m.initialize();
  #endif
  S_f.allocate(styr,endyr,1,nages,"S_f");
  #ifndef NO_AD_INITIALIZE
    S_f.initialize();
  #endif
  S_m.allocate(styr,endyr,1,nages,"S_m");
  #ifndef NO_AD_INITIALIZE
    S_m.initialize();
  #endif
  S_f_mid.allocate(styr,endyr,1,nages,"S_f_mid");
  #ifndef NO_AD_INITIALIZE
    S_f_mid.initialize();
  #endif
  S_m_mid.allocate(styr,endyr,1,nages,"S_m_mid");
  #ifndef NO_AD_INITIALIZE
    S_m_mid.initialize();
  #endif
  q_srv1.allocate("q_srv1");
  q_srv2.allocate("q_srv2");
  q_srv3.allocate("q_srv3");
  #ifndef NO_AD_INITIALIZE
  q_srv3.initialize();
  #endif
  q_srv4.allocate("q_srv4");
  #ifndef NO_AD_INITIALIZE
  q_srv4.initialize();
  #endif
  q_srv5.allocate("q_srv5");
  #ifndef NO_AD_INITIALIZE
  q_srv5.initialize();
  #endif
  q_srv6.allocate("q_srv6");
  q_srv7.allocate("q_srv7");
  q_srv8.allocate("q_srv8");
  q_srv9.allocate("q_srv9");
  #ifndef NO_AD_INITIALIZE
  q_srv9.initialize();
  #endif
  q_LL_fish_recent.allocate("q_LL_fish_recent");
  q_LL_srvy_recent.allocate("q_LL_srvy_recent");
  #ifndef NO_AD_INITIALIZE
  q_LL_srvy_recent.initialize();
  #endif
  M_est.allocate("M_est");
  sam_rec.allocate(styr_rec,endyr,"sam_rec");
  #ifndef NO_AD_INITIALIZE
    sam_rec.initialize();
  #endif
  srm_rec.allocate(styr_rec,endyr,"srm_rec");
  #ifndef NO_AD_INITIALIZE
    srm_rec.initialize();
  #endif
  sigrsq.allocate("sigrsq");
  #ifndef NO_AD_INITIALIZE
  sigrsq.initialize();
  #endif
  alpha.allocate("alpha");
  #ifndef NO_AD_INITIALIZE
  alpha.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  Bzero.allocate("Bzero");
  #ifndef NO_AD_INITIALIZE
  Bzero.initialize();
  #endif
  Rzero.allocate("Rzero");
  #ifndef NO_AD_INITIALIZE
  Rzero.initialize();
  #endif
  phizero.allocate("phizero");
  #ifndef NO_AD_INITIALIZE
  phizero.initialize();
  #endif
  log_Rztemp.allocate("log_Rztemp");
  #ifndef NO_AD_INITIALIZE
  log_Rztemp.initialize();
  #endif
  pred_rec.allocate(styr,endyr,"pred_rec");
  avg_rec.allocate("avg_rec");
  #ifndef NO_AD_INITIALIZE
  avg_rec.initialize();
  #endif
  b_yr_end.allocate("b_yr_end");
  #ifndef NO_AD_INITIALIZE
  b_yr_end.initialize();
  #endif
  b.allocate(styr-nages+2,endyr_rec_est,"b");
  #ifndef NO_AD_INITIALIZE
    b.initialize();
  #endif
  rec_like_bias_adj.allocate(styr-nages+2,endyr_rec_est,"rec_like_bias_adj");
  #ifndef NO_AD_INITIALIZE
    rec_like_bias_adj.initialize();
  #endif
  sigr.allocate("sigr");
  #ifndef NO_AD_INITIALIZE
  sigr.initialize();
  #endif
  natage_m.allocate(styr,endyr,1,nages,"natage_m");
  #ifndef NO_AD_INITIALIZE
    natage_m.initialize();
  #endif
  natage_f.allocate(styr,endyr,1,nages,"natage_f");
  #ifndef NO_AD_INITIALIZE
    natage_f.initialize();
  #endif
  num_len_m.allocate(styr,endyr,1,nlenbins,"num_len_m");
  #ifndef NO_AD_INITIALIZE
    num_len_m.initialize();
  #endif
  num_len_f.allocate(styr,endyr,1,nlenbins,"num_len_f");
  #ifndef NO_AD_INITIALIZE
    num_len_f.initialize();
  #endif
  tot_biom.allocate(styr,endyr,"tot_biom");
  #ifndef NO_AD_INITIALIZE
    tot_biom.initialize();
  #endif
  spbiom_trend.allocate("spbiom_trend");
  #ifndef NO_AD_INITIALIZE
  spbiom_trend.initialize();
  #endif
  Depletion.allocate("Depletion");
  #ifndef NO_AD_INITIALIZE
  Depletion.initialize();
  #endif
  spawn_biom.allocate(styr,endyr,"spawn_biom");
  #ifndef NO_AD_INITIALIZE
    spawn_biom.initialize();
  #endif
  Sp_Biom.allocate(styr_sp,endyr,"Sp_Biom");
  #ifndef NO_AD_INITIALIZE
    Sp_Biom.initialize();
  #endif
  ssbsd.allocate(styr,endyr,"ssbsd");
  catage_fish1_f.allocate(styr,endyr,1,nages,"catage_fish1_f");
  #ifndef NO_AD_INITIALIZE
    catage_fish1_f.initialize();
  #endif
  catage_fish1_m.allocate(styr,endyr,1,nages,"catage_fish1_m");
  #ifndef NO_AD_INITIALIZE
    catage_fish1_m.initialize();
  #endif
  catage_fish3_f.allocate(styr,endyr,1,nages,"catage_fish3_f");
  #ifndef NO_AD_INITIALIZE
    catage_fish3_f.initialize();
  #endif
  catage_fish3_m.allocate(styr,endyr,1,nages,"catage_fish3_m");
  #ifndef NO_AD_INITIALIZE
    catage_fish3_m.initialize();
  #endif
  pred_catch_fish1.allocate(styr,endyr,"pred_catch_fish1");
  #ifndef NO_AD_INITIALIZE
    pred_catch_fish1.initialize();
  #endif
  pred_catch_fish3.allocate(styr,endyr,"pred_catch_fish3");
  #ifndef NO_AD_INITIALIZE
    pred_catch_fish3.initialize();
  #endif
  pred_srv1.allocate(styr,endyr,"pred_srv1");
  #ifndef NO_AD_INITIALIZE
    pred_srv1.initialize();
  #endif
  pred_srv2.allocate(styr,endyr,"pred_srv2");
  #ifndef NO_AD_INITIALIZE
    pred_srv2.initialize();
  #endif
  pred_srv3.allocate(styr,endyr,"pred_srv3");
  #ifndef NO_AD_INITIALIZE
    pred_srv3.initialize();
  #endif
  pred_srv4.allocate(styr,endyr,"pred_srv4");
  #ifndef NO_AD_INITIALIZE
    pred_srv4.initialize();
  #endif
  pred_srv5.allocate(styr,endyr,"pred_srv5");
  #ifndef NO_AD_INITIALIZE
    pred_srv5.initialize();
  #endif
  pred_srv6.allocate(styr,endyr,"pred_srv6");
  #ifndef NO_AD_INITIALIZE
    pred_srv6.initialize();
  #endif
  pred_srv7.allocate(styr,endyr,"pred_srv7");
  #ifndef NO_AD_INITIALIZE
    pred_srv7.initialize();
  #endif
  pred_srv8.allocate(styr,endyr,"pred_srv8");
  #ifndef NO_AD_INITIALIZE
    pred_srv8.initialize();
  #endif
  eac_fish1.allocate(1,nyrs_fish1_age,1,nages,"eac_fish1");
  #ifndef NO_AD_INITIALIZE
    eac_fish1.initialize();
  #endif
  eac_srv1.allocate(1,nyrs_srv1_age,1,nages,"eac_srv1");
  #ifndef NO_AD_INITIALIZE
    eac_srv1.initialize();
  #endif
  eac_srv2.allocate(1,nyrs_srv2_age,1,nages,"eac_srv2");
  #ifndef NO_AD_INITIALIZE
    eac_srv2.initialize();
  #endif
  esc_fish1_m.allocate(1,nyrs_fish1_size,1,nlenbins,"esc_fish1_m");
  #ifndef NO_AD_INITIALIZE
    esc_fish1_m.initialize();
  #endif
  esc_fish1_f.allocate(1,nyrs_fish1_size,1,nlenbins,"esc_fish1_f");
  #ifndef NO_AD_INITIALIZE
    esc_fish1_f.initialize();
  #endif
  esc_fish2.allocate(1,nyrs_fish2_size,1,nlenbins,"esc_fish2");
  #ifndef NO_AD_INITIALIZE
    esc_fish2.initialize();
  #endif
  esc_fish3_m.allocate(1,nyrs_fish3_size,1,nlenbins,"esc_fish3_m");
  #ifndef NO_AD_INITIALIZE
    esc_fish3_m.initialize();
  #endif
  esc_fish3_f.allocate(1,nyrs_fish3_size,1,nlenbins,"esc_fish3_f");
  #ifndef NO_AD_INITIALIZE
    esc_fish3_f.initialize();
  #endif
  esc_fish4.allocate(1,nyrs_fish4_size,1,nlenbins,"esc_fish4");
  #ifndef NO_AD_INITIALIZE
    esc_fish4.initialize();
  #endif
  esc_srv1_m.allocate(1,nyrs_srv1_size,1,nlenbins,"esc_srv1_m");
  #ifndef NO_AD_INITIALIZE
    esc_srv1_m.initialize();
  #endif
  esc_srv1_f.allocate(1,nyrs_srv1_size,1,nlenbins,"esc_srv1_f");
  #ifndef NO_AD_INITIALIZE
    esc_srv1_f.initialize();
  #endif
  esc_srv2_m.allocate(1,nyrs_srv2_size,1,nlenbins,"esc_srv2_m");
  #ifndef NO_AD_INITIALIZE
    esc_srv2_m.initialize();
  #endif
  esc_srv2_f.allocate(1,nyrs_srv2_size,1,nlenbins,"esc_srv2_f");
  #ifndef NO_AD_INITIALIZE
    esc_srv2_f.initialize();
  #endif
  esc_srv7_m.allocate(1,nyrs_srv7_size,1,nlenbins,"esc_srv7_m");
  #ifndef NO_AD_INITIALIZE
    esc_srv7_m.initialize();
  #endif
  esc_srv7_f.allocate(1,nyrs_srv7_size,1,nlenbins,"esc_srv7_f");
  #ifndef NO_AD_INITIALIZE
    esc_srv7_f.initialize();
  #endif
  eac_srv7.allocate(1,nyrs_srv7_age,1,nages,"eac_srv7");
  #ifndef NO_AD_INITIALIZE
    eac_srv7.initialize();
  #endif
  surv_like.allocate(1,8,"surv_like");
  #ifndef NO_AD_INITIALIZE
    surv_like.initialize();
  #endif
  age_like.allocate(1,16,"age_like");
  #ifndef NO_AD_INITIALIZE
    age_like.initialize();
  #endif
  sel_like.allocate(1,12,"sel_like");
  #ifndef NO_AD_INITIALIZE
    sel_like.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  ssqcatch.allocate("ssqcatch");
  #ifndef NO_AD_INITIALIZE
  ssqcatch.initialize();
  #endif
  F_mort_regularity.allocate("F_mort_regularity");
  #ifndef NO_AD_INITIALIZE
  F_mort_regularity.initialize();
  #endif
  M_mort_regularity.allocate("M_mort_regularity");
  #ifndef NO_AD_INITIALIZE
  M_mort_regularity.initialize();
  #endif
  avg_sel_penalty.allocate("avg_sel_penalty");
  #ifndef NO_AD_INITIALIZE
  avg_sel_penalty.initialize();
  #endif
  priors.allocate(1,13,"priors");
  #ifndef NO_AD_INITIALIZE
    priors.initialize();
  #endif
  Like.allocate("Like");
  #ifndef NO_AD_INITIALIZE
  Like.initialize();
  #endif
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  xdum2.allocate(styr,endyr,"xdum2");
  #ifndef NO_AD_INITIALIZE
    xdum2.initialize();
  #endif
  pred_catch.allocate(styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  fratio.allocate("fratio");
  #ifndef NO_AD_INITIALIZE
  fratio.initialize();
  #endif
  ABC3.allocate("ABC3");
  #ifndef NO_AD_INITIALIZE
  ABC3.initialize();
  #endif
  N_proj_f.allocate(endyr+1,endyr+15,1,nages,"N_proj_f");
  #ifndef NO_AD_INITIALIZE
    N_proj_f.initialize();
  #endif
  N_proj_m.allocate(endyr+1,endyr+15,1,nages,"N_proj_m");
  #ifndef NO_AD_INITIALIZE
    N_proj_m.initialize();
  #endif
  FABC_proj.allocate("FABC_proj");
  #ifndef NO_AD_INITIALIZE
  FABC_proj.initialize();
  #endif
  FABC_tot_proj_f.allocate(1,nages,"FABC_tot_proj_f");
  #ifndef NO_AD_INITIALIZE
    FABC_tot_proj_f.initialize();
  #endif
  FABC_tot_proj_m.allocate(1,nages,"FABC_tot_proj_m");
  #ifndef NO_AD_INITIALIZE
    FABC_tot_proj_m.initialize();
  #endif
  FOFL_proj.allocate("FOFL_proj");
  #ifndef NO_AD_INITIALIZE
  FOFL_proj.initialize();
  #endif
  FOFL_tot_proj_f.allocate(1,nages,"FOFL_tot_proj_f");
  #ifndef NO_AD_INITIALIZE
    FOFL_tot_proj_f.initialize();
  #endif
  FOFL_tot_proj_m.allocate(1,nages,"FOFL_tot_proj_m");
  #ifndef NO_AD_INITIALIZE
    FOFL_tot_proj_m.initialize();
  #endif
  ABC.allocate("ABC");
  B40.allocate("B40");
  OFL.allocate("OFL");
  #ifndef NO_AD_INITIALIZE
  OFL.initialize();
  #endif
  Z_proj_f.allocate(1,nages,"Z_proj_f");
  #ifndef NO_AD_INITIALIZE
    Z_proj_f.initialize();
  #endif
  Z_proj_m.allocate(1,nages,"Z_proj_m");
  #ifndef NO_AD_INITIALIZE
    Z_proj_m.initialize();
  #endif
  ZOFL_proj_f.allocate(1,nages,"ZOFL_proj_f");
  #ifndef NO_AD_INITIALIZE
    ZOFL_proj_f.initialize();
  #endif
  ZOFL_proj_m.allocate(1,nages,"ZOFL_proj_m");
  #ifndef NO_AD_INITIALIZE
    ZOFL_proj_m.initialize();
  #endif
  S_proj_f.allocate(1,nages,"S_proj_f");
  #ifndef NO_AD_INITIALIZE
    S_proj_f.initialize();
  #endif
  S_proj_m.allocate(1,nages,"S_proj_m");
  #ifndef NO_AD_INITIALIZE
    S_proj_m.initialize();
  #endif
  catage_proj_f.allocate(endyr+1,endyr+15,1,nages,"catage_proj_f");
  #ifndef NO_AD_INITIALIZE
    catage_proj_f.initialize();
  #endif
  catage_proj_m.allocate(endyr+1,endyr+15,1,nages,"catage_proj_m");
  #ifndef NO_AD_INITIALIZE
    catage_proj_m.initialize();
  #endif
  catage_proj_OFL_f.allocate(endyr+1,endyr+15,1,nages,"catage_proj_OFL_f");
  #ifndef NO_AD_INITIALIZE
    catage_proj_OFL_f.initialize();
  #endif
  catage_proj_OFL_m.allocate(endyr+1,endyr+15,1,nages,"catage_proj_OFL_m");
  #ifndef NO_AD_INITIALIZE
    catage_proj_OFL_m.initialize();
  #endif
  pred_catch_proj_OFL_f.allocate(endyr+1,endyr+15,"pred_catch_proj_OFL_f");
  #ifndef NO_AD_INITIALIZE
    pred_catch_proj_OFL_f.initialize();
  #endif
  pred_catch_proj_OFL_m.allocate(endyr+1,endyr+15,"pred_catch_proj_OFL_m");
  #ifndef NO_AD_INITIALIZE
    pred_catch_proj_OFL_m.initialize();
  #endif
  spawn_biom_proj.allocate(endyr+1,endyr+15,"spawn_biom_proj");
  tot_biom_proj.allocate(endyr+1,endyr+15,"tot_biom_proj");
  pred_catch_proj.allocate(endyr+1,endyr+15,"pred_catch_proj");
  pred_catch_proj_OFL.allocate(endyr+1,endyr+15,"pred_catch_proj_OFL");
  stdev_rec.allocate("stdev_rec");
  #ifndef NO_AD_INITIALIZE
  stdev_rec.initialize();
  #endif
  FOFL.allocate("FOFL");
  #ifndef NO_AD_INITIALIZE
  FOFL.initialize();
  #endif
  FABC.allocate("FABC");
  #ifndef NO_AD_INITIALIZE
  FABC.initialize();
  #endif
  FOFL2.allocate("FOFL2");
  #ifndef NO_AD_INITIALIZE
  FOFL2.initialize();
  #endif
  FABC2.allocate("FABC2");
  #ifndef NO_AD_INITIALIZE
  FABC2.initialize();
  #endif
  size_age_f.allocate(styr,endyr,1,nages,1,nlenbins,"size_age_f");
  #ifndef NO_AD_INITIALIZE
    size_age_f.initialize();
  #endif
  size_age_m.allocate(styr,endyr,1,nages,1,nlenbins,"size_age_m");
  #ifndef NO_AD_INITIALIZE
    size_age_m.initialize();
  #endif
  weight_f.allocate(styr,endyr,1,nages,"weight_f");
  #ifndef NO_AD_INITIALIZE
    weight_f.initialize();
  #endif
  weight_m.allocate(styr,endyr,1,nages,"weight_m");
  #ifndef NO_AD_INITIALIZE
    weight_m.initialize();
  #endif
  maturity.allocate(styr,endyr,1,nages,"maturity");
  #ifndef NO_AD_INITIALIZE
    maturity.initialize();
  #endif
  weight_maturity_prod_f.allocate(styr,endyr,1,nages,"weight_maturity_prod_f");
  #ifndef NO_AD_INITIALIZE
    weight_maturity_prod_f.initialize();
  #endif
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  ofstream& evalout= *pad_evalout;
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
}

dvar_vector model_parameters::SRecruit(const dvar_vector& Stmp)
{
  ofstream& evalout= *pad_evalout;
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
}

dvariable model_parameters::SRecruit(const double& Stmp)
{
  ofstream& evalout= *pad_evalout;
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
}

dvariable model_parameters::SRecruit(CONST dvariable& Stmp)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Get_Biologicals(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Get_Bzero(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Get_Selectivity(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Get_Mortality_Rates(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Get_Numbers_At_Age(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Get_Catch_at_Age(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Get_Dependent_Vars(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=styr;i<=endyr;i++)
  {
    pred_rec(i) = natage_f(i,1)+natage_m(i,1);                  // Setting up results based on estimated paramters
    tot_biom(i) = weight_f(i) * natage_f(i)+weight_m(i)*natage_m(i);          // Total biomass results
    spawn_biom(i) = weight_maturity_prod_f(i)*natage_f(i) ;                     // Spawning biomass result
  }
  avg_rec        = mean(pred_rec);
  Depletion      = spawn_biom(endyr)/spawn_biom(styr);                         // Depletion
  spbiom_trend   = spawn_biom(endyr)/spawn_biom(endyr-1);
}

void model_parameters::Get_Predicted_Values(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::compute_spr_rates(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Calc_priors(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Surv_Likelihood(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Multinomial_Likelihood(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Sel_Like(void)
{
  ofstream& evalout= *pad_evalout;
  sel_like.initialize();
  if (active(log_srv1_sel_coffs_f) ) {
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
 //       sel_like(16) +=wt_sel_dome_srv4 *square(log_srv4_sel(j)-log_srv4_sel(j+1));
 }
}

double model_parameters::round(double r)
{
  ofstream& evalout= *pad_evalout;
    return double((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5));
}

void model_parameters::Get_Population_Projection(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::Evaluate_Objective_Function(void)
{
  ofstream& evalout= *pad_evalout;
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
  obj_fun           += ssqcatch ;
  obj_fun           += sum(surv_like);
  obj_fun           += sum(age_like);
  Like = obj_fun;                     // Put here to capture the data likelihood
  obj_fun           += rec_like;
    cout<<"monitoring SSB "<<sum(elem_prod(natage_f(endyr),weight_maturity_prod_f(endyr)))<<endl;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  if (last_phase()) {
    write_projout();
    write_newproj();
    write_sarareport();
    write_HQreport();
    }
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
}

void model_parameters::write_sarareport(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::write_HQreport(void)
{
  ofstream& evalout= *pad_evalout;
   ofstream HQreport("SABLEAK_HQ.dat");
 HQreport<<"#ASSESSMENT_SUMMARY -------------------------------------------------------------------------------------------------"<<endl;
 HQreport << "#STOCK  " << endl;
 HQreport << "SABLE " << endl;
 HQreport << "#STOCK_NAME"<<endl;
 HQreport << "Sablefish - Eastern Bering Sea / Aleutian Islands / Gulf of Alaska" << endl;
  HQreport << "#REGION" <<endl;
 HQreport << "AK      " << endl;
  HQreport << "#ASMT_TYPE"<<endl;
 HQreport << "Operational  " << endl;
  HQreport << "#ASMT_YEAR"<<endl;
 HQreport << endyr << endl;
  HQreport <<"#ASMT_MONTH"<<endl;
  HQreport << "Dec  " << endl;
   HQreport <<"#TIER"<<endl;
 HQreport << "3a " << endl;
 HQreport <<"#NUM_SEXES"<<endl;
 HQreport << 2  << endl;
  HQreport <<"#NUM_FISHERIES"<<endl;
 HQreport <<2 << endl;
  HQreport <<"#REC_MULT"<<endl;
 HQreport << 1000000 << endl;
  HQreport <<"#RECAGE"<<endl;
 HQreport << 2  << endl;
  HQreport <<"#COMPLEX"<<endl;
  HQreport << "NA " << endl;
   HQreport <<"#LAST_DATA_YEAR"<<endl;
   HQreport <<  endyr  <<endl;
  HQreport <<"#ASMT_MODEL_CATEGORY"<<endl;
   HQreport <<"6 - Statistical Catch-at-Age"<<endl;
   HQreport <<"#ASMT_MODEL"<<endl;
   HQreport <<"Custom SCAA: Custom Statistical Catch-at-Age"<<endl;
   HQreport <<"#MODEL_VERSION"<<endl;
   HQreport <<"23.5"<<endl;
   HQreport <<"#ENSEMBLE"<<endl;
   HQreport <<"NA"<<endl;
   HQreport <<"#LEAD_LAB"<<endl;
   HQreport <<"AFSC ABL"<<endl;
   HQreport <<"#POC_EMAIL"<<endl;
   HQreport <<"daniel.goethel@noaa.gov"<<endl;
   HQreport <<"#REVIEW_RESULT"<<endl;
   HQreport <<"Full acceptance"<<endl;
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
   HQreport<<"F for Fully-Selected Fish"<<endl;
   HQreport <<"#F_UNIT"<<endl;
   HQreport<<"Fully-Selected F"<<endl;
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
   HQreport<<"F from 2023 asmt corresponding to specified 2022 OFL"<<endl;
   HQreport <<"#F_MSY"<<endl;
   HQreport<<FOFL<<endl;
   HQreport <<"#F_MSY_BASIS"<<endl;
   HQreport<<"F35% as proxy"<<endl;
 HQreport<<"#BIOMASS_ESTIMATES --------------------------------------------------------------------------------------------------"<<endl;
   HQreport <<"#B_YEAR"<<endl;
   HQreport <<endyr<<endl;
   HQreport <<"#B_BASIS"<<endl;
   HQreport <<"Mature Female Biomass"<<endl;
   HQreport <<"#B_UNIT"<<endl;
   HQreport <<"Kilotons"<<endl;
   HQreport <<"#BEST_B_ESTIMATE"<<endl;
   HQreport <<spawn_biom(endyr)<<endl;
   HQreport <<"#LOWER_B_ESTIMATE"<<endl;
  HQreport << "!!!!REPLACE WITH LOW 2.5% SSB FROM MCMC" <<endl;
   HQreport <<"#UPPER_B_ESTIMATE"<<endl;
 HQreport << "!!!!!REPLACE WITH high 97.5% SSB FROM MCMC" <<endl;
    HQreport <<"#ESTIMATE_METHOD"<<endl;
   HQreport <<"Credible"<<endl;
    HQreport <<"#INTERVAL_SIZE"<<endl;
   HQreport <<95<<endl;
    HQreport <<"#B_MSY"<<endl;
   HQreport << "!!!!Replace with B_35 from Summary Table"<<endl;
    HQreport <<"#B_MSY_BASIS"<<endl;
   HQreport << "B35%"<<endl;
 HQreport<<"#TIME_SERIES_ESTIMATES ----------------------------------------------------------------------------------------------"<<endl;
 HQreport <<"#FISHERYYEAR" << endl;
 HQreport << yy << endl;
  HQreport<<"#AGE"<<endl;
 HQreport << aa << endl;
 HQreport <<"#RECRUITMENT" << endl;
 HQreport  << pred_rec << endl;     
 HQreport <<"#SPAWNBIOMASS" << endl;
 HQreport  << spawn_biom << endl;  
 HQreport <<"#TOTALBIOMASS" << endl;
 HQreport  << tot_biom << endl;
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
 HQreport  << pred_catch_fish1 + pred_catch_fish3 << endl;
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
}

double model_parameters::sdnr(const dvar_vector& pred,const dvector& obs,double m)
{
  ofstream& evalout= *pad_evalout;
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred)+0.000001;
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs+0.000001-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;
}

void model_parameters::write_projout(void)
{
  ofstream& evalout= *pad_evalout;
 ofstream projout("projold.dat");
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
}

void model_parameters::write_newproj(void)
{
  ofstream& evalout= *pad_evalout;
 ofstream newproj("proj.dat");
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
}

void model_parameters::write_fullrep(void)
{
  ofstream& evalout= *pad_evalout;
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
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e-7  }");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{20000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::final_calcs()
{
 write_fullrep();
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
 gradient_structure::set_MAX_NVAR_OFFSET(1000);
   gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
   gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
   gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  arrmblsize=3900000;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
