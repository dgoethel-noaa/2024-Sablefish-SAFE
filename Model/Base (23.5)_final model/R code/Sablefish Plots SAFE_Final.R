###########################################################################################################################
# Purpose: Graphics for NOAA AK Sablefish Assessment
# Creator: Daniel Goethel and Matthew LH. Cheng (UAF-CFOS) (building on code developed by D. Hanselman)
# Date 10/17/23
###########################################################################################################################

#########################################################################
###### Required Data Files #############################################
########################################################################
# Each of these should be in the current model directory (one up from R directory)
# tem.rdat                     // current year SA report file output as rdat file
# tem_PREVYEAR.rdat                // previous year SA report file output as rdat file, "2022" should be replaced with year of prev. SA
# sable.rep                    // ADMB output report file of current year SA
# tem.par                      // ADMB output parameter file of current year SA
# tem.std                      // ADMB output standard deviation file from current year SA
# tem.ctl                      // ADMB control file of current year SA, for residual plots
# tem_PREVYEAR.ctl                 // ADMB control file of previous year SA, for residual plots
# evalout.sdat                 // MCMC ADMB output from current year SA, need to run Bayesian version of assessment to obtain
# Sablefish_Data_CurrYR.RData    // The Rdata file containing all of the assessment inputs and data pulls for the current year
# management history_PREVYR.csv   // table of the management history and catch as updated in the previous year
# hist age size coop srv.csv      // table of the historic age and length comp sample sizes for the JPN-US Coop LL survey
# alt_ABC_util.rep              // alternate data file from a run that uses a different % ABC utilization...needed if use_alt_proj_bio <- 1  
##################################################################
###### Required R Scripts  ############################################
#####################################################################
# residual_functions.r     // used to get supplemental residual plots (Pearson and OSA; developed by M. Cheng)
# newswath.r               // (NOT USED FOR MAIN PLOTS)  script used to make additional MCMC plot, old plot code; plot quantiles of timeseries (from MCMC), year indicates last time updated (does not need to be updated yearly)

###########################################################################################################################

rm(list=(ls()))

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(reshape2)
library(data.table)
library(R2admb)
library(ggsci)
library(scales)
library(gplots)
library(gridExtra)
library(cowplot)
library(janitor)
library(stringi)

###########################################################################################################################
# Install R function for residual plots, if not already installed
###########################################################################################################################

# Installation
# https://github.com/fishfollower/compResidual
# TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
# remotes::install_github("fishfollower/compResidual/compResidual", force=TRUE)
library(compResidual) 

############### INPUTS TO BE ALTERED ###############################################################
model_name_curr<-"23.5"                                  # final model name for model used in previous year/SAFE
model_name_prev<-"23.5"                                 # final model name for model used in previous year/SAFE
SA_curr_YR<-2024                                         # enter terminal year for previous stock assessment model
SA_prev_YR<-2023                                         # enter terminal year for previous stock assessment model
OFL_ABC_ACT_prev <- as.matrix(c(55084,	47146,	39032))  # the OFL, ABC, and ACT for the curr year (terminal year of SAFE)...take from the federal register final harvest spec
use_alt_proj_bio <- 1                                   # ==1 then read in a second .rep file that uses the % ABC utilized same as the projection model (only terminal year ratio, not 3 year average as calc in the .dat file)
#####################################################################################################

# Set directories
dir_R<-getwd()
setwd('..')
dir_master<-getwd()
dir_results<-paste0(dir_master,"//Results",sep='')
dir_extra<-paste0(dir_results,"//extra stuff",sep='')
dir_tables<-paste0(dir_results,"//Tables",sep='')
dir_tables_input<-paste0(dir_results,"//Tables//Inputs",sep='')
dir_tables_hist<-paste0(dir_R,"//Table Hist Inputs",sep='')


# Create directories
dir.create(dir_results)
dir.create(dir_extra)
dir.create(dir_tables)
dir.create(dir_tables_input)

# Read in data files
sab_curr <- dget(paste0(dir_master,"//tem.rdat",sep='')) 
sab_prev <- dget(paste0(dir_master,"//tem_",SA_prev_YR,".rdat",sep='')) 
sab_rep <- readLines(paste0(dir_master,"//sable.rep",sep=''))
rep_file<-readLines(paste0(dir_master,"//sable.rep",sep=''))
par_file<-readLines(paste0(dir_master,"//tem.par",sep=''))
ctl_wts <- readLines(paste0(dir_master,"//tem.ctl",sep=''))
ctl_wts_prev <- readLines(paste0(dir_master,"//tem_",SA_prev_YR,".ctl",sep=''))
std_file = readLines(paste0(dir_master,"//tem.std",sep=''))

load(paste0(dir_master,"\\Sablefish_Data_",SA_curr_YR,".RData",sep=''))        # for data 

# set indices
ages = 2:31
lengths = seq(41,99,2)
iss = 20                    # input sample size
all_years = as.numeric(row.names(sab_curr$t.series))    # number of years of assessment
curr_yr = max(all_years)
years = length(as.numeric(row.names(sab_curr$t.series)))    # number of years of assessment

# Inputs for OSA resid plots
plot_lens = lengths         # lens for plotting
plot_ages = ages            # ages for plotting
osa_lens = seq(43, 99, 2)   # for OSA resistart at 43 bc dropping first bin
osa_ages = seq(3, 31, 1)    # start at 3 bc dropping first bin


####################################################################################################################################################
# Get a bunch of random values to add to plot (nll, stock status, etc.)
####################################################################################################################################################

B40<-as.numeric(unlist(strsplit(rep_file[grep("B_40",rep_file)+1],split=" ")))                              #Get B40 from report file
mean_recruit<-as.numeric(unlist(strsplit(rep_file[grep("Mean_Recruitment",rep_file)+1],split=" ")))         #Get mean recruitment used for B40 calc from report file

B35<-B40*0.35/0.40
num_par<-as.numeric(unlist(strsplit(rep_file[grep("Number parameters estimated",rep_file)+1],split=" ")))   #Get B40 from report file
F40<-as.numeric(unlist(strsplit(rep_file[grep("F_40",rep_file)+1],split=" "))) 
F40<-F40[1]
F35<-0.40/0.35*F40
F.ABC<-as.numeric(unlist(strsplit(rep_file[grep("F_ABC",rep_file)+1],split=" "))) 
F.ABC<-F.ABC[1]
ABC<-as.numeric(unlist(strsplit(rep_file[grep("ABC for",rep_file)+1],split=" "))) 
ABC<-ABC[3] # first 2 are the F_ABC values, this grabs ABC for first projection year

max_grad<-unlist(strsplit(par_file[grep("Maximum gradient",par_file)],split=" "))
max_grad<-as.numeric(max_grad[length(max_grad)])                                                             # grab the max grad value which is the last value of first line of par file, and turn into numeric value
max_grad<-scientific(max_grad, digits = 3)
model_conv<-file.exists(paste0(dir_master,"//tem.std",sep=''))

term_SSB_tem<-as.numeric(unlist(sab_curr$t.series["spbiom"]))
years<-length(as.numeric(row.names(sab_curr$t.series)))                                                      # number of years of assessment
term_SSB<-term_SSB_tem[years]
B_B40<-term_SSB/B40
term_F_tem<-as.numeric(unlist(sab_curr$t.series["fmort"]))
term_F<-term_F_tem[years]
F_F40<-term_F/F40

likes<-as.numeric(unlist(sab_curr$likecomp))
total_lik<-likes[23]
total_lik<-round(total_lik,digits=2)

####################################################################################################################################################
# Read in MCMC results
# lot of BS in here as way that MCMC outputs are created in ADMB code (evalout.exe) is dumb; column names are off by 1, so rename to enable manipulating dataframe
####################################################################################################################################################

mcmc_mess<-dget(paste0(dir_master,"//evalout.sdat",sep=''))
mcmc_names <- names(mcmc_mess$mcmc)
mcmc <- as.data.frame(mcmc_mess$mcmc[c(-1,-length(mcmc_names))])
colnames(mcmc) <- mcmc_names[c(-(length(mcmc_names)-1),-length(mcmc_names))]
write.csv(mcmc,paste0(dir_extra,'//mcmc.csv',sep=''))

# Get 95% CIs from MCMC
uci_tem <- length(mcmc)
lci_tem <- length(mcmc)

for (i in 1:length(mcmc)){
  uci_tem[i]<-quantile(as.numeric(unlist(mcmc[[i]])) ,0.975)
  lci_tem[i]<-quantile(as.numeric(unlist(mcmc[[i]])),0.025)}

lci<-data.frame( lci = lci_tem,
                 names = mcmc_names[c(-(length(mcmc_names)-1),-length(mcmc_names))])
uci<-data.frame( uci = uci_tem,
                 names = mcmc_names[c(-(length(mcmc_names)-1),-length(mcmc_names))])
write.csv(uci,paste0(dir_extra,"//MCMC_uci.csv", sep=''))
write.csv(lci,paste0(dir_extra,"//MCMC_lci.csv", sep=''))

####################################################################################################################################################
# Set up theme for ggplot
####################################################################################################################################################

theme_reg = function() {
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 14),
        legend.background = element_blank(),
        strip.text = element_text(size = 12))
}


####################################################################################################################################################
####################################################################################################################################################
#############################################################################################################################################
#
# Plots
#
#############################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

#############################################################################################################################################
# catch by gear (3.1)
#############################################################################################################################################

catch_area_plot <- raw_catch %>%
  dplyr::mutate(Gear = dplyr::case_when(fmp_gear == "TRW" ~ "Trawl",
                                        fmp_gear == "BTR" ~ "Trawl",
                                        fmp_gear == "PTR" ~ "Trawl",
                                        fmp_gear == "NPT" ~ "Trawl",
                                        fmp_gear == "HAL" ~ "HAL",
                                        fmp_gear == "POT" ~ "Pot",
                                        fmp_gear == "OTH" ~ "HAL",
                                        fmp_gear == "JIG" ~ "HAL")) %>%   # for graphing purposes just lump limited unknown catch with HAL
  dplyr::group_by(year, area, Gear) %>%
  dplyr::summarise(catch = sum(weight_posted) / 1000) %>%
  dplyr::bind_rows(tibble(year = 1990, fixed90_complete, Gear ="HAL")) %>%
  dplyr::bind_rows(tibble(year = 1990, trawl90_complete, Gear ="Trawl")) %>%
  dplyr::arrange(year, Gear) %>%
  dplyr::filter(!is.na(Gear))

sabl_fixed_abund_subset <- sabl_fixed_abundance %>%
  dplyr::filter(fleet == "domestic", variable == "catch") %>%
  dplyr::select(year, catch= value, gear) %>%
  dplyr::mutate(Gear = recode(gear, 'llf' = 'HAL', 'tf' = 'Trawl')) 

catch_gear_historical <- catch_area_plot %>%
  dplyr::group_by(year, Gear) %>%
  dplyr::summarize(catch = sum(catch, na.rm=T)) %>%
  dplyr::bind_rows(sabl_fixed_abund_subset) %>%
  dplyr::arrange(year, Gear) %>%
  dplyr::rename(Year = year)


par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

catch_gear_plot <- catch_gear_historical %>% 
  ggplot(aes(x=Year,y=catch,fill=Gear)) +
  geom_bar(stat='identity',position="dodge", width=1) + 
  #facet_wrap(~gear) + 
  ggtitle("Catch by Gear Type")+
  ylab("Catch (kilotons)")+
  theme(legend.text=element_text(size=6))+ scale_fill_jco()+ scale_color_jco()

ggsave(paste0(dir_results , "//Fig. 3.1. Catch by gear.png"),width=8,height=5,dpi=300,plot=catch_gear_plot) 


#############################################################################################################################################
# catch by area (3.2)
#############################################################################################################################################

catch_area_hist_filter <- catch_area_hist %>%
  dplyr::select(!Grand.total) %>%
  dplyr::rename(year = Year,BS = Bering.Sea, AI = Aleu.tians, WG = Western , CG = Central, EG = Eastern) 

catch_area_plot_final <- catch_area_plot %>%
  dplyr::group_by(year, area) %>%
  dplyr::summarize(catch = 1000*sum(catch, na.rm=T)) %>%
  tidyr::pivot_wider(names_from= area, values_from= catch, names_expand = TRUE, values_fill = 0) %>%
  dplyr::select(year,BS, AI,WG,CG,EG) %>%
  dplyr::bind_rows(catch_area_hist_filter) %>%
  dplyr::arrange(year) %>%
  tidyr::pivot_longer(cols=!year, values_to= "catch", names_to ="area") %>%
  dplyr::mutate(catch = catch/1000) %>%
  dplyr::rename(Year = year, Area = area)


par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

catch_area_graph <- catch_area_plot_final %>% 
  ggplot(aes(x=Year,y=catch,fill=Area)) +
  geom_bar(stat='identity') + 
  #facet_wrap(~gear) + 
  ggtitle("Catch by NPFMC Area")+
  ylab("Catch (kilotons)")+
  scale_fill_jco()

ggsave(paste0(dir_results , "//Fig. 3.2. Catch by area.png"),width=8,height=5,dpi=300,plot=catch_area_graph) 


#############################################################################################################################################
# survey comparison plot (3.3)
#############################################################################################################################################

par(omi=c(0,0,0,0),mfrow=c(2,1),mgp=c(2.5,1,0))

LL_srvy = data.frame(year = as.numeric(row.names(sab_curr$obssrv3)), 
                     Value = as.numeric(unlist(sab_curr$obssrv3["obssrv3"])),
                     Survey = "NOAA Domestic LL Survey RPNs") %>%
          mutate(Value=Value/mean(Value))

trwl_srvy = data.frame(year = as.numeric(row.names(sab_curr$obssrv7)), 
                     Value = as.numeric(unlist(sab_curr$obssrv7["obssrv7"])),
                     Survey = "NOAA GOA Trawl Survey RPWs") %>%
          mutate(Value=Value/mean(Value))

cpue_srvy = data.frame(year = as.numeric(row.names(sab_curr$obssrv5)), 
                       Value = as.numeric(unlist(sab_curr$obssrv5["obssrv5"])),
                       Survey = "Fishery CPUE RPWs") %>%
          mutate(Value=Value/mean(Value))

surveys = LL_srvy %>% full_join(trwl_srvy) %>% full_join(cpue_srvy)

survey_plot <- ggplot(surveys, aes(x = year)) +
  geom_line(aes(y=Value,color=Survey), linewidth = 1.1) +
  geom_point(aes(y=Value,color=Survey),size=2) +
  geom_line(aes(y=1.0),colour="black",linetype=2,linewidth=.75)+ 
  labs(x="Year", y="Relative Index", title = "Index Population Trends Used in the Assessment Model") +
  theme_reg()+
  theme(legend.position = "bottom", legend.text=element_text(size=8))

ggsave(paste0(dir_results , "//Fig. 3.3. survey_trends.png"),width=8,height=5,dpi=300,plot=survey_plot) 


#############################################################################################################################################
# lls by area (3.4)
#############################################################################################################################################

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

lls_area_fig <- lls_area_rpw_apport %>%
  dplyr::select(year,area,rpn) %>%
  dplyr::rename(Year = year, Area = area) %>%
  dplyr::group_by(Year,Area) %>%
  dplyr::summarize(RPN = rpn/1000) %>%
  dplyr::mutate(Area = recode_factor(Area, 'Aleutians' = "AI" , 'Bering Sea' = 'BS', 'Western Gulf of Alaska' = 'WG',
                                     'Central Gulf of Alaska' = 'CG', 'West Yakutat' = 'WY', 'East Yakutat/Southeast' = 'EY/SE')) %>%
  ggplot(aes(x=Year,y=RPN,fill=Area)) +
  geom_bar(stat='identity') + 
  #facet_wrap(~gear) + 
  ggtitle("AFSC Longline Survey Relative Population Numbers (RPNs) by NPFMC Area")+
  ylab("RPNs (1000s)")+
  scale_fill_jco()

ggsave(paste0(dir_results , "//Fig. 3.4. LLS by area.png"),width=8,height=5,dpi=300,plot=lls_area_fig) 


#############################################################################################################################################
# lls by station (3.5)
#############################################################################################################################################

# this comes from Katy Echave every year, could probably just get her code and imbed it directly here

#############################################################################################################################################
# lls depredation (3.5)
#############################################################################################################################################

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))


### NOTE: the area_RPN tables do not include interpolated values for off-year areas, so will be less than the AK_wide index by the amount of the interpolated area 
###       RPN for whichever area is not surveyed in the current year (AI or BS)


lls_dep <- dplyr::full_join(lls_rpn_dep_alt,
                                lls_rpn_no_dep) %>%
  select(year,RPN,RPN_no_dep) %>%
  tidyr::pivot_longer(cols=!year, values_to= "rpn", names_to ="type") %>%
  dplyr::mutate(type = recode(type, 'RPN' = 'Corrected', 'RPN_no_dep' = 'Uncorrected'))

lls_dep_fig <- lls_dep %>%
  ggplot(aes(x=year,y=rpn, color=type))+
  geom_line(lwd= 1.4)+
  labs(y = "RPNs (1000s)", x = "Year", color = 'Whale Correction')+
  ggtitle("Impact of Whale Depredation Corrections on Longline Survey RPNs")+
  scale_color_jco()

ggsave(paste0(dir_results , "//Fig. 3.6. LLS depredation correction.png"),width=8,height=5,dpi=300,plot=lls_dep_fig) 


#############################################################################################################################################
# fishery depredation (3.7)
#############################################################################################################################################

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

fish_whale_dep_plot <- final_dep_area %>% 
  ggplot(aes(year, pred_dep)) +
  geom_col() + 
  facet_wrap(~factor(fmp_area,levels=c('AI','BS','WG','CG','WY','SE'))) + 
  ggtitle("Fishery Whale Depredation")+
  labs(y = "Depredation (tons)", x = "Year")+
  scale_fill_jco()

ggsave(paste0(dir_results , "//Fig. 3.7. Fishery Whale Dep by area.png"),width=8,height=5,dpi=300,plot=fish_whale_dep_plot) 


#############################################################################################################################################
# other assessment (3.8)
#############################################################################################################################################

# these are provided by DFO, NOAA WC, and ADFG assessment authors, could probably read in time series and plot all on same graph as AK assessment


#############################################################################################################################################
# Model Bridging (3.9)
############################################################################################################################################

# This is made in the mdoel bridging script 'model_comp_RUN.r'



#############################################################################################################################################
# Jitter (3.10)
############################################################################################################################################

# This is made in the jitter script 'run_jitter.r'




#############################################################################################################################################
# Likelihood components (3.11)
#############################################################################################################################################

# Likelihood component names
like_names<-c("LL Fish Age","LL Srv Age", "LL Fish Size_F","LL Fish Size_M","TRWL Fish Size_F","TRWL Fish Size_M",
              "LL Srv Size_F","LL Srv Size_M","Coop Srv Size_F","Coop Srv Size_M", "TRWL Srv Size_F","TRWL Srv Size_M",
              "LL Srv RPN","Coop Srv RPN","LL CPUE RPN", "JPN CPUE RPN", "Trawl Survey RPW","Catch",
              "Recruit_Pen","F_Pen","M_Prior")     

likecomps = sab_curr$likecomp # extract out likelihood components
likecomps_df = data.frame(nLL = likecomps, component = names(likecomps)) # dataframe for plotting

# Do some residual munging (getting rid of data not fit to and objfun, and assign groups to components)
likecomps_df = likecomps_df %>% 
  filter(nLL != 0, component != "obj.fun") %>% 
  mutate(component = like_names,
         # Setting up grouping structure for different likelihood types
         groups = case_when(
           str_detect(component, "Catch") ~ "Catch",
           str_detect(component, "Age") ~ "Age Comps",
           str_detect(component, "Size") ~ "Length Comps",
           str_detect(component, "RPW") ~ "Survey Indices",
           str_detect(component, "RPN") ~ "Survey Indices",
           str_detect(component, "CPUE") ~ "CPUE",
           str_detect(component, "Pen") ~ "Penalties",
           str_detect(component, "Prior") ~ "Penalties"),
         groups = factor(groups, levels = c("Age Comps", "Length Comps", "Survey Indices", "CPUE",
                                            "Catch", "Penalties"))) # set up order of plotting here

summary<-matrix(c(as.character(model_conv), max_grad,total_lik,num_par),ncol=4,byrow=TRUE) 
colnames(summary)<-c("Converged?","Maximum Gradient","Negative Log-Likehood","# Parameters")
summary<-as.table(summary)

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0),mai=c(0.8,0.7,0.1,0.1))

# Plot likelihood components
  like_plot<- ggplot(likecomps_df, aes(x = component, y = nLL, fill = groups)) +
  geom_col() +
  labs(x = "Likelihood Component", y = "Likelihood", fill = "") +
  scale_x_discrete(guide = guide_axis(angle = 90), limits = unique(likecomps_df$component))  +
  theme_reg() +
  theme(legend.position = c(0.85, 0.85))

  like_table <- tableGrob(summary, rows=NULL) 
  
  final_likes_plot<-grid.arrange(like_plot, like_table, nrow = 2, heights = c(2, 0.2))
  
ggsave(paste0(dir_results , "//Fig. 3.11. Like_Comps.png"),width=8,height=8,dpi=300,plot=final_likes_plot) 


#############################################################################################################################################
# Fits to indices (3.12)
#############################################################################################################################################

idx_names = c("Domestic LL Survey Relative Population Weight", "Japanese LL Survey Relative Population Weight",
              "Domestic LL Survey Relative Population Numbers", "Japanese LL Survey Relative Population Numbers",
              "Domestic Fishery CPUE Index", "Japanese Fishery CPUE Index", "GOA Trawl Survey Biomass (kt)") # index names

# extract out index data
idx_df = data.frame()
idx_list = sab_curr[str_detect(names(sab_curr), "obssrv")]

# Loop through to extract stuff out
for(i in 1:length(idx_list)) {
  # extract out index dataframe (just extracting out components from a list)
  idx_tmp = data.frame(year = as.numeric(rownames(idx_list[[i]])), obs = idx_list[[i]][[1]],
                       lci = idx_list[[i]][[2]], uci = idx_list[[i]][[3]],
                       pred = idx_list[[i]][[4]], type = idx_names[i])
  # bind together
  idx_df = rbind(idx_df, idx_tmp)
} # end i loop

# relvel by index name
idx_df$type = factor(idx_df$type, levels = idx_names)

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

# Now plot index data
index_plots <- ggplot(idx_df) +
  geom_line(mapping = aes(x = year, y = pred), linewidth = 1.25, col = "red") +
  geom_pointrange(mapping = aes(x = year, y = obs, ymin = lci, ymax = uci), 
                  alpha = 0.75, size = 0.5, pch=1,col = "blue") +
  facet_wrap(.~factor(type,levels=c("Domestic LL Survey Relative Population Numbers","Domestic LL Survey Relative Population Weight", 
                                    "Japanese LL Survey Relative Population Numbers","Japanese LL Survey Relative Population Weight",
                                    "Domestic Fishery CPUE Index", "Japanese Fishery CPUE Index", 
                                    "GOA Trawl Survey Biomass (kt)")), scales = "free",ncol=2) +
  scale_y_continuous(limits = c(0, NA)) +  # NA for no upper limit
  labs(x = "Year", y = "Index") +
  theme_reg() +
  theme(legend.position = "top") 

ggsave(paste0(dir_results , "//Fig. 3.12. Index_Fits.png"),width=13,height=13,dpi=300,plot=index_plots) 


#############################################################################################################################################
# # Composition Data Fits (3.13-3.31)
############################################################################################################################################

#########################################################################################
### Japanese LL Survey Age Compositions -------------------------------
#########################################################################################

obs_ac_srv2 = data.frame(reshape2::melt(sab_curr$oac.srv2), type = "obs")
pred_ac_srv2 = data.frame(reshape2::melt(sab_curr$eac.srv2), type = "pred")

# Put these into a dataframe
ac_srv2 = rbind(obs_ac_srv2, pred_ac_srv2)
names(ac_srv2) = c("year", "age", "prop", "type")
ac_srv2$broodYear = ac_srv2$year - ac_srv2$age                               # create brood year to track cohort over time
ac_srv2 = ac_srv2 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ac_srv2_mean = ac_srv2 %>% group_by(type, age) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Age"="age","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))

ac_srv2_mean_plot <- ggplot(ac_srv2_mean, aes(x = Age)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs( title = "Japanese Coop LL Survey Aggregated Age Compositions") +
  theme_reg()

ggsave(paste0(dir_results , "//Fig. 3.13. Aggregated Age Comps Japan LL Survey.png"),width=8,height=5,dpi=300,plot=ac_srv2_mean_plot) 


ac_srv2_plot <- ggplot() +
  geom_col(ac_srv2 %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_srv2 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Ages", y = "Proportion", title = "Japanese Coop LL Survey Age Compositions") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.14. Age Comps Japan LL Survey.png"),width=20,height=20,dpi=300,plot=ac_srv2_plot) 


#########################################################################################
### Domestic LL Survey Age Compositions -------------------------------
#########################################################################################

obs_ac_srv1 = data.frame(reshape2::melt(sab_curr$oac.srv1), type = "obs")
pred_ac_srv1 = data.frame(reshape2::melt(sab_curr$eac.srv1), type = "pred")

# Put these into a dataframe
ac_srv1 = rbind(obs_ac_srv1, pred_ac_srv1)
names(ac_srv1) = c("year", "age", "prop", "type")
ac_srv1$broodYear = ac_srv1$year - ac_srv1$age # create brood year to track cohort over time
ac_srv1 = ac_srv1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear)) # make plus group consistent color

ac_srv1_mean = ac_srv1 %>% group_by(type, age) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Age"="age","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))

ac_srv1_mean_plot <- ggplot(ac_srv1_mean, aes(x = Age)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs( title = "NOAA Domestic LL Survey Aggregated Age Compositions") +
  theme_reg()

ggsave(paste0(dir_results , "//Fig. 3.15. Aggregated Age Comps NOAA Domestic LL Survey.png"),width=8,height=5,dpi=300,plot=ac_srv1_mean_plot) 


ac_srv1_plot <- ggplot() +
  geom_col(ac_srv1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_srv1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Ages", y = "Proportion", title = "NOAA Domestic LL Survey Age Compositions") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.16. Age Comps NOAA Domestic LL Survey.png"),width=20,height=20,dpi=300,plot=ac_srv1_plot) 


#########################################################################################
### Domestic LL Fishery Age Compositions -------------------------------
#########################################################################################

obs_ac_fish1 = data.frame(reshape2::melt(sab_curr$oac.fish1), type = "obs")
pred_ac_fish1 = data.frame(reshape2::melt(sab_curr$eac.fish1), type = "pred")

# Put these into a dataframe
ac_fish1 = rbind(obs_ac_fish1, pred_ac_fish1)
names(ac_fish1) = c("year", "age", "prop", "type")
ac_fish1$broodYear = ac_fish1$year - ac_fish1$age                               # create brood year to track cohort over time
ac_fish1 = ac_fish1 %>% mutate(broodYear = ifelse(age == 31, "31", broodYear))  # make plus group consistent color

ac_fish1_mean = ac_fish1 %>% group_by(type, age) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Age"="age","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


ac_fish1_mean_plot <- ggplot(ac_fish1_mean, aes(x = Age)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs( title = "Fixed Gear Fishery Aggregated Age Compositions") +
  theme_reg()

ggsave(paste0(dir_results , "//Fig. 3.17. Aggregated Age Comps Fixed Gear Fishery.png"),width=8,height=5,dpi=300,plot=ac_fish1_mean_plot) 


ac_fish1_plot <- ggplot() +
  geom_col(ac_fish1 %>% filter(type == "obs"), mapping = aes(x = age, y = prop, fill = factor(broodYear))) +
  geom_line(ac_fish1 %>% filter(type == "pred"), mapping = aes(x = age, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Ages", y = "Proportion", title = "Domestic LL Fishery Age Compositions") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.18. Age Comps Fixed Gear Fishery.png"),width=20,height=20,dpi=300,plot=ac_fish1_plot) 


#########################################################################################
### Japanese LL Survey Length Compositions (Female) -------------------------------
#########################################################################################

obs_lc_srv2_female = data.frame(reshape2::melt(sab_curr$olc.srv2.f), type = "obs")
pred_lc_srv2_female = data.frame(reshape2::melt(sab_curr$elc.srv2.f), type = "pred")

# Put these into a dataframe
lc_srv2_female = rbind(obs_lc_srv2_female, pred_lc_srv2_female)
names(lc_srv2_female) = c("year", "length", "prop", "type")

lc_srv2_female_mean = lc_srv2_female %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


### Japanese LL Survey Length Compositions (Male) -------------------------------

obs_lc_srv2_male = data.frame(reshape2::melt(sab_curr$olc.srv2.m), type = "obs")
pred_lc_srv2_male = data.frame(reshape2::melt(sab_curr$elc.srv2.m), type = "pred")

# Put these into a dataframe
lc_srv2_male = rbind(obs_lc_srv2_male, pred_lc_srv2_male)
names(lc_srv2_male) = c("year", "length", "prop", "type")

lc_srv2_male_mean = lc_srv2_male %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


lc_srv2_female_mean_plot <- ggplot(lc_srv2_female_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)",  title = "Japanese Coop LL Survey Aggregated Length Compositions (Female)") +
  theme_reg()

lc_srv2_male_mean_plot <- ggplot(lc_srv2_male_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)", title = "Japanese Coop LL Survey Aggregated Length Compositions (Male)") +
  theme_reg()

lc_srv2_mean_plot<-grid.arrange(lc_srv2_female_mean_plot, lc_srv2_male_mean_plot, nrow = 2)

ggsave(paste0(dir_results , "//Fig. 3.19. Aggregated Length Comps Japanese LL Survey.png"),width=8,height=5,dpi=300,plot=lc_srv2_mean_plot) 


lc_srv2_female_plot <- ggplot() +
  geom_col(lc_srv2_female %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv2_female %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "Japanese Coop LL Survey Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.21. Length Comps Japanese LL Survey (Female).png"),width=20,height=20,dpi=300,plot=lc_srv2_female_plot) 


lc_srv2_male_plot <- ggplot() +
  geom_col(lc_srv2_male %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv2_male %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "Japanese Coop LL Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.20. Length Comps Japanese LL Survey (Male).png"),width=20,height=20,dpi=300,plot=lc_srv2_male_plot) 


#########################################################################################
### Domestic LL Survey Length Compositions  -------------------------------
#########################################################################################

#Female
obs_lc_srv1_female = data.frame(reshape2::melt(sab_curr$olc.srv1.f), type = "obs")
pred_lc_srv1_female = data.frame(reshape2::melt(sab_curr$elc.srv1.f), type = "pred")

# Put these into a dataframe
lc_srv1_female = rbind(obs_lc_srv1_female, pred_lc_srv1_female)
names(lc_srv1_female) = c("year", "length", "prop", "type")

lc_srv1_female_mean = lc_srv1_female %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


### Domestic LL Survey Length Compositions (Male) -------------------------------
obs_lc_srv1_male = data.frame(reshape2::melt(sab_curr$olc.srv1.m), type = "obs")
pred_lc_srv1_male = data.frame(reshape2::melt(sab_curr$elc.srv1.m), type = "pred")

# Put these into a dataframe
lc_srv1_male = rbind(obs_lc_srv1_male, pred_lc_srv1_male)
names(lc_srv1_male) = c("year", "length", "prop", "type")

lc_srv1_male_mean = lc_srv1_male %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


lc_srv1_female_mean_plot <- ggplot(lc_srv1_female_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)",  title = "NOAA Domestic LL Survey Aggregated Length Compositions (Female)") +
  theme_reg()

lc_srv1_male_mean_plot <- ggplot(lc_srv1_male_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)", title = "NOAA Domestic LL Survey Aggregated Length Compositions (Male)") +
  theme_reg()

lc_srv1_mean_plot<-grid.arrange(lc_srv1_female_mean_plot, lc_srv1_male_mean_plot, nrow = 2)

ggsave(paste0(dir_results , "//Fig. 3.22. Aggregated Length Comps NOAA LL Survey.png"),width=8,height=5,dpi=300,plot=lc_srv1_mean_plot) 


lc_srv1_female_plot <- ggplot() +
  geom_col(lc_srv1_female %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv1_female %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "NOAA Domestic LL Survey Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.24. Length Comps NOAA LL Survey (Female).png"),width=20,height=20,dpi=300,plot=lc_srv1_female_plot) 


lc_srv1_male_plot <- ggplot() +
  geom_col(lc_srv1_male %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv1_male %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "NOAA Domestic LL Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.23. Length Comps NOAA LL Survey (Male).png"),width=20,height=20,dpi=300,plot=lc_srv1_male_plot) 


#########################################################################################
### Domestic Trawl Survey Length Compositions  -------------------------------
#########################################################################################

#Female
obs_lc_srv7_female = data.frame(reshape2::melt(sab_curr$olc.srv7.f), type = "obs")
pred_lc_srv7_female = data.frame(reshape2::melt(sab_curr$elc.srv7.f), type = "pred")

# Put these into a dataframe
lc_srv7_female = rbind(obs_lc_srv7_female, pred_lc_srv7_female)
names(lc_srv7_female) = c("year", "length", "prop", "type")

lc_srv7_female_mean = lc_srv7_female %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


### Domestic Trawl Survey Length Compositions (Male) -------------------------------
obs_lc_srv7_male = data.frame(reshape2::melt(sab_curr$olc.srv7.m), type = "obs")
pred_lc_srv7_male = data.frame(reshape2::melt(sab_curr$elc.srv7.m), type = "pred")

# Put these into a dataframe
lc_srv7_male = rbind(obs_lc_srv7_male, pred_lc_srv7_male)
names(lc_srv7_male) = c("year", "length", "prop", "type")

lc_srv7_male_mean = lc_srv7_male %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


lc_srv7_female_mean_plot <- ggplot(lc_srv7_female_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)",  title = "NOAA GOA Trawl Survey Aggregated Length Compositions (Female)") +
  theme_reg()

lc_srv7_male_mean_plot <- ggplot(lc_srv7_male_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)", title = "NOAA GOA Trawl Survey Aggregated Length Compositions (Male)") +
  theme_reg()

lc_srv7_mean_plot<-grid.arrange(lc_srv7_female_mean_plot, lc_srv7_male_mean_plot, nrow = 2)

ggsave(paste0(dir_results , "//Fig. 3.25. Aggregated Length Comps NOAA Trawl Survey.png"),width=8,height=5,dpi=300,plot=lc_srv7_mean_plot) 


lc_srv7_female_plot <- ggplot() +
  geom_col(lc_srv7_female %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv7_female %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "NOAA GOA Trawl Survey Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.27. Length Comps NOAA Trawl Survey (Female).png"),width=20,height=20,dpi=300,plot=lc_srv7_female_plot) 


lc_srv7_male_plot <- ggplot() +
  geom_col(lc_srv7_male %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_srv7_male %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "NOAA GOA Trawl Survey Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.26. Length Comps NOAA Trawl Survey (Male).png"),width=20,height=20,dpi=300,plot=lc_srv7_male_plot) 


#########################################################################################
### Domestic LL Fishery Length Compositions  -------------------------------
#########################################################################################

#Female
obs_lc_fish1_female = data.frame(reshape2::melt(sab_curr$olc.fish1.f), type = "obs")
pred_lc_fish1_female = data.frame(reshape2::melt(sab_curr$elc.fish1.f), type = "pred")

# Put these into a dataframe
lc_fish1_female = rbind(obs_lc_fish1_female, pred_lc_fish1_female)
names(lc_fish1_female) = c("year", "length", "prop", "type")

lc_fish1_female_mean = lc_fish1_female %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))

### Domestic LL Fishery Length Compositions (Male) -------------------------------
obs_lc_fish1_male = data.frame(reshape2::melt(sab_curr$olc.fish1.m), type = "obs")
pred_lc_fish1_male = data.frame(reshape2::melt(sab_curr$elc.fish1.m), type = "pred")

# Put these into a dataframe
lc_fish1_male = rbind(obs_lc_fish1_male, pred_lc_fish1_male)
names(lc_fish1_male) = c("year", "length", "prop", "type")

lc_fish1_male_mean = lc_fish1_male %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


lc_fish1_female_mean_plot <- ggplot(lc_fish1_female_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)",  title = "Fixed Gear Fishery Aggregated Length Compositions (Female)") +
  theme_reg()

lc_fish1_male_mean_plot <- ggplot(lc_fish1_male_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)", title = "Fixed Gear Fishery Aggregated Length Compositions (Male)") +
  theme_reg()

lc_fish1_mean_plot<-grid.arrange(lc_fish1_female_mean_plot, lc_fish1_male_mean_plot, nrow = 2)

ggsave(paste0(dir_results , "//Fig. 3.28. Aggregated Length Comps Fixed Gear Fishery.png"),width=8,height=5,dpi=300,plot=lc_fish1_mean_plot) 


lc_fish1_female_plot <- ggplot() +
  geom_col(lc_fish1_female %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish1_female %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "Fixed Gear Fishery Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.30. Length Comps Fixed Gear Fishery (Female).png"),width=20,height=20,dpi=300,plot=lc_fish1_female_plot) 

lc_fish1_male_plot <- ggplot() +
  geom_col(lc_fish1_male %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish1_male %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "Fixed Gear Fishery Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.29. Length Comps Fixed Gear Fishery (Male).png"),width=20,height=20,dpi=300,plot=lc_fish1_male_plot) 


#########################################################################################
### Domestic Trawl Fishery Length Compositions -------------------------------
#########################################################################################

#Female
obs_lc_fish3_female = data.frame(reshape2::melt(sab_curr$olc.fish3.f), type = "obs")
pred_lc_fish3_female = data.frame(reshape2::melt(sab_curr$elc.fish3.f), type = "pred")

# Put these into a dataframe
lc_fish3_female = rbind(obs_lc_fish3_female, pred_lc_fish3_female)
names(lc_fish3_female) = c("year", "length", "prop", "type")

lc_fish3_female_mean = lc_fish3_female %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))

### Domestic Trawl Fishery Length Compositions (Male) -------------------------------
obs_lc_fish3_male = data.frame(reshape2::melt(sab_curr$olc.fish3.m), type = "obs")
pred_lc_fish3_male = data.frame(reshape2::melt(sab_curr$elc.fish3.m), type = "pred")

# Put these into a dataframe
lc_fish3_male = rbind(obs_lc_fish3_male, pred_lc_fish3_male)
names(lc_fish3_male) = c("year", "length", "prop", "type")

lc_fish3_male_mean = lc_fish3_male %>% group_by(type, length) %>% 
  summarize(UCI = quantile(prop, probs = 0.05), 
            LCI = quantile(prop, probs = 0.95),
            Agg = mean(prop)) %>%
  rename("Source"="type","Length"="length","Proportion"="Agg") %>%
  mutate(Source = recode(Source, obs="Obs",pred="Pred"))


lc_fish3_female_mean_plot <- ggplot(lc_fish3_female_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)",  title = "Trawl Fishery Aggregated Length Compositions (Female)") +
  theme_reg()

lc_fish3_male_mean_plot <- ggplot(lc_fish3_male_mean, aes(x = Length)) +
  geom_line(aes(y=Proportion, linetype=Source, color=Source), linewidth = 0.8) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, fill=Source), alpha=0.1) +
  labs(x="Length (cm)", title = "Trawl Fishery Aggregated Length Compositions (Male)") +
  theme_reg()

lc_fish3_mean_plot<-grid.arrange(lc_fish3_female_mean_plot, lc_fish3_male_mean_plot, nrow = 2)

ggsave(paste0(dir_results , "//Fig. 3.31. Aggregated Length Comps Trawl Fishery.png"),width=8,height=5,dpi=300,plot=lc_fish3_mean_plot) 


lc_fish3_female_plot <- ggplot() +
  geom_col(lc_fish3_female %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish3_female %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "Trawl Fishery Length Compositions (Female)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.33. Length Comps Trawl Fishery (Female).png"),width=20,height=20,dpi=300,plot=lc_fish3_female_plot) 


lc_fish3_male_plot <- ggplot() +
  geom_col(lc_fish3_male %>% filter(type == "obs"), mapping = aes(x = length, y = prop), fill = "darkgreen", alpha = 0.85) +
  geom_line(lc_fish3_male %>% filter(type == "pred"), mapping = aes(x = length, y = prop), linewidth = 1) +
  facet_wrap(~year, scales = "free_x", ncol = 3, dir = "h") +
  labs(x = "Length (cm)", y = "Proportion", title = "Trawl Fishery Length Compositions (Male)") +
  theme_reg() +
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.32. Length Comps Trawl Fishery (Male).png"),width=20,height=20,dpi=300,plot=lc_fish3_male_plot) 


#############################################################################################################################################
# Biomass and SSB with MCMC CIs (3.34)
############################################################################################################################################

par(omi=c(0,0,0,0),mfrow=c(2,1),mgp=c(2.5,1,0))

bio_ssb = data.frame(year = as.numeric(row.names(sab_curr$t.series)), 
                     bio = as.numeric(unlist(sab_curr$t.series["totbiom"])),
                     ssb = as.numeric(unlist(sab_curr$t.series["spbiom"]))) %>% 
          pivot_longer(!year, names_to = "type", values_to = "est")

bio_ci = data.frame(year = as.numeric(row.names(sab_curr$t.series)),
                 uci= uci %>% filter(names=="totBio") %>% select(uci), 
                 lci = lci %>% filter(names=="totBio") %>% select(lci),
                 type = 'bio')
                 
ssb_ci = data.frame(year = as.numeric(row.names(sab_curr$t.series)),
                    uci= uci %>% filter(names=="spbio") %>% select(uci), 
                    lci = lci %>% filter(names=="spbio") %>% select(lci),
                    b40 = B40,
                    b35 = B35,
                    type = 'ssb')

bio_ssb_ci <- bio_ci %>% full_join(ssb_ci) %>% right_join(bio_ssb)

ann_text_b40 <- data.frame(year=1965,type='ssb',b40=B40+.1*B40)
ann_text_b35 <- data.frame(year=1965,type='ssb',b35=B35-.1*B35)


  bio_ssb_plot <- ggplot(bio_ssb_ci, aes(x = year)) +
  geom_line(aes(y=est), linewidth = 1.1) +
  geom_line(aes(y=b40), linetype=2, col='red')+
  geom_line(aes(y=b35), linetype=1, col='black')+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.1) +
  geom_text(data = ann_text_b35,aes(y=b35),label="B35%", color='black')+
  geom_text(data = ann_text_b40, aes(y=b40),label="B40%", color='red')+
  facet_wrap(~type, scales = "free", labeller = as_labeller(c('bio'="Biomass",'ssb'="Spawning Stock Biomass"))) +
  labs(x="Year", y="Biomass (kt)", title = "Biomass and SSB with MCMC 95% Credible Intervals") +
  theme_reg()+
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.34. Biomass and SSB with MCMC CIs.png"),width=20,height=10,dpi=300,plot=bio_ssb_plot) 

#############################################################################################################################################
# Numbers at age (3.35-3.37)
############################################################################################################################################

# Females
n_at_age_f = data.frame(year = rownames(sab_curr$natage.female), sab_curr$natage.female) %>% 
  pivot_longer(!year, names_to = "age", values_to = "N") %>% 
  mutate(age = as.numeric(str_remove(age, "X")), year = as.numeric(year),
         sex = "Female") 

# Males
n_at_age_m = data.frame(year = rownames(sab_curr$natage.male),  sab_curr$natage.male) %>% 
  pivot_longer(!year, names_to = "age", values_to = "N") %>% 
  mutate(age = as.numeric(str_remove(age, "X")), year = as.numeric(year),
         sex = "Male")

n_at_age = rbind(n_at_age_f, n_at_age_m) # bind together

# Create proportions to look at age structure
n_at_age = n_at_age %>% 
  group_by(year) %>% 
  mutate(prop = N / sum(N))

# Plot numbers at age (females)
n_at_age_F_plot <- ggplot(n_at_age %>% filter(sex == "Female"), aes(x = year, y = N)) +
  geom_line(linewidth = 1) +
  facet_wrap(~age, scales = "free",ncol = 3, dir = "h") +
  theme_reg() +
  labs(x = "Year", y = "Numbers-at-Age (millions)", title="Population Numbers-at-Age Females")

ggsave(paste0(dir_results , "//Fig. 3.35. Numbers at age Females.png"),width=20,height=20,dpi=300,plot=n_at_age_F_plot) 

# Plot numbers at age (males)
n_at_age_M_plot <- ggplot(n_at_age %>% filter(sex == "Male"), aes(x = year, y = N)) +
  geom_line(linewidth = 1) +
  facet_wrap(~age, scales = "free",ncol = 3, dir = "h") +
  theme_reg() +
  labs(x = "Year", y = "Numbers-at-Age (millions)", title="Population Numbers-at-Age Males")

ggsave(paste0(dir_results , "//Fig. 3.36. Numbers at age Males.png"),width=20,height=20,dpi=300,plot=n_at_age_M_plot) 

# Plot age-structure of the population (filter to recent years)
N_age_rec_plot <- ggplot(n_at_age %>% filter(year %in% c((curr_yr-20):curr_yr)), aes(x = age, y = prop, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~year,ncol = 3, dir = "h") +
  theme_reg() +
  theme(legend.position = "top") +
  labs(x = "Age", y = "Proportion", fill = "Sex", title = "Population Proportion of Abundance by Age and Sex")

ggsave(paste0(dir_results , "//Fig. 3.37. Population Proportion at age Recent.png"),width=20,height=20,dpi=300,plot=N_age_rec_plot) 


#############################################################################################################################################
# Recruitment with MCMC CIs (3.38)
############################################################################################################################################

par(omi=c(0,0,0,0),mfrow=c(2,1),mgp=c(2.5,1,0))

rec_mean_rec<-mean(as.numeric(unlist(sab_curr$t.series["Recr"]))[20:(years-1)])         #average recruitment since 1978 year class (historic largest)

recruit = data.frame(cohort = as.numeric(row.names(sab_curr$t.series))[1:(years-1)]-2, 
                     rec = as.numeric(unlist(sab_curr$t.series["Recr"]))[1:(years-1)],
                     uci= uci %>% filter(names=="rec") %>% select(uci) %>% filter(row_number() <= n()-1), 
                     lci = lci %>% filter(names=="rec") %>% select(lci) %>% filter(row_number() <= n()-1),
                     rec_mean = rec_mean_rec)  %>% 
  mutate(mean_rec = mean(rec))

recruit_plot <- ggplot(recruit, aes(x = cohort)) +
  geom_line(aes(y=rec), linewidth = 0.9) +
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.5) +
  geom_line(aes(y=mean_rec),colour="red",linetype=2,linewidth=.75)+
  geom_line(aes(y=rec_mean),colour="black",linetype=2,linewidth=.75)+ 
  labs(x="Cohort", y="Recruitment (millions of fish)", title = "Recruitment (Age-2) with MCMC 95% Credible Intervals") +
  theme_reg()+
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.38. Recruitment with MCMC CIs.png"),width=8,height=5,dpi=300,plot=recruit_plot) 


#############################################################################################################################################
# Recruitment and SSB (3.39)
############################################################################################################################################

rec_ssb_df = data.frame(year = sab_curr$t.series$year[-c(1:2)], rec = sab_curr$t.series$Recr[-c(1:2)],
                        ssb = sab_curr$t.series$spbiom[-c(1:2)]) # removing first 2 years from time series

rec_ssb_plot <- ggplot(rec_ssb_df, aes(x = ssb, y = rec, label = year - 2)) +
  geom_text(color = "blue") +
  scale_x_continuous(limits = c(0, NA)) +  # NA for no upper limit
  theme_reg() +
  labs(y = "Age-2 Recruits (millions)", x = "SSB (kt)", title = "SSB and Subsequent Year Class Size")

ggsave(paste0(dir_results , "//Fig. 3.39. Recruitment and SSB.png"),width=8,height=5,dpi=300,plot=rec_ssb_plot) 


#############################################################################################################################################
# Recruitment Compared with Previous Year Model (3.40)
############################################################################################################################################

recruit_comp <- recruit %>% mutate(Model = paste0(SA_curr_YR," Assessment (Model ", model_name_curr,")",sep='')) %>%
                  full_join(data.frame(cohort = as.numeric(row.names(sab_curr$t.series))[1:(years-1)]-2,
                                       rec = as.numeric(unlist(sab_prev$t.series["Recr"]))[1:years-1],
                                       Model = paste0(SA_prev_YR," Assessment (Model ", model_name_prev,")",sep='')))


recruit_plot_comp <- ggplot(recruit_comp, aes(x = cohort, col=Model)) +
  geom_line(aes(y=rec), linewidth = 0.9) +
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.1,linetype=0) +
  labs(x="Cohort", y="Recruitment (millions of fish)", title = paste0("Model ", model_name_curr," Recruitment Compared to Previous SAFE")) +
  theme_reg()+
  theme(legend.position = c(.1,.9),
        legend.text = element_text(color = "black", size = 4.5),
        legend.title = element_text(color = "black", size = 5),
        legend.background = element_rect(fill="white",linetype=0))

ggsave(paste0(dir_results , "//Fig. 3.40. Recruit Compare with Prev SAFE.png"),width=8,height=5,dpi=300,plot=recruit_plot_comp) 


#############################################################################################################################################
# Selectivity (3.41)
############################################################################################################################################

selex_names = c("Derby fishery female","Derby fishery male","Trawl fishery female" ,"Trawl fishery male" ,"IFQ fishery female" ,"IFQ fishery male", "IFQ Recent fishery female" ,
                "IFQ Recent fishery male", "Domestic LL survey female","Domestic LL survey male", "Domestic LL Recent survey female","Domestic LL Recent survey male","Cooperative LL survey female" ,
                "Cooperative LL survey male" ,"GOA trawl survey female","GOA trawl survey male" ) # selectivity names

# pivot longer for plotting
selex = data.frame(sab_curr$agesel, ages = ages) 
names(selex) = c(selex_names, "ages") # rename columns for plotting

# pivot longer for plotting purposes
selex_df = selex %>%  
  pivot_longer(!ages, names_to = "type", values_to = "selex") %>% 
  mutate(type = factor(type, levels = selex_names),
         Sex = case_when( # differentiate sexes
           str_detect(type, "female") ~ "Female",
           str_detect(type, "male") ~ "Male" )) %>%
  mutate(fleet = case_when(
    str_detect(type, "Derby") ~ "Fixed Gear Fishery (1960-1995)",
    str_detect(type, "IFQ fishery") ~ "Fixed Gear Fishery (1995-2016)",
    str_detect(type, "IFQ Recent") ~ "Fixed Gear Fishery (2016+)",
    str_detect(type, "Trawl fishery") ~ "Trawl Fishery",
    str_detect(type, "Domestic LL survey") ~ "NOAA Domestic LL Survey (1990-2016))",
    str_detect(type, "Domestic LL Recent") ~ "NOAA Domestic LL Survey (2016+))",
    str_detect(type, "Cooperative LL") ~ "Japanese LL Survey",
    str_detect(type, "GOA trawl") ~ "NOAA GOA Trawl Survey"))


# plot selectivities!
sel_plot <- ggplot(selex_df, aes(x = ages, y = selex, color = Sex)) +
  geom_line(linewidth=1.1) +
  geom_point(size=2) +
  facet_wrap(~fleet) +
  labs(x = "Ages", y = "Selectivity", title = "Estimated Selectivity for Each Modeled Fleet") +
  theme_reg() +
  theme(legend.position = c(.82,.18),
        legend.text = element_text(color = "black", size = 20),
        legend.title = element_text(color = "black", size = 24))

ggsave(paste0(dir_results , "//Fig. 3.41. Selectivity.png"),width=20,height=20,dpi=300,plot=sel_plot) 

#############################################################################################################################################
# Fishing Mortality (3.42)
############################################################################################################################################

par(omi=c(0,0,0,0),mfrow=c(2,1),mgp=c(2.5,1,0))

fishing_mort = data.frame(year = as.numeric(row.names(sab_curr$t.series)), 
                     f = as.numeric(unlist(sab_curr$t.series["fmort"]))) %>%
                mutate(mean_f = mean(f))

f_plot <- ggplot(fishing_mort, aes(x = year)) +
  geom_bar(aes(y=f),stat="identity", position=position_dodge())+  
  geom_line(aes(y=f)) +
  geom_line(aes(y=mean_f),colour="red",linetype=2,linewidth=.75)+ 
  labs(x="Year", y="Fishing Mortality", title = "Fully-Selected Summary Fishing Mortality") +
  theme_reg()+
  theme(legend.position = "none")

ggsave(paste0(dir_results , "//Fig. 3.42. Fishing Mortality.png", sep=''),width=8,height=5,dpi=300,plot=f_plot) 


#############################################################################################################################################
# Management phase plane plot (3.43)
############################################################################################################################################

mean_ssb_proj1 <- mean(mcmc[,c(as.vector(which(colnames(mcmc)=="SSB.proj"))[1])])         # get the projected SSB from MCMC and take mean
mean_ssb_proj2 <- mean(mcmc[,c(as.vector(which(colnames(mcmc)=="SSB.proj"))[2])])

f_ssb_phase <- data.frame(year = all_years,
                    ssb = sab_curr$t.series[,5],
                    f = sab_curr$t.series[,6],
                    f35 = F35,
                    b35 = B35) %>%
                 full_join(data.frame(year = c(curr_yr+1,curr_yr+2),
                                 ssb = c(mean_ssb_proj1,mean_ssb_proj2),
                                 f = c(F.ABC, F.ABC),
                                 f35 = F35,
                                 b35 = B35)) %>%
                  mutate(f_stat = f/f35, b_stat=ssb/b35)


phase_plot <- ggplot(f_ssb_phase, aes(x = b_stat, y= f_stat)) +
  geom_hline(aes(yintercept=1),linetype=2,linewidth=0.5)+
  geom_vline(aes(xintercept=1),linetype=2,linewidth=0.5)+ 
  geom_segment(aes(x=0.05,y=0,xend=0.4/0.35,yend=1,color='blue'),linewidth=1.2)+
  geom_segment(aes(x=0.05,y=0,xend=0.4/0.35,yend=F.ABC/F35,color='red'),linewidth=1.2)+
  geom_segment(aes(x=0.4/0.35,y=1,xend=max(b_stat)+.1*max(b_stat),yend=1,color='blue'),linewidth=1.2)+
  geom_segment(aes(x=0.4/0.35,y=F.ABC/F35,xend=max(b_stat)+.1*max(b_stat),yend=F.ABC/F35,color='red'),linewidth=1.2)+
  geom_point(aes(x = b_stat, y = f_stat))+
  geom_path(aes(x = b_stat, y = f_stat), arrow=arrow())+
  geom_text(aes(label=year,fontface = ifelse(year == curr_yr, "bold", "italic")),size=2,nudge_x = 0.075)+
  labs(x="SSB/B35%", y="F/F35%", title = "Phase-Plane Diagram of Management Trajectory") +
  theme_reg()+
  theme(legend.position = "none")+
  xlim(0,3)

ggsave(paste0(dir_results , "//Fig. 3.43. Phase Plane.png", sep=''),width=8,height=5,dpi=300,plot=phase_plot) 

#############################################################################################################################################
# MCMC Marginal distribution (3.44)
############################################################################################################################################

# get this from MCMC folder and adnuts package

#############################################################################################################################################
# SSB Projected Density Distribution (3.45)
############################################################################################################################################

ProjectedSSB = mcmc[,c(as.vector(which(colnames(mcmc)=="SSB.proj"))[1:3])]                                      # take only the first 3 projection years values
names(ProjectedSSB)<-c(curr_yr+1,curr_yr+2,curr_yr+3)
  PSSB <- ProjectedSSB %>% pivot_longer(everything(),names_to = "Year", values_to = "Proj.SSB") 

plot_proj_ssb <- ggplot(PSSB,aes(x=Proj.SSB)) + 
  geom_density(aes(group=Year,fill=Year), alpha=0.3)+ 
  ylab("Density")+
  geom_vline(xintercept=B40,linetype=5,color="green",linewidth=1.2,alpha=0.7)+
  geom_vline(xintercept=B35,linetype=5,color="red",linewidth=1.2,alpha=0.7)+
  xlab("Projected Female Spawning Biomass (kt)")+
  geom_text(aes(x=B40+15,y=0.035),label="B40%",size=3,color="dark green")+
  geom_text(aes(x=B35-15,y=0.035),label="B35%",size=3,color="dark red")+
  labs(x="SSB (kt)", y="Probality Density", title = "Posterior Probability Distribution of Projected SSB") 

ggsave(paste0(dir_results , "//Fig. 3.45. Projected_SSB_Density.png", sep=''),width=8,height=5,dpi=300,plot=plot_proj_ssb)


#############################################################################################################################################
# Retrospective SSB (3.46)
############################################################################################################################################

# This is made in the retrospective script 'Sab_Retro.R'


#############################################################################################################################################
# Retrospective Recruitment (3.47)
############################################################################################################################################

# This is made in the retrospective script 'Sab_Retro.R'


#############################################################################################################################################
# All Model Historical Retrospective SSB (3.48)
############################################################################################################################################

# This is made in the retrospective script 'All Model Historical Retro (Fig 3.47).r'


#############################################################################################################################################
# Current Model Historical Retrospective SSB (3.49)
############################################################################################################################################

# This is made in the retrospective script 'Current Model Historical Retro (Fig 3.47).r


#############################################################################################################################################
# Profile Likelihood R0 (3.50)
############################################################################################################################################

# This is made in the retrospective script prolile_likelihood.r


#############################################################################################################################################
# Incremental Data Addition R and SSB (3.51)
############################################################################################################################################

# This is made in the data addition script model_comp_RUN.r


#############################################################################################################################################
# Index Sensitivity Analysis R and SSB (3.52)
############################################################################################################################################

# This is made in the index sensitivity run script model_comp_RUN.r


#############################################################################################################################################
# Sensitivity Runs R and SSB (3.53)
############################################################################################################################################

# This is made in the sensitivity run script model_comp_RUN.r


#############################################################################################################################################
# Contribution of Recent Yearclass to SSB (3.54)
############################################################################################################################################


#### Calculate spawners per year class 
p_mature<-sab_curr$growthmat[,5]
wt_m<-sab_curr$growthmat[,2]
wt_f<-sab_curr$growthmat[,1]

abund<-as.numeric(unlist(strsplit(rep_file[grep("N_proj_f",rep_file)+1],split=" "))[2:31])           # adjust number in strsplit by+3
mature<-abund*wt_f*p_mature
years_coh_plot<-seq(endyr-1,endyr-30,by=-1)
prop_mat<-mature/sum(mature)

cohort<-data.frame(cbind(years_coh_plot,ages,abund,p_mature,mature,prop_mat))
names(cohort)<-c("Year_Class","Age",paste0(curr_yr+1, " Female Abundance (millions)",sep=''),"Proportion Mature","SSB (kt)","Proportion SSB")

write.csv(cohort, paste0(dir_extra,"//percent_contribution_SSB.csv", sep=''))

yc_ssb_graph <- cohort %>% pivot_longer(cols = 3:6, names_to = "Source", values_to = "Value")

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

ssb_cont_plot<-ggplot(yc_ssb_graph,aes(x=Year_Class,y= Value,color= Source)) +
  geom_point(cex=1.9,pch=19)+
  facet_grid(factor(Source,levels=c('Proportion Mature',paste0(curr_yr+1, " Female Abundance (millions)",sep=''),'SSB (kt)','Proportion SSB'))~., scales="free")+
  geom_text(data=yc_ssb_graph %>% filter(Year_Class=='2014' | Year_Class=='2016' | Year_Class=='2017' |  Year_Class=='2019'), 
            aes(label= Year_Class), cex=2.,vjust=0,hjust=-0.45)+
  geom_segment( aes(xend=Year_Class,yend=0))+
  ylab("")+ 
  scale_x_continuous(breaks = seq(min(yc_ssb_graph$Year_Class),max(yc_ssb_graph$Year_Class), by = 2))+
  xlab("Year Class")+
  ggtitle(paste0("Contribution to ", curr_yr+1, " SSB by Year Class",sep=''))+ 
  theme(legend.position ="none",axis.line=element_line(),strip.text.y = element_text(size = 7),axis.text = element_text(size = 7))+ 
  scale_fill_jco()+ scale_color_jco()

ggsave(paste0(dir_results , "//Fig. 3.54. Percent Contr to SSB by YC.png"),plot=ssb_cont_plot,width=8,height=8,dpi=400)


#############################################################################################################################################
# Recruitment, SSB, Catch (3.55)
############################################################################################################################################

rec_proj <- rep(2*abund[1],times=2)        
ssb_proj1 <-as.numeric(unlist(strsplit(rep_file[grep(paste0("Female_Spawning Biomass for ",curr_yr+1,sep=''),rep_file)+1],split=" ")))
ssb_proj2 <-as.numeric(unlist(strsplit(rep_file[grep(paste0("Female_Spawning_Biomass for ",curr_yr+2,sep=''),rep_file)+1],split=" ")))

catch_proj1 <-as.numeric(unlist(strsplit(rep_file[grep(paste0("ABC for ",curr_yr+1,sep=''),rep_file)+1],split=" ")))[2]    # use the catch not F (the grep grabs 2 values, first is catch 2nd Fabc)
catch_proj2 <-as.numeric(unlist(strsplit(rep_file[grep(paste0("ABC for ",curr_yr+2,sep=''),rep_file)+1],split=" ")))[2]    # use the catch not F

proj_rec_ssb_catch = data.frame(year = c(curr_yr+1,curr_yr+2), rec = rec_proj,
                                 ssb = c(ssb_proj1,ssb_proj2), catch = c(catch_proj1,catch_proj2))

rec_ssb_catch = data.frame(year = sab_curr$t.series$year, rec = sab_curr$t.series$Recr,
                           ssb = sab_curr$t.series$spbiom, catch = sab_curr$t.series$Catch_HAL + sab_curr$t.series$Catch_TWL) %>%
                  full_join(proj_rec_ssb_catch)

write.csv(rec_ssb_catch, paste0(dir_extra,"//Recruit_SSB_Catch Timeseries.csv", sep=''))


# Plot!
rec_ssb_catch_plot <- ggplot(rec_ssb_catch) +
  geom_col(mapping = aes(x = year, y = rec, color = "Recruitment")) +
  geom_line(mapping = aes(x = year, y = ssb /3, color = "SSB"), linewidth = 1.) +
  geom_line(mapping = aes(x = year, y = catch/3 , color = "Catch"), linewidth = 1.) +
  scale_y_continuous(sec.axis = sec_axis(~.*3, name = "SSB or Catch (kt)") ) +
  scale_color_manual(values = c("Recruitment" = "lightblue3", "SSB" = "orange", 
                                "Catch" = "yellow2"),  name = "") + 
  labs(x = "Year", y = "Recruitment (millions of fish)", title="Population and Catch Trends") +
  theme_reg() + theme(legend.position = "bottom") 

ggsave(paste0(dir_results , "//Fig. 3.55. SSB, Recruitment, Catch.png"),plot=rec_ssb_catch_plot,width=8,height=5,dpi=300)


#############################################################################################################################################
############################################################################################################################################
#############################################################################################################################################
#
# Tables (mostly saved as csv, rows of which can be pasted into SAFE doc Tables directly)
#
###########################################################################################################################################
#############################################################################################################################################
############################################################################################################################################
#############################################################################################################################################
############################################################################################################################################


#############################################################################################################################################
# Main Summary Table (1st in exec summary)
############################################################################################################################################

# This is calculated by the projection model and comes out of spreadsheet Table 3.11 Sable_Projections_2023.xlsx


#############################################################################################################################################
# Apportionment Summary Table (2-6 in exec summary): Whale Depredated ABCs
############################################################################################################################################

# This is calculated in the apportionment spreadsheet All proportions_changable ABC_2023.xlsx

# uses the 5 year average survey area RPW values, the area whale depredation values, and the previous and new ABCs and OFLs

write.csv(apportionment,paste0(dir_tables_input,"//_Exec_Summ_Tables_Apportionment 5yr Average.csv",sep=''))

if(LLS_active==1){ #include terminal year survey values when there was a survey
                    lls_area_rpw_apport_out = lls_area_rpw_apport %>% filter(year == curr_yr)
}else{ # do not include terminal year survey values when there was not a survey
  lls_area_rpw_apport_out = lls_area_rpw_apport %>% filter(year == curr_yr-1)  }
  
write.csv(lls_area_rpw_apport_out,paste0(dir_tables_input,"//_Exec_Summ_Tables_Apportionment area RPWs terminal.csv",sep=''))

area_dep_out = final_dep_area %>% group_by(year) %>% filter (year >= curr_yr-3) %>% select(year, fmp_area, pred_dep) %>%
                 group_by(fmp_area) %>% mutate(ave_dep =  mean(pred_dep))

write.csv(area_dep_out,paste0(dir_tables_input,"//_Exec_Summ_Tables_Whale Depredation.csv",sep=''))


catch_term_subarea_out <- raw_catch %>%
  dplyr::mutate(Gear = dplyr::case_when(fmp_gear == "TRW" ~ "Trawl",
                                        fmp_gear == "BTR" ~ "Trawl",
                                        fmp_gear == "PTR" ~ "Trawl",
                                        fmp_gear == "NPT" ~ "Trawl",
                                        fmp_gear == "HAL" ~ "HAL",
                                        fmp_gear == "POT" ~ "Pot")) %>%
  dplyr::group_by(year, fmp_subarea, Gear) %>%
  dplyr::summarise(catch = sum(weight_posted))

catch_term_subarea_apport = catch_term_subarea_out %>% filter(year >= curr_yr-1) %>% group_by(year,fmp_subarea) %>%
                                summarise(catch = sum(catch))

write.csv(catch_term_subarea_apport,paste0(dir_tables_input,"//_Exec_Summ_Tables_Catch_fmp_subarea_term.csv",sep=''))


#############################################################################################################################################
# Apportionment Summary Table (7 in exec summary): Regional Biomass, ABC, Catch
############################################################################################################################################

# this table requires biomass apportioned to region (based on survey RPWs) along with catch and OFL/ABCs; also previous and current year values are based on the previous year models, only use new value for current year+1 and current year+2 (i.e., years that are calculating a new ABC/OFL)


if(use_alt_proj_bio == 1){   # if want to read in a different .rep file that uses an %ABC utilization that matches the projection model
  
  rep_file_ABC_alt <- readLines(paste0(dir_master,"//alt_ABC_util.rep",sep=''))
  
  bio_proj1 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("TotalBiomass for ",curr_yr+1,sep=''),rep_file_ABC_alt)+1],split=" ")))
  bio_proj2 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("TotalBiomass for ",curr_yr+2,sep=''),rep_file_ABC_alt)+1],split=" ")))
  
  
  abund<-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep("N_proj_f",rep_file_ABC_alt)+1],split=" "))[2:31])           # adjust number in strsplit by+3
  rec_proj <- rep(2*abund[1],times=2)        
  
  rec_proj1 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("TotalBiomass for ",curr_yr+1,sep=''),rep_file_ABC_alt)+1],split=" ")))
  rec_proj2 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("TotalBiomass for ",curr_yr+2,sep=''),rep_file_ABC_alt)+1],split=" ")))
  
  ssb_proj1 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("Female_Spawning Biomass for ",curr_yr+1,sep=''),rep_file_ABC_alt)+1],split=" ")))
  ssb_proj2 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("Female_Spawning_Biomass for ",curr_yr+2,sep=''),rep_file_ABC_alt)+1],split=" ")))
  
  catch_proj1 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("ABC for ",curr_yr+1,sep=''),rep_file_ABC_alt)+1],split=" ")))[2]    # use the catch not F (the grep grabs 2 values, first is catch 2nd Fabc)
  catch_proj2 <-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep(paste0("ABC for ",curr_yr+2,sep=''),rep_file_ABC_alt)+1],split=" ")))[2]    # use the catch not F
  
  proj_bio_catch = data.frame(year = c(curr_yr+1,curr_yr+2), rec = rec_proj,
                              ssb = c(ssb_proj1,ssb_proj2), catch = c(catch_proj1,catch_proj2), bio = c(bio_proj1,bio_proj2))
  
  write.csv(proj_bio_catch, paste0(dir_tables,"//_Exec_Summ_Tables_Bio_total.csv", sep=''))
  
  # now assign to region based on survey RPW proportions by region (using 5 year average proportions, same as catch apportionment)
  
  AI_prop <- term_apportionment[1]
  BS_prop <- term_apportionment[2]
  GOA_prop <- sum(term_apportionment[c(3:6)])
  
  
  growth.mat<-sab_curr$growthmat
  wt_f<-growth.mat$wt.f.block1  
  wt_m<-growth.mat$wt.m.block1
  
  n.f.proj.1<-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep("N_proj_f",rep_file_ABC_alt)+1],split=" "))[2:31]) # 1st year of projection
  n.f.proj.2<-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep("N_proj_f",rep_file_ABC_alt)+2],split=" "))[2:31]) # 2nd year of projection
  
  n.m.proj.1<-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep("N_proj_m",rep_file_ABC_alt)+1],split=" "))[2:31])
  n.m.proj.2<-as.numeric(unlist(strsplit(rep_file_ABC_alt[grep("N_proj_m",rep_file_ABC_alt)+2],split=" "))[2:31])
  
  
  age4bio.yr1_proj<-n.f.proj.1*wt_f+n.m.proj.1*wt_m
  age4bio.yr1_proj2<-age4bio.yr1_proj[-c(1,2)]      #only age-4+
  age4bio.yr1_proj.final<-sum(age4bio.yr1_proj2)
  
  age4bio.yr2_proj<-n.f.proj.2*wt_f+n.m.proj.2*wt_m
  age4bio.yr2_proj2<-age4bio.yr2_proj[-c(1,2)]      #only age-4+
  age4bio.yr2_proj.final<-sum(age4bio.yr2_proj2)
  
  
  # now assign to region based on survey proportions
  
  age4_bio_AI_proj_yr1 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+1)) %>% mutate(age_4_bio = age4bio.yr1_proj.final*AI_prop, region = "AI",
                                                                                                         age_2_bio = bio*AI_prop, region = "AI")
  age4_bio_AI_proj_yr2 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+2)) %>% mutate(age_4_bio = age4bio.yr2_proj.final*AI_prop, region = "AI",
                                                                                                         age_2_bio = bio*AI_prop, region = "AI")
  
  age4_bio_BS_proj_yr1 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+1)) %>% mutate(age_4_bio = age4bio.yr1_proj.final*BS_prop, region = "BS",
                                                                                                         age_2_bio = bio*BS_prop, region = "BS")
  age4_bio_BS_proj_yr2 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+2)) %>% mutate(age_4_bio = age4bio.yr2_proj.final*BS_prop, region = "BS",
                                                                                                         age_2_bio = bio*BS_prop, region = "BS")
  
  age4_bio_GOA_proj_yr1 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+1)) %>% mutate(age_4_bio = age4bio.yr1_proj.final*GOA_prop, region = "GOA",
                                                                                                          age_2_bio = bio*GOA_prop, region = "GOA")
  age4_bio_GOA_proj_yr2 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+2)) %>% mutate(age_4_bio = age4bio.yr2_proj.final*GOA_prop, region = "GOA",
                                                                                                          age_2_bio = bio*GOA_prop, region = "GOA")
  
  proj_age4_bio = age4_bio_AI_proj_yr1 %>% full_join(age4_bio_AI_proj_yr2) %>% full_join(age4_bio_BS_proj_yr1) %>% full_join(age4_bio_BS_proj_yr2) %>%
    full_join(age4_bio_GOA_proj_yr1) %>% full_join(age4_bio_GOA_proj_yr2) %>% rename("Pop. Total Bio." = bio)
  
  write.csv(proj_age4_bio, paste0(dir_tables,"//_Exec_Summ_Tables_Age_4_Bio_by_region.csv", sep=''))
  
}else{
  bio_proj1 <-as.numeric(unlist(strsplit(rep_file[grep(paste0("TotalBiomass for ",curr_yr+1,sep=''),rep_file)+1],split=" ")))
  bio_proj2 <-as.numeric(unlist(strsplit(rep_file[grep(paste0("TotalBiomass for ",curr_yr+2,sep=''),rep_file)+1],split=" ")))
  
  proj_bio_catch = data.frame(year = c(curr_yr+1,curr_yr+2), rec = rec_proj,
                              ssb = c(ssb_proj1,ssb_proj2), catch = c(catch_proj1,catch_proj2), bio = c(bio_proj1,bio_proj2))
  
  write.csv(proj_bio_catch, paste0(dir_tables,"//_Exec_Summ_Tables_Bio_total.csv", sep=''))
  
  # now assign to region based on survey RPW proportions by region (using 5 year average proportions, same as catch apportionment)
  
  AI_prop <- term_apportionment[1]
  BS_prop <- term_apportionment[2]
  GOA_prop <- sum(term_apportionment[c(3:6)])
  
  
  growth.mat<-sab_curr$growthmat
  wt_f<-growth.mat$wt.f.block1  
  wt_m<-growth.mat$wt.m.block1
  
  n.f.proj.1<-as.numeric(unlist(strsplit(rep_file[grep("N_proj_f",rep_file)+1],split=" "))[2:31]) # 1st year of projection
  n.f.proj.2<-as.numeric(unlist(strsplit(rep_file[grep("N_proj_f",rep_file)+2],split=" "))[2:31]) # 2nd year of projection
  
  n.m.proj.1<-as.numeric(unlist(strsplit(rep_file[grep("N_proj_m",rep_file)+1],split=" "))[2:31])
  n.m.proj.2<-as.numeric(unlist(strsplit(rep_file[grep("N_proj_m",rep_file)+2],split=" "))[2:31])
  
  
  age4bio.yr1_proj<-n.f.proj.1*wt_f+n.m.proj.1*wt_m
  age4bio.yr1_proj2<-age4bio.yr1_proj[-c(1,2)]      #only age-4+
  age4bio.yr1_proj.final<-sum(age4bio.yr1_proj2)
  
  age4bio.yr2_proj<-n.f.proj.2*wt_f+n.m.proj.2*wt_m
  age4bio.yr2_proj2<-age4bio.yr2_proj[-c(1,2)]      #only age-4+
  age4bio.yr2_proj.final<-sum(age4bio.yr2_proj2)
  
  
  # now assign to region based on survey proportions
  
  age4_bio_AI_proj_yr1 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+1)) %>% mutate(age_4_bio = age4bio.yr1_proj.final*AI_prop, region = "AI",
                                                                                                         age_2_bio = bio*AI_prop, region = "AI")
  age4_bio_AI_proj_yr2 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+2)) %>% mutate(age_4_bio = age4bio.yr2_proj.final*AI_prop, region = "AI",
                                                                                                         age_2_bio = bio*AI_prop, region = "AI")
  
  age4_bio_BS_proj_yr1 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+1)) %>% mutate(age_4_bio = age4bio.yr1_proj.final*BS_prop, region = "BS",
                                                                                                         age_2_bio = bio*BS_prop, region = "BS")
  age4_bio_BS_proj_yr2 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+2)) %>% mutate(age_4_bio = age4bio.yr2_proj.final*BS_prop, region = "BS",
                                                                                                         age_2_bio = bio*BS_prop, region = "BS")
  
  age4_bio_GOA_proj_yr1 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+1)) %>% mutate(age_4_bio = age4bio.yr1_proj.final*GOA_prop, region = "GOA",
                                                                                                          age_2_bio = bio*GOA_prop, region = "GOA")
  age4_bio_GOA_proj_yr2 <-  proj_bio_catch %>% select(year, bio) %>% filter(year==(curr_yr+2)) %>% mutate(age_4_bio = age4bio.yr2_proj.final*GOA_prop, region = "GOA",
                                                                                                          age_2_bio = bio*GOA_prop, region = "GOA")
  
  proj_age4_bio = age4_bio_AI_proj_yr1 %>% full_join(age4_bio_AI_proj_yr2) %>% full_join(age4_bio_BS_proj_yr1) %>% full_join(age4_bio_BS_proj_yr2) %>%
    full_join(age4_bio_GOA_proj_yr1) %>% full_join(age4_bio_GOA_proj_yr2) %>% rename("Pop. Total Bio." = bio)
  
  write.csv(proj_age4_bio, paste0(dir_tables,"//_Exec_Summ_Tables_Age_4_Bio_by_region.csv", sep=''))
  
} 



#############################################################################################################################################
# Final Whale Adjusted Catch Table (8 in exec summary)
############################################################################################################################################

# This is calculated in the apportionment spreadsheet All proportions_changable ABC_2023.xlsx

# TAC can be found in the harvest spec tables, usually just equals ABC, catch by area can be found in Table 3.1



#############################################################################################################################################
# Table 3.1: Catch by area and gear type
############################################################################################################################################

table.3.1_long <- raw_catch %>%
  dplyr::mutate(Gear = dplyr::case_when(fmp_gear == "TRW" ~ "Trawl",
                                        fmp_gear == "BTR" ~ "Trawl",
                                        fmp_gear == "PTR" ~ "Trawl",
                                        fmp_gear == "NPT" ~ "Trawl",
                                        fmp_gear == "HAL" ~ "HAL",
                                        fmp_gear == "POT" ~ "Pot",
                                        fmp_gear == "OTH" ~ "UNK",
                                        fmp_gear == "JIG" ~ "HAL"),
                                        Gear = replace_na(Gear, "UNK")) %>%   #assign the small amount of other catch to HAL just to document it in table
  dplyr::rename("Year" = "year") %>%
  dplyr::group_by(Year, fmp_subarea, Gear) %>%
  dplyr::summarise(catch = sum(weight_posted)) 


table_3.1_gear = catch_gear_historical  %>%
                     select(!gear) %>% full_join(non_fsh_catch_final %>% rename("Year" = "year", "catch" = "non_fsh_catch") %>%
                                                   add_column(Gear = 'Non_Comm_Catch')) %>%
                     mutate(catch = catch*1000) %>% 
                     pivot_wider(names_from = Gear, values_from = catch, values_fill = NA) %>% ungroup() %>%
                     mutate(total = rowSums(.,na.rm=T)-Year) %>% 
                     mutate(prop_trawl = Trawl/total) 
  

table_3.1_area_wide = catch_area_hist %>% pivot_longer(!Year, names_to = "fmp_subarea", values_to = "catch") %>%
                     dplyr::mutate(fmp_subarea = dplyr::case_when(fmp_subarea == "Bering.Sea" ~ "BS",
                                                                  fmp_subarea == "Aleu.tians" ~ "AI",
                                                                  fmp_subarea == "Western" ~ "WG",
                                                                  fmp_subarea == "Central" ~ "CG",
                                                                  fmp_subarea == "Eastern" ~ "EG",
                                                                  fmp_subarea == "Grand.total" ~ "Total")) %>%
  filter(fmp_subarea != "Total" & Year>=1960)



EG = table.3.1_long %>% filter(fmp_subarea %in% c("WY","SE","EY")) %>% group_by(Year) %>% 
  summarize(catch = sum(catch)) %>% mutate(fmp_subarea = "EG") %>%
  full_join(table_3.1_area_wide %>% filter(fmp_subarea == "EG")) %>%
  pivot_wider(names_from = fmp_subarea, values_from = catch, names_expand = TRUE, values_fill = 0) %>%
  arrange(Year) %>% full_join(catch_90 %>% filter(area == 'EG') %>% summarize(catch=sum(catch)) %>%
                                mutate(Year = 1990) %>% rename("EG" = "catch"))
                              


EY_SE = table.3.1_long %>% filter(fmp_subarea %in% c("SE","EY")) %>% group_by(Year) %>% 
  summarize(catch = sum(catch)) %>% mutate(fmp_subarea = "EY_SE") %>%
  pivot_wider(names_from = fmp_subarea, values_from = catch, names_expand = TRUE, values_fill = 0) %>%
  add_row(Year = 1960:1990, EY_SE = NA) %>%
  arrange(Year)


table_3.1_area = table_3.1_area_wide %>%
  pivot_wider(names_from= fmp_subarea, values_from= catch, names_expand = TRUE, values_fill = 0) %>%
  mutate(EY_SE = NA, WY = NA) %>% full_join(catch_90 %>% group_by(area) %>% summarize(catch=sum(catch)) %>%
                                              mutate(Year = 1990) %>% rename("fmp_subarea" = "area") %>%
                                              pivot_wider(names_from = fmp_subarea, values_from = catch, names_expand = TRUE, values_fill = 0)) %>%
  full_join(table.3.1_long %>% select(Year, fmp_subarea, catch) %>% group_by(Year, fmp_subarea) %>% summarize(catch = sum(catch)) %>%
            pivot_wider(names_from = fmp_subarea, values_from =catch, values_fill = NA)) %>% 
  select(-c("EY","SE","EG","EY_SE")) %>%
  full_join(full_join(EG,EY_SE))

table_3.1_num = table_3.1_area %>% full_join(table_3.1_gear) %>% rename("Total" = "total") %>%
              relocate(CG, .after = WG) %>% relocate(WY, .after = EG) %>% relocate(Total, .after = EY_SE) %>% replace(is.na(.), 0) %>%
              mutate(Percent_Trawl = prop_trawl*100) %>%  select(!prop_trawl)  

table_3.1 =  table_3.1_num %>% mutate(across(-Year, round,0),
                                      across(-Year, ~ format(., big.mark = ",", scientific = F, trim=TRUE))) %>%
              relocate("Total", .after = "Year") %>% relocate("Non_Comm_Catch", .after = "EY_SE")
                      

write.csv(table_3.1,paste0(dir_tables,"//Table 3.1. Catch_area_gear.csv",sep=''))


#############################################################################################################################################
# Table 3.2: Management Actions and Catch
############################################################################################################################################

# This is calculated from the harvest specs table and the final catch

# reads in 'management history.csv' which is just table 3.2 (as of 2023); need to provide input of OFL, ABC, ACT for current/terminal year at start of code

total_catch_table_3.2_catch = table_3.1_num %>% mutate(Total = Total-Non_Comm_Catch) %>% 
                            select(Year,Total)

OFL_ABC_ACT_prev_yr = OFL_ABC_ACT_prev %>% cbind(c("OFL","ABC","ACT")) %>% rbind(c(as.numeric(curr_yr),"Year")) %>%
                       as.data.frame() %>% rename("value" = "V1","type" = "V2") %>%
                        pivot_wider(names_from = type, values_from = (value)) %>% mutate(Year = as.numeric(Year)) %>%
                         full_join(total_catch_table_3.2_catch %>% filter(Year == curr_yr) %>% 
                                     rename("Catch"="Total")) %>% mutate_if(is.character,as.numeric)

man_hist = as.data.frame(read_csv(paste0(dir_tables_hist,"//management history_",year-1,".csv",sep=''), col_types = "d")) 

table_3.2 = man_hist %>% select(!Catch) %>% filter(Year < curr_yr) %>%
               full_join(total_catch_table_3.2_catch %>% filter(Year %in% c(1980:(curr_yr-1)))) %>% 
                          rename("Catch" = "Total", "ACT" = "TAC") %>% ungroup() %>%
               full_join(OFL_ABC_ACT_prev_yr) %>% 
               mutate(across(-Year, round,-2),
                      across(-Year, ~ format(., big.mark = ",", scientific = F, trim=TRUE))) %>%
               relocate(Catch, .before = OFL) 
              

write.csv(table_3.2,paste0(dir_tables,"//Table 3.2. Catch and Management.csv",sep=''))


#############################################################################################################################################
# Table 3.3: Data used and years
############################################################################################################################################

# This is edited by hand, just add an additional year for each data source, found in 'Table 3.3 Data used in assessment.xlsx'


#############################################################################################################################################
# Table 3.4: Compositional Data Sample Sizes
############################################################################################################################################

# this reads in 'hist age size coop srv.csv' which is just from Table 3.4 (as of 2023)....no idea where this data comes from, the actual comps were provided as proportions with no actual data

hist_jap_srvy_comp = as.data.frame(read_csv(paste0(dir_tables_hist,"//hist age size coop srv.csv",sep=''), col_types = "d")) %>%
                         rename("Japanese COOP Survey Age" = "age", "Japanese COOP Survey Length" = "length")


table_3.4 = fsh1_age_effN_table %>% rename("US Fixed Gear Fishery Age" = "n") %>% full_join(lls_age_effN_table %>% rename("NOAA LL Survey Age" = "n")) %>%
               full_join(fsh1_length_effN_table %>% group_by(year) %>% summarize(Eff_N = sum(Eff_N)) %>% rename("US Fixed Gear Fishery Length" = "Eff_N")) %>%   
               full_join(fsh2_length_effN_table %>% group_by(year) %>% summarize(Eff_N = sum(Eff_N)) %>% rename("US Trawl Gear Fishery Length" = "Eff_N")) %>%
               full_join(lls_length_effN_table_temp %>% group_by(year) %>% summarize(Eff_N = sum(samples)) %>% rename("NOAA LL Survey Length" = "Eff_N")) %>%
               full_join(ts_length_Eff_N %>% rename("NOAA GOA Trawl Survey Length" = "freq")) %>% 
               full_join(hist_jap_srvy_comp) %>% group_by(year) %>% arrange(year) %>%  
               #mutate_if(is.character,as.numeric) %>%
               mutate(across(everything(), ~ format(., big.mark = ",", scientific = F, trim=TRUE)))  %>%
               relocate('US Fixed Gear Fishery Age', .after = 'NOAA LL Survey Age') %>% relocate('Japanese COOP Survey Age', .before = 'NOAA LL Survey Age') %>%
               relocate('Japanese COOP Survey Length', .before = 'US Fixed Gear Fishery Length') %>% relocate('NOAA LL Survey Length', .before = 'US Fixed Gear Fishery Length') %>%
               relocate('NOAA GOA Trawl Survey Length', .before = 'US Fixed Gear Fishery Length') %>% rename("Year" = "year")
              

write.csv(table_3.4,paste0(dir_tables,"//Table 3.4. length and age comp sample sizes.csv",sep=''))


#############################################################################################################################################
# Table 3.5: Abundance Index Values
############################################################################################################################################


LL_srvy_RPN = lls_rpn_dep %>% filter(survey == "Domestic") %>% select (!c('rpw','rpw_se')) %>%
                              mutate(mean = mean(rpn),
                                     CV = rpn_se/mean) %>%
                     select(!c(survey,mean,rpn_se)) %>% complete(year = full_seq (c(1964:curr_yr),1)) %>%
                     rename("NOAA Domestic LL Survey RPNs" = rpn,"NOAA Domestic LL Survey RPN CV" = CV) 

LL_srvy_RPW = lls_rpn_dep %>% filter(survey == "Domestic") %>% select (!c('rpn','rpn_se')) %>%
                    mutate(mean = mean(rpw),
                            CV = rpw_se/mean) %>%
                    select(!c(survey,mean,rpw_se)) %>%  complete(year = full_seq (c(1964:curr_yr),1)) %>% rename(rpn = rpw) %>%
                     rename("NOAA Domestic LL Survey RPWs" = rpn, "NOAA Domestic LL Survey RPW CV" = CV) 

COOP_srvy_RPN = lls_rpn_dep %>% filter(survey == "Japanese/Cooperative") %>% select (!c('rpw','rpw_se')) %>%
                                mutate(mean = mean(rpn),
                                       CV = rpn_se/mean) %>%
                                select(!c(survey,mean,rpn_se)) %>% complete(year = full_seq (c(1964:curr_yr),1)) %>%
                         rename("Japanese COOP LL Survey RPNs" = rpn, "Japanese COOP LL Survey RPN CV" = CV) 

COOP_srvy_RPW = lls_rpn_dep %>% filter(survey == "Japanese/Cooperative") %>% select (!c('rpn','rpn_se')) %>%
                              mutate(mean = mean(rpw),
                                     CV = rpw_se/mean) %>%
                              select(!c(survey,mean,rpw_se)) %>% complete(year = full_seq (c(1964:curr_yr),1)) %>% rename(rpn = rpw) %>%
                         rename("Japanese COOP LL Survey RPWs" = rpn, "Japanese COOP LL Survey RPW CV" = CV) 

trwl_srvy = ts_bio %>% select (!c('LCI','UCI')) %>%
                  mutate(mean = mean(Biom),
                         CV = SE/mean) %>%
                  select(!c(mean,SE)) %>%  rename(rpn = Biom, year = YEAR) %>% complete(year = full_seq (c(1964:curr_yr),1)) %>%
                      rename("NOAA GOA Trawl Survey RPWs"=rpn, "NOAA GOA Trawl Survey RPW CV"=CV) 

cpue_srvy = CPUE_stand %>% select (!c('ll','ul')) %>%
                    mutate(mean = mean(CPUE),
                           CV = se/mean) %>%
                    select(!c(mean,se)) %>% rename(rpn = CPUE, year = Year) %>%  complete(year = full_seq (c(1964:curr_yr),1)) %>% 
                       rename("Fishery CPUE RPWs" = rpn, "Fishery CPUE RPW CV" = CV) 

JAP_cpue_srvy =  llf_cpue_jpn %>% select (!c('ll','ul')) %>%
                      mutate(mean = mean(cpue),
                             CV = se/mean) %>%
                      select(!c(mean,se)) %>% complete(year = full_seq (c(1964:curr_yr),1)) %>% rename(rpn = cpue) %>%
                       rename("Japanese Fishery CPUE RPWs" = rpn, "Japanese Fishery CPUE RPW CV" = CV) 


table_3.5 = bind_cols(COOP_srvy_RPN, LL_srvy_RPN[,-1],COOP_srvy_RPW[,-1],LL_srvy_RPW[,-1],trwl_srvy[,-1],JAP_cpue_srvy[,-1],cpue_srvy[,-1]) %>%
  mutate(across(-year, ~ format(., big.mark = ",", digits=2,  scientific = F, trim=TRUE))) %>% rename("Year" = "year") %>% arrange(Year)


write.csv(table_3.5,paste0(dir_tables,"//Table 3.5. Indices of abundance.csv",sep=''))


#############################################################################################################################################
# Table 3.6: Estimated Parameters
############################################################################################################################################

# Updated by hand directly 'Table 3.6 Estimated Parameters.xlsx'; usually just need to an additional year for dev vectors (rec, F)


#############################################################################################################################################
# Table 3.7: Francis Weights
############################################################################################################################################

# Updated by hand (_In text table D_data weights.xlsx); just paste in weights from ctl file (use the iter weights)

comp_wts_curr = as.data.frame(as.numeric(unlist(strsplit(ctl_wts[c(grep(" #wt fish1 age comp iter",ctl_wts):
                                                                     (grep(" #wt fish1 age comp iter",ctl_wts)+12))],split=" ")))) %>% 
                 drop_na() %>% rename_all(~paste0(SA_curr_YR," Model ", model_name_curr, sep='')) %>% 
                 mutate(source = c("Fixed Gear Fishery Age Comp", "NOAA Domestic LL Survey Age Comp","Japanese Coop LL Survey Age Comp",
                    "Fixed Gear Fishery Male Length Comp","Fixed Gear Fishery Female Length Comp","NOAA Domestic LL Survey Male Length Comp",
                    "NOAA Domestic LL Survey Female Length Comp", "Japenese Coop LL Survey Male Length Comp","Japanese Coop LL Survey Female Length Comp",
                    "Trawl Fishery Male Length Comp","Trawl Fishery Female Length Comp","NOAA GOA Trawl Survey Male Length Comp",
                    "NOAA GOA Trawl Survey Female Length Comp"),
                    ISS = iss)

index_wts_curr = as.data.frame(as.numeric(unlist(strsplit(ctl_wts[c(grep("#catch LL weight fish1",ctl_wts):
                                                                      (grep("#catch LL weight fish1",ctl_wts)+17))],split=" ")))) %>% 
  drop_na() %>% rename_all(~paste0(SA_curr_YR," Model ", model_name_curr, sep='')) %>% filter_all(all_vars(. != 0)) %>%
  mutate(source = c("Fixed Gear Fishery Catch", "Trawl Gear Fishery Catch","NOAA Domestic LL Survey RPN Index",
                    "Japenese Coop LL Survey RPN Index","Fixed Gear Fishery CPUE Index", "Japanese LL Fishery CPUE Index",
                    "NOAA GOA Trawl Survey Biomass Index"),
         ISS = NA)


comp_wts_prev = as.data.frame(as.numeric(unlist(strsplit(ctl_wts_prev[c(grep(" #wt fish1 age comp iter",ctl_wts_prev):
                                                                          (grep(" #wt fish1 age comp iter",ctl_wts_prev)+12))],split=" ")))) %>% 
                 drop_na() %>% rename_all(~paste0(SA_prev_YR," Model ", model_name_prev, sep='')) %>% 
                 mutate(source = c("Fixed Gear Fishery Age Comp", "NOAA Domestic LL Survey Age Comp","Japanese Coop LL Survey Age Comp",
                                   "Fixed Gear Fishery Male Length Comp","Fixed Gear Fishery Female Length Comp","NOAA Domestic LL Survey Male Length Comp",
                                   "NOAA Domestic LL Survey Female Length Comp", "Japenese Coop LL Survey Male Length Comp","Japanese Coop LL Survey Female Length Comp",
                                   "Trawl Fishery Male Length Comp","Trawl Fishery Female Length Comp","NOAA GOA Trawl Survey Male Length Comp",
                                   "NOAA GOA Trawl Survey Female Length Comp"),
                        ISS = iss)


index_wts_prev = as.data.frame(as.numeric(unlist(strsplit(ctl_wts_prev[c(grep("#catch LL weight fish1",ctl_wts_prev):
                                                                           (grep("#catch LL weight fish1",ctl_wts_prev)+17))],split=" ")))) %>% 
  drop_na() %>% rename_all(~paste0(SA_prev_YR," Model ", model_name_prev, sep='')) %>% filter_all(all_vars(. != 0)) %>%
  mutate(source = c("Fixed Gear Fishery Catch", "Trawl Gear Fishery Catch","NOAA Domestic LL Survey RPN Index",
                    "Japenese Coop LL Survey RPN Index","Fixed Gear Fishery CPUE Index", "Japanese LL Fishery CPUE Index",
                    "NOAA GOA Trawl Survey Biomass Index"),
         ISS = NA)


table_3.7 = comp_wts_curr %>% full_join(index_wts_curr) %>% full_join(comp_wts_prev %>% full_join(index_wts_prev)) %>% 
                mutate(across(-c(source), round, 3)) %>% 
                relocate('source','ISS') %>% rename("Data Source" = "source")

write.csv(table_3.7,paste0(dir_tables,"//Table 3.7. Data weights.csv",sep=''))


#############################################################################################################################################
# Table 3.8: Estimates of Recruits, Bio, SSB with MCMC CIs
############################################################################################################################################

table_3.8 = bio_ssb_ci %>% select(!c(b40,b35)) %>% filter(type == 'ssb') %>% select(!type) %>% rename("SSB" = "est") %>%
              full_join(bio_ssb_ci %>% select(!c(b40,b35)) %>% filter(type == 'bio') %>% select(!type) %>% 
                          rename("Biomass" = "est", "lci_bio" = "lci", "uci_bio" = "uci")) %>%
              full_join(recruit %>% mutate(year = cohort+2) %>% select(!c(rec_mean,mean_rec,cohort)) %>% 
                          rename("uci_rec" = "uci", "lci_rec" = "lci")) %>% mutate(across(-c(year), round, 1)) %>%
              relocate(c('rec','lci_rec','uci_rec')) %>% relocate (c("SSB","lci"), .before = "year") %>% 
              relocate (c("Biomass","lci_bio"), .before = "uci_bio") %>% relocate("year") %>%
              rename("Year" = "year", "Recruitment (millions) MLE" = "rec", "Recruitment 2.5% CI" = "lci_rec", "Recruitment 97.5% CI" = "uci_rec",
                     "SSB 2.5% CI" = "lci", "SSB 97.5% CI" = "uci", "Biomass 2.5% CI" = "lci_bio", "Biomass 97.5% CI" = "uci_bio",
                     "SSB (kt) MLE" = "SSB", "Biomass (kt) MLE" = "Biomass")

write.csv(table_3.8,paste0(dir_tables,"//Table 3.8. MLE and CIs Recruit, SSB, Bio.csv",sep=''))

#############################################################################################################################################
# Table 3.9: MCMC Estimates of key Parameters
############################################################################################################################################

q1 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("q_srv1",std_file)][2],split=" "))) %>% nth(-2),       # want transformed value not log_space so take 2nd entry in list
                mle_std = as.numeric(unlist(strsplit(std_file[grep("q_srv1",std_file)][2],split=" "))) %>% last(),
                mcmc_mean = mcmc %>% select("q1") %>% unlist() %>% as.numeric() %>% mean(na.rm = T),          
                mcmc_med = mcmc %>% select("q1") %>% unlist() %>% as.numeric() %>% median(na.rm = T),          
                mcmc_std = mcmc %>% select("q1") %>% unlist() %>% as.numeric() %>% sd(na.rm = T),          
                uci_mcmc = uci %>% filter(names=="q1") %>% select(uci),          
                lci_mcmc = lci %>% filter(names=="q1") %>% select(lci)) %>%
                transpose() %>% rename("q1" = "V1") %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

q2 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("q_srv2",std_file)][2],split=" "))) %>% nth(-2),       # want transformed value not log_space so take 2nd entry in list
                     mle_std = as.numeric(unlist(strsplit(std_file[grep("q_srv2",std_file)][2],split=" "))) %>% last(),
                     mcmc_mean = mcmc %>% select("q2") %>% unlist() %>% as.numeric() %>% mean(na.rm = T),          
                     mcmc_med = mcmc %>% select("q2") %>% unlist() %>% as.numeric() %>% median(na.rm = T),          
                     mcmc_std = mcmc %>% select("q2") %>% unlist() %>% as.numeric() %>% sd(na.rm = T),          
                     uci_mcmc = uci %>% filter(names=="q2") %>% select(uci),          
                     lci_mcmc = lci %>% filter(names=="q2") %>% select(lci)) %>% 
                     transpose() %>% rename("q2" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

q7 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("q_srv7",std_file)][2],split=" "))) %>% nth(-2),       # want transformed value not log_space so take 2nd entry in list
                mle_std = as.numeric(unlist(strsplit(std_file[grep("q_srv7",std_file)][2],split=" "))) %>% last(),
                mcmc_mean = mcmc %>% select("q7") %>% unlist() %>% as.numeric() %>% mean(na.rm = T),          
                mcmc_med = mcmc %>% select("q7") %>% unlist() %>% as.numeric() %>% median(na.rm = T),          
                mcmc_std = mcmc %>% select("q7") %>% unlist() %>% as.numeric() %>% sd(na.rm = T),          
                uci_mcmc = uci %>% filter(names=="q7") %>% select(uci),          
                lci_mcmc = lci %>% filter(names=="q7") %>% select(lci)) %>% 
                transpose() %>% rename("q7" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

M = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("M_est",std_file)],split=" "))) %>% nth(-2),       
               mle_std = as.numeric(unlist(strsplit(std_file[grep("M_est",std_file)],split=" "))) %>% last(),
               mcmc_mean = mcmc %>% select("M") %>% unlist() %>% as.numeric() %>% mean(na.rm = T),          
               mcmc_med = mcmc %>% select("M") %>% unlist() %>% as.numeric() %>% median(na.rm = T),          
               mcmc_std = mcmc %>% select("M") %>% unlist() %>% as.numeric() %>% sd(na.rm = T),          
               uci_mcmc = uci %>% filter(names=="M") %>% select(uci),          
               lci_mcmc = lci %>% filter(names=="M") %>% select(lci)) %>% 
               transpose() %>% rename("M" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

mean_r = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("log_mean_rec",std_file)],split=" "))) %>% nth(-2),       
                    mle_std = as.numeric(unlist(strsplit(std_file[grep("log_mean_rec",std_file)],split=" "))) %>% last(),
                    mcmc_mean = mcmc %>% select("LMR") %>% unlist() %>% as.numeric() %>% mean(na.rm = T),          
                    mcmc_med = mcmc %>% select("LMR") %>% unlist() %>% as.numeric() %>% median(na.rm = T),          
                    mcmc_std = mcmc %>% select("LMR") %>% unlist() %>% as.numeric() %>% sd(na.rm = T),          
                    uci_mcmc = uci %>% filter(names=="LMR") %>% select(uci),          
                    lci_mcmc = lci %>% filter(names=="LMR") %>% select(lci)) %>% 
                    transpose() %>% rename("LMR" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

sigma_r = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("sigr",std_file)],split=" "))) %>% nth(-2),       
                     mle_std = as.numeric(unlist(strsplit(std_file[grep("sigr",std_file)],split=" "))) %>% last(),
                     mcmc_mean = mcmc %>% select("sigr") %>% unlist() %>% as.numeric() %>% mean(na.rm = T),          
                     mcmc_med = mcmc %>% select("sigr") %>% unlist() %>% as.numeric() %>% median(na.rm = T),          
                     mcmc_std = mcmc %>% select("sigr") %>% unlist() %>% as.numeric() %>% sd(na.rm = T),          
                     uci_mcmc = uci %>% filter(names=="sigr") %>% select(uci),          
                     lci_mcmc = lci %>% filter(names=="sigr") %>% select(lci)) %>% 
                     transpose() %>% rename("sigr" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

Proj_Catch_tab = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("pred_catch_proj",std_file)],split=" ")[[1]])) %>% nth(-2),       # want only first projected catch value (terminal year + 1), so use [[1]]
                            mle_std = as.numeric(unlist(strsplit(std_file[grep("pred_catch_proj",std_file)],split=" ")[[1]])) %>% last(),
                            mcmc_mean = mcmc %>% select("catch.proj") %>% unlist() %>% as.numeric() %>% mean(na.rm = T),          
                            mcmc_med = mcmc %>% select("catch.proj") %>% unlist() %>% as.numeric() %>% median(na.rm = T),          
                            mcmc_std = mcmc %>% select("catch.proj") %>% unlist() %>% as.numeric() %>% sd(na.rm = T),          
                            uci_mcmc = uci %>% filter(names=="catch.proj") %>% first() %>% select(uci),          
                            lci_mcmc = lci %>% filter(names=="catch.proj") %>% first() %>% select(lci)) %>% 
                            transpose() %>% rename("catch" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

SSB_term_tab = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("ssbsd",std_file)],split=" ")[[years]])) %>% nth(-2),       # want only terminal ssb, so use [[years]]
                      mle_std = as.numeric(unlist(strsplit(std_file[grep("ssbsd",std_file)],split=" ")[[years]])) %>% last(),
                      mcmc_mean = mcmc %>% janitor::clean_names() %>% select(paste0("spbio_",years,sep='')) %>% unlist() %>% 
                        as.numeric() %>% mean(na.rm = T),          
                      mcmc_med = mcmc %>% janitor::clean_names() %>% select(paste0("spbio_",years,sep='')) %>% unlist() %>% 
                        as.numeric() %>% median(na.rm = T),          
                      mcmc_std = mcmc %>% janitor::clean_names() %>% select(paste0("spbio_",years,sep='')) %>% unlist() %>% 
                        as.numeric() %>% sd(na.rm = T),          
                      uci_mcmc = uci %>% filter(names=="spbio") %>% last() %>% select(uci),          
                      lci_mcmc = lci %>% filter(names=="spbio") %>% last() %>% select(lci)) %>% 
                      transpose() %>% rename("ssb" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))


rec_2014 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[57]])) %>% nth(-2),       # want only 2014 year class, so use [[57]] since this won't change unless change start year
                      mle_std = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[57]])) %>% last(),
                      mcmc_mean = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",57,sep='')) %>% unlist() %>% 
                        as.numeric() %>% mean(na.rm = T),          
                      mcmc_med = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",57,sep='')) %>% unlist() %>% 
                        as.numeric() %>% median(na.rm = T),          
                      mcmc_std = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",57,sep='')) %>% unlist() %>% 
                        as.numeric() %>% sd(na.rm = T),          
                      uci_mcmc = uci %>% filter(names=="rec") %>% nth(57) %>% select(uci),          
                      lci_mcmc = lci %>% filter(names=="rec") %>% nth(57) %>% select(lci)) %>% 
                      transpose() %>% rename("rec14" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

rec_2016 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[59]])) %>% nth(-2),       # want only 2016 year class, so use [[59]] since this won't change unless change start year
                      mle_std = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[59]])) %>% last(),
                      mcmc_mean = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",59,sep='')) %>% unlist() %>% 
                        as.numeric() %>% mean(na.rm = T),          
                      mcmc_med = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",59,sep='')) %>% unlist() %>% 
                        as.numeric() %>% median(na.rm = T),          
                      mcmc_std = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",59,sep='')) %>% unlist() %>% 
                        as.numeric() %>% sd(na.rm = T),          
                      uci_mcmc = uci %>% filter(names=="rec") %>% nth(59) %>% select(uci),          
                      lci_mcmc = lci %>% filter(names=="rec") %>% nth(59) %>% select(lci)) %>% 
                      transpose() %>% rename("rec16" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))


rec_2017 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[60]])) %>% nth(-2),       # want only 2017 year class, so use [[60]] since this won't change unless change start year
                      mle_std = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[60]])) %>% last(),
                      mcmc_mean = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",60,sep='')) %>% unlist() %>% 
                        as.numeric() %>% mean(na.rm = T),          
                      mcmc_med = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",60,sep='')) %>% unlist() %>% 
                        as.numeric() %>% median(na.rm = T),          
                      mcmc_std = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",60,sep='')) %>% unlist() %>% 
                        as.numeric() %>% sd(na.rm = T),          
                      uci_mcmc = uci %>% filter(names=="rec") %>% nth(60) %>% select(uci),          
                      lci_mcmc = lci %>% filter(names=="rec") %>% nth(60) %>% select(lci)) %>% 
                      transpose() %>% rename("rec17" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

rec_2018 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[61]])) %>% nth(-2),       # want only 2018 year class, so use [[61]] since this won't change unless change start year
                      mle_std = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[61]])) %>% last(),
                      mcmc_mean = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",61,sep='')) %>% unlist() %>% 
                        as.numeric() %>% mean(na.rm = T),          
                      mcmc_med = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",61,sep='')) %>% unlist() %>% 
                        as.numeric() %>% median(na.rm = T),          
                      mcmc_std = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",61,sep='')) %>% unlist() %>% 
                        as.numeric() %>% sd(na.rm = T),          
                      uci_mcmc = uci %>% filter(names=="rec") %>% nth(61) %>% select(uci),          
                      lci_mcmc = lci %>% filter(names=="rec") %>% nth(61) %>% select(lci)) %>% 
                      transpose() %>% rename("rec18" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

rec_2019 = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[62]])) %>% nth(-2),       # want only 2019 year class, so use [[62]] since this won't change unless change start year
                      mle_std = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[62]])) %>% last(),
                      mcmc_mean = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",62,sep='')) %>% unlist() %>% 
                        as.numeric() %>% mean(na.rm = T),          
                      mcmc_med = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",62,sep='')) %>% unlist() %>% 
                        as.numeric() %>% median(na.rm = T),          
                      mcmc_std = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",62,sep='')) %>% unlist() %>% 
                        as.numeric() %>% sd(na.rm = T),          
                      uci_mcmc = uci %>% filter(names=="rec") %>% nth(62) %>% select(uci),          
                      lci_mcmc = lci %>% filter(names=="rec") %>% nth(62) %>% select(lci)) %>% 
                      transpose() %>% rename("rec19" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))


rec_term = data.frame(mle = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[years-1]])) %>% nth(-2),       # want only terminal -1 rec, so use [[years-1]]
                      mle_std = as.numeric(unlist(strsplit(std_file[grep("pred_rec",std_file)],split=" ")[[years-1]])) %>% last(),
                      mcmc_mean = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",years-1,sep='')) %>% unlist() %>% 
                        as.numeric() %>% mean(na.rm = T),          
                      mcmc_med = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",years-1,sep='')) %>% unlist() %>% 
                        as.numeric() %>% median(na.rm = T),          
                      mcmc_std = mcmc %>% janitor::clean_names() %>% select(paste0("rec_",years-1,sep='')) %>% unlist() %>% 
                        as.numeric() %>% sd(na.rm = T),          
                      uci_mcmc = uci %>% filter(names=="rec") %>% nth(years-1) %>% select(uci),          
                      lci_mcmc = lci %>% filter(names=="rec") %>% nth(years-1) %>% select(lci)) %>% 
                      transpose() %>% rename("rec_term" = "V1")  %>% mutate(source = c("mle", "mle_sd","MCMC_mean","MCMC_med","MCMC_sd","MCMC_UCI","MCMC_LCI"))

table_3.9 = q1 %>% full_join(q2) %>% full_join(q7) %>% full_join(M) %>% full_join(mean_r) %>% full_join(sigma_r) %>% full_join(Proj_Catch_tab) %>%
                     full_join(SSB_term_tab) %>% full_join(rec_2014) %>% full_join(rec_2016) %>% full_join(rec_2017) %>% full_join(rec_2018) %>%
                     full_join(rec_2019) %>% full_join(rec_term) %>% 
                     rename("q_NOAA_LL_Survey" = "q1","q_Jap-Coop_LL_Survey" = "q2",
                             "q_NOAA_Trawl_Survey" = "q7","Natural Mortality" = "M", "ln_Mean_Recruitment" = "LMR", "Recruitment Variance" = "sigr",
                             "Terminal SSB (kt)" = "ssb", "Projected Catch (kt) (Terminal Year+1)" = "catch", 
                             "2014 Year Class (millions)" = "rec14", "2016 Year Class (millions)" = "rec16",
                             "2017 Year Class (millions)" = "rec17", "2018 Year Class (millions)" = "rec18", "2019 Year Class (millions)" = "rec19",
                             "Last Estimated Year Class (millions)" = "rec_term") %>%  
                      mutate(across(-c(source), round, 2)) %>%
                      relocate('source')


write.csv(table_3.9,paste0(dir_tables,"//Table 3.9. MLE and MCMC estimates key parameters.csv",sep=''))


#############################################################################################################################################
# Table 3.10: Comparison of Current and Previous Assessment estimates of Recruitment, SSB, and  Bio
############################################################################################################################################

rec_bio_ssb_prev = data.frame(Year = as.numeric(row.names(sab_prev$t.series)), 
                     bio = as.numeric(unlist(sab_prev$t.series["totbiom"])),
                     ssb = as.numeric(unlist(sab_prev$t.series["spbiom"])),
                     rec = as.numeric(unlist(sab_prev$t.series["Recr"])),
                     F = as.numeric(unlist(sab_prev$t.series["fmort"])))

table_3.10 = table_3.8 %>% select(Year, `Recruitment (millions) MLE`, `SSB (kt) MLE`, `Biomass (kt) MLE`) %>% #filter(Year < curr_yr) %>%
                      full_join(fishing_mort %>% select(year,f) %>% rename(Year = year, 'Fishing Mortality MLE' = f)) %>%
                      full_join(rec_bio_ssb_prev) %>% group_by(Year) %>% 
                      mutate(rec_diff=(((`Recruitment (millions) MLE`-rec)/`Recruitment (millions) MLE`)*100),
                             ssb_diff=(((`SSB (kt) MLE`-ssb)/`SSB (kt) MLE`)*100),
                             bio_diff=(((`Biomass (kt) MLE`-bio)/`Biomass (kt) MLE`)*100),
                             f_diff=(((`Fishing Mortality MLE` - F)/`Fishing Mortality MLE`)*100)) %>%
                      relocate('rec', .after = Year) %>% relocate('rec_diff', .before = `SSB (kt) MLE`) %>% 
                      relocate ("ssb", .after = rec_diff) %>% relocate("ssb_diff", .before = `Biomass (kt) MLE`) %>% 
                      relocate(bio, .before = `Biomass (kt) MLE`) %>% relocate(`Fishing Mortality MLE`, .after = bio_diff) %>%
                      relocate(F, .before = `Fishing Mortality MLE`) %>%
                      rename("Current SAFE Recruitment" = `Recruitment (millions) MLE`, "Previous SAFE Recruitment" = "rec",
                             "Current SAFE SSB" = `SSB (kt) MLE`, "Previous SAFE SSB" = "ssb",
                             "Current SAFE Biomass" = `Biomass (kt) MLE`, "Previous SAFE Biomass" = "bio",
                             "Recruit Difference (%)" = "rec_diff","SSB Difference (%)" = "ssb_diff", "Biomass Difference (%)" = "bio_diff",
                             "Current SAFE F" = `Fishing Mortality MLE`, "Previous SAFE F" = F,
                             "F Difference (%)" = f_diff) %>%
                      ungroup() %>% mutate(across(-Year, ~ format(., big.mark = ",", digits=2,  scientific = F, trim=TRUE))) #~round(.,digits = 2))) 


write.csv(table_3.10,paste0(dir_tables,"//Table 3.10. Comparison Recruit, SSB, Biom Current and Prev Assessment.csv",sep=''))

#############################################################################################################################################
# Table 3.11: Projection Results
############################################################################################################################################

# This is all done in the projections folder with the final table found in the  'Table 3.11 Sable_Projections_2023.xlsx' spreadsheet

#NOTE: this can and should all be automated, currently there is a thorough step-by-step readme for the projections, but others have automated
#        the process and it would prevent any errors if it were done....this is probabyl a priority for near future



#############################################################################################################################################
# Appendix Table 3B.1: Discards by area and gear  (can just pull this directly from AKFIN, prev received from Cara)
############################################################################################################################################

write.csv(discard_catch,paste0(dir_tables,"//Table 3B.1. Discards.csv",sep=''))

#############################################################################################################################################
# Appendix Table 3B.2: Bycatch of FMP species in targeted sablefish by gear  (can just pull this directly from AKFIN, prev received from Cara)
############################################################################################################################################

write.csv(incident_catch,paste0(dir_tables,"//Table 3B.2. incidental catch.csv",sep=''))

#############################################################################################################################################
# Appendix Table 3B.3: Bycatch of nontarget species in targeted sablefish by gear  (can just pull this directly from AKFIN, prev received from Cara)
############################################################################################################################################

write.csv(nontar,paste0(dir_tables,"//Table 3B.3. nontarget bycatch.csv",sep=''))

#############################################################################################################################################
# Appendix Table 3B.4: prohibited species catch by gear and region  (can just pull this directly from AKFIN, prev received from Cara)
############################################################################################################################################

write.csv(PSC_GOA,paste0(dir_tables,"//Table 3B.4a. PSC GOA.csv",sep=''))
write.csv(PSC_BSAI,paste0(dir_tables,"//Table 3B.4b. PSC BSAI.csv",sep=''))


#############################################################################################################################################
# Appendix Table 3C.1: biological inputs
############################################################################################################################################

# this is a fixed input table, no updates needed:  Table 3C.1. Biological Inputs


#############################################################################################################################################


#############################################################################################################################################
# Appendix Table 3C.2: table of parameter estimates 
############################################################################################################################################

par_est <- matrix(std_file, ncol = 1)[-1,]
par_est_temp <- matrix(nrow = length(par_est), ncol = 4)

for(y in 1:length(par_est)){
             par_est_temp[y,] <- unlist(strsplit(par_est[y],split=" ")) %>% stri_remove_empty() 
}

par_est_final <- as.data.frame(par_est_temp) %>%  slice_head(n = num_par) %>% select(-V1) %>% 
                         rename("Par. Name" = V2, "MLE Estimate" = V3, "MLE Standard Dev" = V4) %>%
                         mutate("Par. Type" = "Estimated", 
                                "Par. Purpose" = NA, "Log Space" = "Yes") %>%
                         mutate(`Par. Purpose` = case_when(`Par. Name` == "log_rec_dev" ~ "Recr. Deviation",
                                                           `Par. Name` == "log_mean_rec" ~ "Average Recruitment",
                                                           `Par. Name` == "log_sigr" ~ "Recruitment Variance",
                                                           `Par. Name` == "logm" ~ "Natural Mortality",
                                                            grepl("F_devs", `Par. Name`) ~ "F Deviation",
                                                            grepl("avg_F", `Par. Name`) ~ "Average F",
                                                            grepl("q",`Par. Name`) ~ "Catchability",
                                                            grepl("a50_fish", `Par. Name`) ~ "Fishery Age at 50% Selectivity",
                                                            grepl("a50_srv", `Par. Name`) ~ "Survey Age at 50% Selectivity",
                                                            grepl("delta_fish", `Par. Name`) ~ "Fishery Diff in 50% and 95% Selectivity",
                                                            grepl("delta_srv",`Par. Name`) ~ "Survey Diff in 50% and 95% Selectivity"))

write.csv(par_est_final,paste0(dir_tables,"//Table 3C.2. Parameter Estimates.csv",sep=''))


#############################################################################################################################################

#############################################################################################################################################
# Appendix Table 3C.3: female number at age 
############################################################################################################################################

NAA_f <- matrix(nrow=nyrs,ncol=nages+2)
NAA_m <- matrix(nrow=nyrs,ncol=nages+2)

for(y in 1:nyrs){
  NAA_f[y,]<-as.numeric(unlist(strsplit(rep_file[grep("Numbers Females",rep_file)+y],split=" "))) # 1st year of projection
  NAA_m[y,]<-as.numeric(unlist(strsplit(rep_file[grep("Numbers Males",rep_file)+y],split=" "))) # 1st year of projection
}

NAA_f = as.data.frame(NAA_f[,-2])
NAA_m = as.data.frame(NAA_m[,-2])

colnames(NAA_f) <- c("Year", "Age-2","Age-3","Age-4","Age-5","Age-6","Age-7","Age-8","Age-9","Age-10","Age-11","Age-12","Age-13","Age-14",
                     "Age-15","Age-16","Age-17","Age-18","Age-19","Age-20","Age-21","Age-22","Age-23","Age-24","Age-25","Age-26","Age-27","Age-28",
                     "Age-29","Age-30","Age-31")

colnames(NAA_m) <- c("Year", "Age-2","Age-3","Age-4","Age-5","Age-6","Age-7","Age-8","Age-9","Age-10","Age-11","Age-12","Age-13","Age-14",
                     "Age-15","Age-16","Age-17","Age-18","Age-19","Age-20","Age-21","Age-22","Age-23","Age-24","Age-25","Age-26","Age-27","Age-28",
                     "Age-29","Age-30","Age-31")


write.csv(NAA_f,paste0(dir_tables,"//Table 3C.3. NAA_Female.csv",sep=''))


#######################################################################################################################################

#############################################################################################################################################
# Appendix Table 3C.4: male number at age 
############################################################################################################################################

write.csv(NAA_m,paste0(dir_tables,"//Table 3C.4. NAA_Male.csv",sep=''))

#############################################################################################################################################


#############################################################################################################################################
############################################################################################################################################
#############################################################################################################################################
#
# Supplementary Plots (address SSC/PT and Residuals)
#
#############################################################################################################################################
############################################################################################################################################
#############################################################################################################################################
############################################################################################################################################


# Plots for comparing Age and Length comps by HAL vs. Pot gear, 
# From M. Cheng, using .csvs, so need to adjust data pull to provide comps by gear when doing the raw pull, haven't gotten to that yet so leaving this all commented out since won't run with existing inputs


# 
# # Set up ------------------------------------------------------------------
# 
# library(here)
# library(tidyverse)
# library(ggrepel)
# library(ggridges)
# 
# # theme
# theme_trevor_jordan <- function() {
#   theme_bw() +
#     theme(panel.grid.minor = element_blank(),
#           plot.background = element_rect(fill = "white", color = NA),
#           axis.text = element_text(size = 15, color = "black"),
#           strip.text = element_text(size = 19),
#           axis.title = element_text(size = 19),
#           legend.position = "top",
#           legend.title = element_text(size = 17),
#           legend.text = element_text(size = 14))
# }
# 
# 
# # Ages --------------------------------------------------------------------
# 
# # Read in age data
# age_v0 <- read.csv(file.path(dir.dat, "norpac_age_report_flattened.csv"), skip = 7) %>% 
#   filter(Gear.Description %in% c("POT OR TRAP", "LONGLINER")) # 1999 - 2020
# age_v1 <- read.csv(file.path(dir.dat, "norpac_age_report_flattened_2021_2022.csv"), skip = 6) %>% 
#   filter(Gear.Description %in% c("POT OR TRAP", "LONGLINER")) # 2021 - 2022
# 
# # Bind these
# age_df <- rbind(age_v0, age_v1)
# 
# # create + group and clean up data
# age_plot_df <- age_df %>% 
#   rename_all(tolower) %>%  # Change to lower case
#   # Filtering certai things out
#   filter(
#     gear.description %in% c("LONGLINER", "POT OR TRAP"),
#     year > 1998, 
#     !is.na(fmp.subarea), 
#     !is.na(age), 
#     age > 1, 
#     sex != "U", 
#     fmp.area != 'INSID',
#   ) %>% 
#   # Create plus group
#   mutate(age = ifelse( age > 31, 31, age)) 
# 
# # Create df plot for sample sizes
# n_age_samp_df <- age_plot_df %>%
#   group_by(year, gear.description, sex) %>% summarize(n = n()) 
# 
# # set year levels to plot year y axis from top to bottomt
# age_year_levels <- as.factor(rev(seq(1999, 2022, 1))) 
# 
# # Plot
# ggplot() +
#   ggrepel::geom_label_repel(n_age_samp_df,
#                             mapping = aes(x = 25, y = factor(year), 
#                             label = n, fill = gear.description),
#                             alpha = 0.7, size = 5, force = 0.1
#   ) + # sample size labels
#   ggridges::geom_density_ridges(age_plot_df,
#                                 mapping = aes(x = age, y = factor(year), 
#                                 fill = gear.description), alpha = 0.5
#   ) + # density ridges
#   scale_y_discrete(limits = age_year_levels) +
#   facet_wrap(~sex) +
#   theme_trevor_jordan() +
#   xlim(2, 31) +
#   labs(x = "Age", y = "Year", fill = "Gear")
# 
# 
# # Lengths -----------------------------------------------------------------
# 
# # data prior to 2020
# length_v0 <- read.csv(file.path(dir.dat, "norpac_length_report.csv"), skip = 6) 
# # data after 2020
# length_v1 <- read.csv(file.path(dir.dat, "norpac_length_report_2021_2022.csv"), skip = 6) 
# 
# # bind together
# length_df <- rbind(length_v0, length_v1)
# 
# # some quick data stuff here
# length_dat <- length_df %>% 
#   filter(!is.na(Length..cm.),
#          Performance.Description == "NO PROBLEM",
#          Gear.Description %in% c("LONGLINER", "POT OR TRAP"),
#          Sex != "U",
#          !Length..cm. < 40) %>%  # Removing < 40 length bins
#   mutate(lengthadj = ifelse(Length..cm. < 41, 41,
#                             ifelse(Length..cm. > 99, 99, Length..cm.)),
#          # Cutting and binning to what the assessment uses
#          length_bins = cut(lengthadj, breaks = seq(39.9, 99.9, 2),
#                            labels = paste(seq(41, 99, 2))))  
# 
# # more data stuff
# plot_len_df <- length_dat %>% 
#   rename_all(tolower) %>%  # Change to lower case
#   rename(length = length..cm.) %>% 
#   filter(!is.na(fmp.area), !fmp.subarea %in% c("PWSI","SEI","WOC"))  %>% 
#   group_by(year, gear.description) 
# 
# # Create df plot for sample sizes
# n_len_samp_df <- plot_len_df %>%
#   group_by(year, gear.description, sex) %>% summarize(n = n()) 
# 
# len_year_levels <- as.factor(rev(seq(1990, 2022, 1))) 
# 
# 
# # Plot!
# ggplot() +
#   geom_density_ridges(plot_len_df %>% mutate(length = ifelse(length > 99, 99, length)),
#     mapping = aes(x = length, y = factor(year), fill = gear.description), alpha = 0.5) +
#   ggrepel::geom_label_repel(
#     n_len_samp_df,
#     mapping = aes(x = 94, y = factor(year), label = n, fill = gear.description),
#     alpha = 0.7, size = 5, force = 0.1
#   ) +
#   scale_y_discrete(limits = len_year_levels) +
#   facet_grid(~sex) +
#   theme_matt() +
#   labs(x = "Length (cm)", y = "Year", fill = "Gear")



#############################################################################################################################################
# Pearson and OSA Residuals
############################################################################################################################################

# Purpose: To employ OSA residuals, and compare to Pearson residuals
# Creator: Matthew LH. Cheng
# Date 9/27/23

# Set up theme for ggplot
theme_reg = function() {
  theme_bw() +
    theme(axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(color = "black", size = 14),
          legend.text = element_text(color = "black", size = 12),
          legend.title = element_text(color = "black", size = 14),
          legend.background = element_blank(),
          strip.text = element_text(size = 12))
}

# Run function call
source(paste0(dir_R,"//residual_functions.r",sep=''))

##################################################################################
# Fishery LL Age Comps ----------------------------------------------------
##################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish1 age comp iter")]) 

# Get fishery osa_ages OSA
fish1_age_osa = get_osa_res(obs = sab_curr$oac.fish1, pred = sab_curr$eac.fish1, 
                            iss = iss, iter_wt = iter_wt, index = osa_ages, drop_bin = 1)

# Get fishery osa_ages Pearson
fish1_age_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                    obs = sab_curr$oac.fish1, pred = sab_curr$eac.fish1)

# Plot OSA and Pearson Comparison
fish1_age_osa_plot = res_plot(data = fish1_age_osa[[1]], 
                    comp_type = "Age", res_type = "OSA Residuals for Fixed Gear Fishery Age Compositions (Remove Age-2)",
                    ymin = min(plot_ages))

fish1_age_pearson_plot = res_plot(data = fish1_age_pearson, comp_type = "Age", 
                        res_type = "Pearson Residuals for Fixed Gear Fishery Age Compositions",
                        ymin = min(plot_ages))
  
fish1_age_resid_comb_plot = plot_grid(fish1_age_osa_plot , fish1_age_pearson_plot, ncol = 1) 

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Fixed Gear Fishery Age Comp.png"),plot=fish1_age_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Fixed Gear Fishery Age Comp.png"),plot=fish1_age_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.3. Resid Comparison Fixed Gear Fishery Age Comp.png"),plot=fish1_age_resid_comb_plot,width=12,height=10,dpi=300)


##################################################################################
# Domestic Survey LL Age Comps ----------------------------------------------------
##################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv1 age comp iter")]) 

# Get domestic survey osa_ages OSA
srv1_age_osa = get_osa_res(obs = sab_curr$oac.srv1, pred = sab_curr$eac.srv1, 
                           iss = iss, iter_wt = iter_wt, index = osa_ages, drop_bin = 1)

# Get domestic survey osa_ages Pearson
srv1_age_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                   obs = sab_curr$oac.srv1, pred = sab_curr$eac.srv1)

# Plot OSA and Pearson Comparison
srv1_age_osa_plot = res_plot(data = srv1_age_osa[[1]], comp_type = "Age", 
                    res_type = "OSA Residuals for NOAA Domestic LL Survey Age Compositions (Remove Age-2)",
                    ymin = min(plot_ages))

srv1_age_pearson_plot = res_plot(data = srv1_age_pearson, comp_type = "Age",
                        res_type = "Pearson Residuals for NOAA Domestic LL Survey Age Compositions",
                        ymin = min(plot_ages))

srv1_age_comb_plot = plot_grid(srv1_age_osa_plot , srv1_age_pearson_plot, ncol = 1) 


ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Dom LL Survey Age Comp.png"),plot=srv1_age_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Dom LL Survey Age Comp.png"),plot=srv1_age_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.2. Resid Comparison Dom LL Survey Age Comp.png"),plot=srv1_age_comb_plot,width=12,height=10,dpi=300)


##################################################################################
# Japanese Survey LL Age Comps ----------------------------------------------------
##################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv2 age comp iter")]) 

# Get Japanese Survey osa_ages OSA
srv2_age_osa = get_osa_res(obs = sab_curr$oac.srv2, pred = sab_curr$eac.srv2, 
                           iss = iss, iter_wt = iter_wt, index = osa_ages, drop_bin = 1)

# Get Japanese Survey osa_ages Pearson
srv2_age_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                   obs = sab_curr$oac.srv2, pred = sab_curr$eac.srv2)

# Plot OSA and Pearson Comparison
srv2_age_osa_plot = res_plot(data = srv2_age_osa[[1]], comp_type = "Age", 
                    res_type = "OSA Residuals for Japanese LL Survey Age Compositions (Remove Age-2)",
                    ymin = min(plot_ages))

srv2_age_pearson_plot = res_plot(data = srv2_age_pearson, comp_type = "Age", 
                        res_type = "Pearson Residuals for Japanese LL Survey Age Compositions",
                        ymin = min(plot_ages))

srv2_age_comb_plot = plot_grid(srv2_age_osa_plot , srv2_age_pearson_plot, ncol = 1) 

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Japanese LL Survey Age Comp.png"),plot=srv2_age_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Japanese LL Survey Age Comp.png"),plot=srv2_age_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.1. Resid Comparison Japanese LL Survey Age Comp.png"),plot=srv2_age_comb_plot,width=12,height=10,dpi=300)


##################################################################################
# Domestic LL Survey (Length Females) -------------------------------------------
##################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv1 size comp female iter")]) 

# Get Survey length females residuals (OSA)
srv1_fem_len_osa = get_osa_res(obs = sab_curr$olc.srv1.f, pred = sab_curr$elc.srv1.f, 
                               iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length females residuals (Pearson)
srv1_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.srv1.f, pred = sab_curr$elc.srv1.f)

# Plot OSA and Pearson Comparison
srv1_fem_len_osa_plot = res_plot(data = srv1_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for NOAA Domestic LL Survey Female Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

srv1_fem_len_pearson_plot = res_plot(data = srv1_fem_len_pearson,
                        comp_type = "Length", res_type = "Pearson Residuals for NOAA Domestic LL Survey Female Length Compositions",
                        ymin = min(plot_lens))

srv1_fem_len_comb_plot = plot_grid(srv1_fem_len_osa_plot, srv1_fem_len_pearson_plot, ncol = 1)


ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Dom LL Survey Length Comp Female.png"),plot=srv1_fem_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Dom LL Survey Length Comp Female.png"),plot=srv1_fem_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.7. Resid Comparison Dom LL Survey Length Comp Female.png"),plot=srv1_fem_len_comb_plot,width=12,height=10,dpi=300)


##################################################################################
# Domestic LL Survey (Length Males) -------------------------------------------
##################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv1 size comp male iter")]) 

# Get Survey length males residuals (OSA)
srv1_male_len_osa = get_osa_res(obs = sab_curr$olc.srv1.m, pred = sab_curr$elc.srv1.m, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length males residuals (Pearson)
srv1_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.srv1.m, pred = sab_curr$elc.srv1.m)

# Plot OSA and Pearson Comparison
srv1_male_len_osa_plot = res_plot(data = srv1_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for NOAA Domestic LL Survey Male Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

srv1_male_len_pearson_plot = res_plot(data = srv1_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for NOAA Domestic LL Survey Male Length Compositions",
                        ymin = min(plot_lens))

srv1_male_len_comb_plot = plot_grid(srv1_male_len_osa_plot, srv1_male_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Dom LL Survey Length Comp Male.png"),plot=srv1_male_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Dom LL Survey Length Comp Male.png"),plot=srv1_male_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.6. Resid Comparison Dom LL Survey Length Comp Male.png"),plot=srv1_male_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Japanese LL Survey (Length Females) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv2 size comp female iter")]) 

# Get Survey length females residuals (OSA)
srv2_fem_len_osa = get_osa_res(obs = sab_curr$olc.srv2.f, pred = sab_curr$elc.srv2.f, 
                               iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length females residuals (Pearson)
srv2_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.srv2.f, pred = sab_curr$elc.srv2.f)

# Plot OSA and Pearson Comparison
srv2_fem_len_osa_plot = res_plot(data = srv2_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for Japanese LL Survey Female Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

srv2_fem_len_pearson_plot = res_plot(data = srv2_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for Japanese LL Survey Female Length Compositions",
                        ymin = min(plot_lens))

srv2_fem_len_comb_plot = plot_grid(srv2_fem_len_osa_plot, srv2_fem_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Japanese LL Survey Length Comp Female.png"),plot=srv2_fem_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Japanese LL Survey Length Comp Female.png"),plot=srv2_fem_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.5. Resid Comparison Japanese LL Survey Length Comp Female.png"),plot=srv2_fem_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Japanese LL Survey (Length Males) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt surv2 size comp male iter")]) 

# Get Survey length males residuals (OSA)
srv2_male_len_osa = get_osa_res(obs = sab_curr$olc.srv2.m, pred = sab_curr$elc.srv2.m, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length males residuals (Pearson)
srv2_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.srv2.m, pred = sab_curr$elc.srv2.m)

# Plot OSA and Pearson Comparison
srv2_male_len_osa_plot = res_plot(data = srv2_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for Japanese LL Survey Male Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

srv2_male_len_pearson_plot = res_plot(data = srv2_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for Japanese LL Survey Male Length Compositions",
                        ymin = min(plot_lens))

srv2_male_len_comb_plot = plot_grid(srv2_male_len_osa_plot, srv2_male_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Japanese LL Survey Length Comp Male.png"),plot=srv2_male_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Japanese LL Survey Length Comp Male.png"),plot=srv2_male_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.4. Resid Comparison Japanese LL Survey Length Comp Male.png"),plot=srv2_male_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Domestic LL Fishery (Length Females) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish1 size comp female iter")]) 

# Get fishery length females residuals (OSA)
fish1_fem_len_osa = get_osa_res(obs = sab_curr$olc.fish1.f, pred = sab_curr$elc.fish1.f, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length females residuals (Pearson)
fish1_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.fish1.f, pred = sab_curr$elc.fish1.f)

# Plot OSA and Pearson Comparison
fish1_fem_len_osa_plot = res_plot(data = fish1_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for Fixed Gear Fishery Female Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

fish1_fem_len_pearson_plot = res_plot(data = fish1_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for Fixed Gear Fishery Female Length Compositions",
                        ymin = min(plot_lens))

fish1_fem_len_comb_plot = plot_grid(fish1_fem_len_osa_plot, fish1_fem_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Fixed Gear Fishery Length Comp Female.png"),plot=fish1_fem_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Fixed Gear Fishery Length Comp Female.png"),plot=fish1_fem_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.11. Resid Comparison Fixed Gear Fishery Length Comp Female.png"),plot=fish1_fem_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Domestic LL Fishery (Length Males) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish1 size comp male iter")]) 

# Get fishery length males residuals (OSA)
fish1_male_len_osa = get_osa_res(obs = sab_curr$olc.fish1.m, pred = sab_curr$elc.fish1.m, 
                                 iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length males residuals (Pearson)
fish1_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                         obs = sab_curr$olc.fish1.m, pred = sab_curr$elc.fish1.m)

# Plot OSA and Pearson Comparison
fish1_male_len_osa_plot = res_plot(data = fish1_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for Fixed Gear Fishery Male Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

fish1_male_len_pearson_plot = res_plot(data = fish1_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for Fixed Gear Fishery Male Length Compositions",
                        ymin = min(plot_lens))

fish1_male_len_comb_plot = plot_grid(fish1_male_len_osa_plot, fish1_male_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Fixed Gear Fishery Length Comp Male.png"),plot=fish1_male_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Fixed Gear Fishery Length Comp Male.png"),plot=fish1_male_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.10. Resid Comparison Fixed Gear Fishery Length Comp Male.png"),plot=fish1_male_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Domestic Trawl Fishery (Length Females) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish3 size comp female iter")]) 

# Get fishery length females residuals (OSA)
fish3_fem_len_osa = get_osa_res(obs = sab_curr$olc.fish3.f, pred = sab_curr$elc.fish3.f, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length females residuals (Pearson)
fish3_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.fish3.f, pred = sab_curr$elc.fish3.f)

# Plot OSA and Pearson Comparison
fish3_fem_len_osa_plot = res_plot(data = fish3_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for Trawl Fishery Female Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

fish3_fem_len_pearson_plot = res_plot(data = fish3_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for Trawl Fishery Female Length Compositions",
                        ymin = min(plot_lens))

fish3_fem_len_comb_plot = plot_grid(fish3_fem_len_osa_plot, fish3_fem_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Trawl Fishery Length Comp Female.png"),plot=fish3_fem_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Trawl Fishery Length Comp Female.png"),plot=fish3_fem_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.13. Resid Comparison Trawl Fishery Length Comp Female.png"),plot=fish3_fem_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Domestic Trawl Fishery (Length Males) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt fish3 size comp male iter")]) 

# Get fishery length males residuals (OSA)
fish3_male_len_osa = get_osa_res(obs = sab_curr$olc.fish3.m, pred = sab_curr$elc.fish3.m, 
                                 iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get fishery length males residuals (Pearson)
fish3_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                         obs = sab_curr$olc.fish3.m, pred = sab_curr$elc.fish3.m)

# Plot OSA and Pearson Comparison
fish3_male_len_osa_plot = res_plot(data = fish3_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for Trawl Fishery Male Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

fish3_male_len_pearson_plot = res_plot(data = fish3_male_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for Trawl Fishery Male Length Compositions",
                        ymin = min(plot_lens))

fish3_male_len_comb_plot = plot_grid(fish3_male_len_osa_plot, fish3_male_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Trawl Fishery Length Comp Male.png"),plot=fish3_male_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Trawl Fishery Length Comp Male.png"),plot=fish3_male_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.12. Resid Comparison Trawl Fishery Length Comp Male.png"),plot=fish3_male_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Domestic Trawl Survey (Length Females) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt srv7 size comp female iter")]) 

# Get Survey length females residuals (OSA)
srv7_fem_len_osa = get_osa_res(obs = sab_curr$olc.srv7.f, pred = sab_curr$elc.srv7.f, 
                               iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length females residuals (Pearson)
srv7_fem_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                       obs = sab_curr$olc.srv7.f, pred = sab_curr$elc.srv7.f)

# Plot OSA and Pearson Comparison
srv7_fem_len_osa_plot = res_plot(data = srv7_fem_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for NOAA GOA Trawl Survey Female Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

srv7_fem_len_pearson_plot = res_plot(data = srv7_fem_len_pearson, 
                        comp_type = "Length", res_type = "Pearson Residuals for NOAA GOA Trawl Survey Female Length Compositions",
                        ymin = min(plot_lens))

srv7_fem_len_comb_plot = plot_grid(srv7_fem_len_osa_plot, srv7_fem_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Trawl Survey Length Comp Female.png"),plot=srv7_fem_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Trawl Survey Length Comp Female.png"),plot=srv7_fem_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.9. Resid Comparison Trawl Survey Length Comp Female.png"),plot=srv7_fem_len_comb_plot,width=12,height=10,dpi=300)


###########################################################################################
# Domestic Trawl Survey (Length Males) -----------------------------------------------------
###########################################################################################

# get effective sample size weights
iter_wt = parse_number(ctl_wts[str_detect(ctl_wts, "#wt srv7 size comp male iter")]) 

# Get Survey length males residuals (OSA)
srv7_male_len_osa = get_osa_res(obs = sab_curr$olc.srv7.m, pred = sab_curr$elc.srv7.m, 
                                iss = iss, iter_wt = iter_wt, index = osa_lens, drop_bin = 1)

# Get Survey length males residuals (Pearson)
srv7_male_len_pearson = get_pearson_res(iter_wt = iter_wt, iss = iss,
                                        obs = sab_curr$olc.srv7.m, pred = sab_curr$elc.srv7.m)

# Plot OSA and Pearson Comparison
srv7_male_len_osa_plot = res_plot(data = srv7_male_len_osa[[1]], 
                    comp_type = "Length", res_type = "OSA Residuals for NOAA GOA Trawl Survey Male Length Compositions (Remove 41cm Bin)",
                    ymin = min(plot_lens))

srv7_male_len_pearson_plot = res_plot(data = srv7_male_len_pearson, 
                        comp_type = "Length", 
                        res_type = "Pearson Residuals for NOAA GOA Trawl Survey Male Length Compositions",
                        ymin = min(plot_lens))

srv7_male_len_comb_plot = plot_grid(srv7_male_len_osa_plot, srv7_male_len_pearson_plot, ncol = 1)

ggsave(paste0(dir_extra , "//ALT_Fig Pearson Resid Trawl Survey Length Comp Male.png"),plot=srv7_male_len_pearson_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_extra , "//ALT_Fig OSA Resid Trawl Survey Length Comp Male.png"),plot=srv7_male_len_osa_plot,width=10,height=8,dpi=300)
ggsave(paste0(dir_results , "//Fig. 3C.8. Resid Comparison Trawl Survey Length Comp Male.png"),plot=srv7_male_len_comb_plot,width=12,height=10,dpi=300)


#############################################################################################################################################
# size-age transition plots
############################################################################################################################################

clr.pos = rgb(255, 100, 100, maxColorValue=255)    # color for positive residuals
clr.neg = rgb(200, 100, 200, maxColorValue=255)    # color for positive residuals
clr.neg2 = rgb(155, 100, 255, maxColorValue=255)    # color for positive residuals
#
png(file=paste0(dir_extra,"//ALT_Fig Age-length conversion matrix.png", sep=''),res=500,width=8,height=8.,units="in")
par(omi=c(0,0,0,0),mfrow=c(2,2),mgp=c(2.5,1,0),mai=c(0.8,0.7,.5,0.1))

z<-as.matrix(sab_curr$sizeage.f.block1)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.f)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.f.block1)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))
plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T)
title(main="Size-Age Transition Block 1, Female")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }
z<-as.matrix(sab_curr$sizeage.f.block2)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.f)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.f.block2)[[2]])
#y<-as.numeric(dimnames(z)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))
plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T, main="Size-Age Transition Block 2, Female")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }
#
z<-as.matrix(sab_curr$sizeage.m.block1)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.m)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.m.block1)[[2]])
#y<-as.numeric(dimnames(z)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))
plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T, main="Size-Age Transition Block 1, Male")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }
z<-as.matrix(sab_curr$sizeage.m.block2)*10
x1<-as.numeric(dimnames(sab_curr$olc.fish1.m)[[2]])
y1<-as.numeric(dimnames(sab_curr$sizeage.m.block2)[[2]])
#y<-as.numeric(dimnames(z)[[2]])
x2 <- rep(y1,ncol(z))
y2 <- sort(rep(x1,nrow(z)))


plot(x2,y2,xlab="Age (yrs)",ylab="Length (cm)",type="n", las=T, main="Size-Age Transition Block 2, Male")
for(p in 1:length(x2))
{ if(z[p]>0.5*max(z)) {  clr=clr.pos}
  else {if(z[p]>0.25*max(z)) {clr=clr.neg} else {clr=clr.neg2}}
  points(x2[p],y2[p],cex=z[p],col=clr,pch=16)
  points(x2[p],y2[p],cex=z[p],col=1,pch=1) }

dev.off()


#############################################################################################################################################
# lls abund old fish
############################################################################################################################################

par(omi=c(0,0.1,0,0.25),mgp=c(2.5,1,0))

lls_age_dist_fig <- lls_age_filter %>%
  dplyr::select(year,area,Age,obs) %>%
  dplyr::group_by(year,Age) %>%
  dplyr::summarize(Freq = sum(obs)) %>%
  dplyr::mutate(Bin = dplyr::case_when(Age < 12 ~ '2-11',
                                       Age >= 12 & Age < 22 ~ '12-21',
                                       Age >= 22 ~ '22+') ,
                Bin = factor(Bin, levels = c('2-11', '12-21', '22+'), ordered = TRUE)) %>%
  dplyr::group_by(year, Bin) %>%
  dplyr::summarize(Freq = sum(Freq, na.rm=T)) %>%
  ggplot(aes(x = year, y = Freq, color = Bin)) +
  geom_line(lwd= 1.5)+
  ggtitle("Trends in Abundance by Age Group on the Domestic Longline Survey")+
  labs(y = "Number of Otoliths Aged", x = "Year", color = 'Age Group')+
  scale_color_jco()


ggsave(paste0(dir_extra , "//ALT_Fig LLS freq of old fish.png"),width=8,height=5,dpi=300,plot=lls_age_dist_fig) 



#############################################################################################################################################
# MCMC Correlations

# NOTE: FOLLOWING MCMC PLOTS MAY NOT BE ACCURATELY ADJUSTED FOR THE MESSY MCMC FILE AND/OR RECENT CHANGES TO MCMC APPROACH
############################################################################################################################################

### Plot correlation of MCMC matrix
## panel.cor function for making loess and numbers on graph
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  text(0.5, 0.5, txt, cex = cex * r) 
  #text(.8, .8, Signif, cex=cex, col=2) # for adding stars for significance
}


b<-cbind(mcmc_mess$mcmc[,4],mcmc_mess$mcmc[,5],mcmc_mess$mcmc[,8],mcmc_mess$mcmc[,9],mcmc_mess$mcmc[,(2*years+14)],mcmc_mess$mcmc[,(3*years+30)],
         mcmc_mess$mcmc[,(3*years+4)],mcmc_mess$mcmc[,12])   
colnames(b)<-c("Dom.LL.q","Jap.LL.q","Dom.Fish.q","Jap.Fish.q","Ending SSB","ABC","2000 YC","F40%")

png(file=paste0(dir_extra,"//ALT_Fig MCMC Correlations.png", sep=''),res=500,width=8,height=8.,units="in")

pairs(b,lower.panel=panel.smooth, upper.panel=panel.cor)

dev.off()


#############################################################################################################################################
# SSB MCMC with Projections
############################################################################################################################################

source(paste0(dir_R,"//newswath.r",sep=''))
### MCMC projection
png(file=paste0(dir_extra,"//ALT_Fig SSB MCMC w Proj.png", sep=''),res=500,width=8,height=5.,units="in")
par(omi=c(0,0,0,0),mgp=c(2.5,1,0),mai=c(0.8,0.7,.5,0.1))

newswath(mcmc_mess,min(all_years),curr_yr)
dev.off()


#############################################################################################################################################
# SSB Projected Probability Below BRPs 
############################################################################################################################################

d<-mcmc_mess$mcmc[,(3*years+15):(3*years+29)]

belowb40<-d/B40
belowb35<-d/B35
belowmmst<-d/(0.5*B35)
pbelowb35<-seq(1,14)
pbelowb40<-seq(1,14)
pbelowmsst<-seq(1,14)
for(i in 1:14) {
  pbelowb40[i]<- length(which(belowb40[,i]<1))/length(belowb40[,i])
  pbelowb35[i]<- length(which(belowb35[,i]<1))/length(belowb40[,i])
  pbelowmsst[i]<- length(which(belowmmst[,i]<1))/length(belowb40[,i]) }

png(file=paste0(dir_extra,"//ALT_Fig MCMC prob below BRP.png", sep=''),res=500,width=8,height=8.,units="in")
par(omi=c(0,0,0,0),mgp=c(2.5,1,0),mai=c(0.8,0.7,.5,0.1))

plot(seq(curr_yr+1,curr_yr+14),pbelowmsst,ylim=c(0,1),col=4,xlab="Year",ylab="Probability",pch="")
lines(seq(curr_yr+1,curr_yr+14),pbelowmsst,lty=2,col=4,lwd=3)
points(seq(curr_yr+1,curr_yr+14),pbelowb35,col=3)
lines(seq(curr_yr+1,curr_yr+14),pbelowb35,col=3,lty=3,lwd=3)
points(seq(curr_yr+1,curr_yr+14),pbelowb40,col=9)
lines(seq(curr_yr+1,curr_yr+14),pbelowb40,col=9,lty=1,lwd=3)
legend(curr_yr+2,0.48,c("P<B40%","P<B35%","P<B17.5%"),lty=c(1,3,2),col=c(1,3,4),cex=0.9,lwd=c(3,3,3))

dev.off()


###################################################################################################################################################
# Write PDF of all plots
###################################################################################################################################################

pdf(paste0(dir_extra,"//Model_",model_name_curr,"_Outputs(",format(Sys.time(), "%d-%b-%Y %H.%M"),").pdf", sep=''))

  print(catch_gear_plot)
  print(catch_area_graph)
  print(survey_plot)
  print(lls_area_fig)
  print(lls_dep_fig)
  print(fish_whale_dep_plot)
  print(like_plot)
  print(index_plots)
  print(ac_srv2_mean_plot)
  print(ac_srv2_plot)
  print(ac_srv1_mean_plot)
  print(ac_srv1_plot)
  print(ac_fish1_mean_plot)
  print(ac_fish1_plot)
  print(lc_srv2_mean_plot)
  print(lc_srv2_male_plot)
  print(lc_srv2_female_plot)
  print(lc_srv1_mean_plot)
  print(lc_srv1_male_plot)
  print(lc_srv1_female_plot)
  print(lc_srv7_mean_plot)
  print(lc_srv7_male_plot)
  print(lc_srv7_female_plot)
  print(lc_fish1_mean_plot)
  print(lc_fish1_male_plot)
  print(lc_fish1_female_plot)
  print(lc_fish3_mean_plot)
  print(lc_fish3_male_plot)
  print(lc_fish3_female_plot)
  print(bio_ssb_plot)
  print(n_at_age_F_plot)
  print(n_at_age_M_plot)
  print(N_age_rec_plot)
  print(recruit_plot)
  print(rec_ssb_plot)
  print(recruit_plot_comp)
  print(sel_plot)
  print(f_plot)
  print(phase_plot)
  print(plot_proj_ssb)
  print(ssb_cont_plot)
  print(rec_ssb_catch_plot)
  print(fish1_age_resid_comb_plot)
  print(srv1_age_comb_plot)
  print(srv2_age_comb_plot)
  print(srv1_fem_len_comb_plot)
  print(srv1_male_len_comb_plot)
  print(srv2_fem_len_comb_plot)
  print(srv2_male_len_comb_plot)
  print(fish1_fem_len_comb_plot)
  print(fish1_male_len_comb_plot)
  print(fish3_fem_len_comb_plot)
  print(fish3_male_len_comb_plot)
  print(srv7_fem_len_comb_plot)
  print(srv7_male_len_comb_plot)
  print(lls_age_dist_fig)

  dev.off()

  
######################################################################################################################################################################  
  save.image(file=paste0(dir_master,"//Sablefish_Results_",year,".RData",sep=''))        # Save all the data frames in case want to reload data without doing full data pulls
######################################################################################################################################################################
