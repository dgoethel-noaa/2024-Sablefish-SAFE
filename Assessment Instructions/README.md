Steps for running the sablefish assessment and producing all needed outputs (Figure and Tables) for SAFE.

NOTE: Each step is in a clearly marked folder. Within most folders is a README describing how to use the associated R code to run the given analysis.

NOTE: The spreadsheet '2023 SAFE Assignments and Sources for Tables_Figs.xlsx' provides additional guidance on developing SAFE doc (where to get Tables/Figs and contacts for certain sections not done by lead author).

1. Copy folders from previous year, including '_Final Model' and '_Final Data'.
2. Pull the data and get the final .dat file:
	a. In the '_Final Data' folder, open the 'Sable_Data_Pull_Final.R'.
	b. Update the inputs to reflect current year, ABCs, etc. and make sure all necessary .csvs are in the correct folders.
	c. Paste last year's tem.ctl, tem.dat, and tem.pin into the working directory.
	d. Run this code to pull the raw data from AKFIN and create the final tem_2023_na_wh.dat used for the admb tem.tpl model; it will also update the ctl and pin files to be ready to run with the additional model year
	e. Perform comparitive analysis with last year's .dat file to ensure no major changes and that all new data was actually pulled (i.e., is up on AKFIN).
3. Setup the initial model run and find the final Base Model:
	a. Copy the final .dat, .ctl, and .pin files into the model building directory for the continuity model along with the most recent tem.tpl, tem.exe, sable-r-report.cxx, mhp-s-funcs.cpp   (latter 2 are additional libraries for functions, etc.)
	b. Run model and debug focusing on making sure right number of estimated parameters for the new data year, etc. (e.g., does .dat, .pin reflect added year of data and associated additional recruit and F dev parameters).
	c. Perform any model updates to perform as part of the model bridging, or skip ahead if update assessment and no changes.
	d. Perform Francis reweighting to update the compositional data weights to account for new data.
		i. This code should be explored a bit more in depth as there is a tendency to upweight fishery lengths over ages, which is odd.
			It might be a factor of fitting both ages and lengths, which should also be explored--may be better to drop lengths in years with ages.
	e. Once new weights are defined and put in ctl, perform a jitter to find the model fit with the lowest negative log-likelihood (best fit to data).
		i. This code might be updated to base jitter comparison only on data weights, so not just fitting to penalty/prior funtion nLL.
	f. If a model with a lower nLL is found then use this as the base model.
		i. May want to iterate through using .par as a new .pin a few times just to hone in on the best nLL and gradient. Sometimes can get better gradient by updating .pin, etc. Could also redo Francis reweight, but unlikely to change much.
	g. This is the final base model for the current SAFE.
	h. Create the final model bridging plot (Fig. 3.9) by pasting in .par and .rep files from various runs then by running 'model_comp_RUN.r' to get output graphics.
		i. Check that time series have not greatly changed and ABCs appear appropriate. 
		NOTE: the F40 values are no longer estimated, so the internal ABC calculation will not be correct for each model run. The F40 phase can be set to positive for each run to get the ABC, then turned off again and rerun for the final model. This will ensure the .par uses the correct F40 and output ABC even though not estimated in final run.
4. Run MCMC:
	a. Copy the final base model files into the MCMC folder.
	b. Open the tem.tpl and use ADMB 'run with args' (in AD Studio) then type -mcmc 5000000 -mcsave 1000 (or just enter run time arguments if in command line)
	c. Once MCMCs are run, open tem.tpl and run with args -mceval (there is a line in the tpl telling ADMB which parameters to store/track for MCMC).
	d. An output file called 'evalout.prj' will be created, rename this 'evalout.dat'.
	e. Open 'evalout.dat' in a text editor and insert a 3 line header to the file, which reads:
			5000 #(number of thinned chain)
			64   #(number of years) .... NOTE UPDATE THIS EVERY YEAR FOR NEW NUMBER OF YEARS (64 is for terminal year 2023)
			30   #(number of ages)
	f. In the MCMC folder there should be an MCMC executable to parse the outputs, called 'evalout.exe'; run this in the command line (or just open .exe and it will run)
	g. An output file 'evalout.sdat' will be created which is used to create MCMC plots and tables; copy this to the plotting directly or main base model folder.
    NOTE: C. Monnahan has updated MCMC files that should be integrated in future, probably much quicker and better. TMB probably has a built in MCMC, but need to check if switch to TMB.

5. Retrospective Analysis with Reweighting
    This script works much like the Francis reweight script, but also removes data automatically for the retro. The reweight adds time to retro run, but makes it more realistic to avoid over inflating apparent retro if just used current weights with old data. 
    NOTE: the .ctl file must have a "#-" marker after any values that need to change due to a change in terminal model year
	(mainly, end yr for rec devs est, ph_IFQ_block2, ph_LL_block2, yr_sel_chg_fish, yr_sel_chg_srv3, ph_q_LL_srv_rec, ph_q_IFQ_rec)
    NOTE: Same is true for .dat file with begin marker ("#-") and end marker ("#!") for each section (this is already coded in .dat file creator)
	a. Make sure 'Run iterative ESS_retro.R' is in the Code folder.
	b. Copy base model files into Model folder, and copy .pin, .ctl, and .dat files into the Input Files folder; SEE NOTE above about structure of control and data files.
	c. Open the 'Sab_Retro.R' file and update inputs for the current year and number of retro peels want to do.
	d. Run 'Sab_Retro.R', this run the retro analysis, compute Mohn's rho, and produce figures for SAFE.

6. 'Historical' Retrospectives
	a. 'All Model'--this produces a retro analysis showing the SSB time series and two year projections from previous assessments used for management advice.
		i. Copy the sable.rep files for previous assessments into the folder and rename as sable_YEAR.rep.
		ii. Copy the final base model report file and name as above.
		iii. Run the 'All Model Historical Retro (Fig 3.47).r' script.
		iv. This should produce the SAFE figure 3.47 and associated Mohn's rho for the 2 year projections.
	b. 'Current Model'--this produces an alternate retro analysis that applies the current year model to the actual data available during previous years (as opposed to just peeling off a year of data as in a normal retro).
	     NOTE: this is mostly done by hand (see associated README), and should/could be automated in the future.
		i. Create YEAR folder for each peel want to complete.
		ii. Copy the current year base model files into the Model folder and associated .ctl and .pin files into the Input Files folder (if doing reweighting).
		iii. Copy the .dat file from terminal year of the given retro run into the Model folder and the Input Files folder (if doing reweighting).
		iv. Update the .ctl, .pin, and .dat files to represent the given end year (i.e., adjust the .ctl to reflect the .dat file name, number of years of recruit devs, etc...in .pin make sure dev pars represent the number of years of estimated parameters).
		v. Run the 'Run iterative ESS.R' script in each YEAR folder to get the final model run for each peel (again the reweighting isn't strictly necessary, but helps produce more realistic results).
		vi. Copy the final sable.rep file into main folder and rename as sable_YEAR.rep.
		vii. Run the 'Current Model Historical Retro (Fig 3.48).r' script to get Mohn's rho for 2 year projections and final figure for SAFE doc.

7. Profile likelihoods
	a. Copy final base model files into the Model folder.
	b. Open 'profile_likelihood.r' script and update inputs, if needed (don't really need to change from year to year, except to fine tune the dimensions of the search).
	c. Run 'profile_likelihood.r' script, this produces a figure for the SAFE doc.
	d. If profile figure is a bit off, may need to adjust bounds of the runs or step size; generally these plots are a bit weird for sablefish partly because the response surface is not very smooth (probably need to improve handling of comps).

8. Index Sensitivity
    NOTE: this is done by brute force and each .ctl file needs to be updated manually (to turn data weight to 0 and any associated estimated 
          parameters to a negative phase, ie not estimated).
    NOTE: the tpl and exe files in each folder of this analysis are model 23.5 updated to turn off specific parameters to that model run (i.e., the selectivity parameters associated with a given survey). Thus, don't need to 
          worry about the selectivity parameters that are hardcoded (need to update ctl file to let all parameters be switched on and off from it...unless you update model and use a new tpl then need to recode these runs).
             *********** The specific steps for each index run are provided below. Mostly, don't need to worry about selectivity phases for each run, because these are hardcoded now for the given folder. If model changes and need to update tpl for each folder, then will need to rehardcode the selectivity phases for each run in this analysis.
	a. Set up folders for each model run, BUT DO NOT COPY IN A NEW TEM.tpl (there is some hard coding for each run that will be lost if tpl is replaced).
	b. Turn off estimation and set data weight to 0 in .ctl for given index being removed.
		i. LL survey:
			ph_q_srvy1 (and ph_q_LL_srv_rec if currently pos) set to -1
			Leave the ph_surv_sel and ph_surv_sel_delt as positive (the tpl in the LL survey folder has the phases for LL survey sel hardcoded to neg for this analysis)
			Set #wt DOM LL Srvy RPN #wt surv1 age comp iter #wt surv1 size comp male iter #wt surv1 size comp female iter all to 0
			Hard code phases for just the LL survey a50 parameters to negative (leave the delta par positive because this is shared across all fleets, so need to use it for JAP LL survey).
			Need to turn q_srvy1 from sd_report_number to number (and remove any report sections that try to report out sd value...these will all cause indefinite covariance matrix)
		ii. Trawl survey:
			ph_q_srvy7 set to -1
			Leave the ph_surv_sel and ph_surv_sel_delt as positive (the tpl in Trawl LL survey folder has the phases for Trawlsurvey sel hardcoded to neg for this analysis)
			Set #wt GOA Trawl bio #wt srv7 size comp #wt srv7 size comp male iter #wt srv7 size comp female iter all to 0
			Hard code phases for just the Trawl survey a50 parameters to negative (leave the delta par positive because this is shared across all fleets, so need to use it for LLsurvey and JAP LL survey).
			Need to turn q_srvy7 from sd_report_number to number (and remove any report sections that try to report out sd value...these will all cause indefinite covariance matrix)
		iii. Fishery CPUE:	
			#ph surv8 q  #ph_q_IFQ_rec  #ph surv5 q all set to 0 (srvy5 no longer used in 23.5 so should already be set to 0)
			#wt surv5 DOM CPUE RPW set to 0
			No changes to selectivity parameters needed since no separate selectivity fleet for the CPUE
			Need to turn q_srvy8 and q_LL_fish_recent from sd_report_number to number (and remove any report sections that try to report out sd value...these will all cause indefinite covariance matrix)
	c. Run each of the new models.			 
	d. Copy results to '_Model Comp' folder.
	e. Once all runs are complete and results copied to model comp, run model comp, and compare graphically.
	
9. Stepwise Data Addition
    NOTE: this is done by brute force and each .dat file needs to be updated manually (i.e., adding one data source at a time starting from previous year .dat file).
    NOTE: when update from previous year .dat file, make sure TO ADD A YEAR FOR SEX RATIO and update ENDYR input, then add new data
    NOTE FOR 2023: all model runs start with standardized CPUE (as implemented for 2023) and not nominal CPUE (used in 2022 SAFE)
	a. Set up folders for each model run.
		Runs are (note survey index added with a set of comp data only, not independently; catch is added in first step and included in all subsequent models; otherwise only the noted data is updated):
			"Add Catch",
			"Add Catch+Fixed Gear Fish Age",
			"Add Catch+Fixed Gear Fish Length",
			"Add Catch+Fixed Gear Fish Age+Length",
			"Add Catch+Fixed Gear Fish Age+Length+CPUE",
                	"Add Catch+Trawl Gear Fish Length", 
			"Add Catch+LL Srvy Age", 
			"Add Catch+LL Srvy Length", 
			"Add Catch+LL Srvy Age+Length"
                	"Add Catch+Trawl Srvy Length"
	b. Add given data set (replace all years, not just terminal year) to previous year .dat file and run the model.
	c. Copy results to '_Model Comp' folder.
	d. Once all runs are complete and results copied to model comp, run model comp, and compare graphically.

10. Sensitivity Runs
	a. Setup folders and add base model files.
	b. Make changes to .tpl, .ctl., or .dat file and run model.
	c. Add final model files to main folder and update and Run 'model_comp_RUN.r'.

11. Final Projections
	a. See the 'PROJ Steps README.docx' for complete instructions, but this essentially involves updating the PROJ model inputs for the current year and catch values.
	b. Once run copy results into the 'Table 3.11 Sable_Projections_2023.xlsx' spreadsheet as instructed.
	c. Once complete, the spreadsheet will contain the final SAFE Summary Table (except whale depredation values) and Table 3.11 projection outputs.
	d. The final ABCs and OFLs will then be input into the '_Exec summary tables_All proportions_changable ABC_2023' spreadsheet to do whale depredation calcs and apportionment.

12. Final Apportionment
	a. Open the '_Exec summary tables_All proportions_changable ABC_2023' spreadsheet and associated README file for full instructions.
	b. Copy over the 'apportionment.csv' and 'total_depredation_area.csv' produced by the 'Sable_Data_Pull_Final.r' code.
	c. Populate the cells in red on the 'Apportionment' tab with appropriate values from the projection spreadsheet, the .CSVs noted above, and last year catch and harvest recs (ABC, OFL, etc....these should already be here from prev year SAFE).
		i. Put in new ABC and OFL from final projections.
	 	ii. Copy in the new 5 year average survey apportionments from 'apportionment.csv' output from the data pull code.
		iii. Copy last year's final ABCw and OFLw from the associated harvest spec sheets (downloaded from AKFIN Answers).
		iv. Copy in other area apportionment rates from previous year that want to use for comparisons (5 yr ave., SSC rec., etc.)
		v. Copy whale depredation by area for last 3 years (excluding partial year estimates for current year) from 'total_depredation_area.csv'. 
		vi. Copy in current year catch (partial year) and TAC.
		vii. Copy in RPW by area for terminal year.
		viii. Copy in biomass for terminal year and projected biomass from projections spreadsheet (for harvest rate calcs).
	d. Increment each year reference in all tabs by +1 to account for new year. 
	e. Tables should automatically populate.
	f. For 2023 and future, will use full 5 year survey average proportions (no more stair step).
	g. All of the executive summary tables should now be created, but may need some cleaning (or can just copy/paste values into tables already in SAFE).

13. Final Tables and Figures
	a. Open the 'Sablefish Plots SAFE_Final.R' in the 'R code' folder of the base model.
	b. Copy the final output from the data pull code 'Sablefish_Data_YEAR.RData' into the 'R code' folder.
	c. Update the inputs to the r code, including model names, terminal year, last year ABCs, etc.
	d. Make sure all necessary files are in the main folder, including MCMC results and previous year .rdat results file (see below for full list of files used in this script, those without YEAR are the results files from current base model assessment)
	e. Run the script.
	f. All of the figures for the SAFE (excluding the secondary analysis figs) and the inputs for each table will be produced by this script.
		#########################################################################
		###### Required Data Files #############################################
		########################################################################
			# tem.rdat                     // current year SA report file output as rdat file
			# tem_YEAR.rdat                // previous year SA report file output as rdat file, "2022" should be replaced with year of prev. SA
			# sable.rep                    // ADMB output report file of current year SA
			# tem.par                      // ADMB output parameter file of current year SA
			# tem.std                      // ADMB output standard deviation file from current year SA
			# tem.ctl                      // ADMB control file of current year SA, for residual plots
			# tem_YEAR.ctl                 // ADMB control file of previous year SA, for residual plots
			# evalout.sdat                 // MCMC ADMB output from current year SA, need to run Bayesian version of assessment to obtain
			# Sablefish_Data_2023.RData    // The Rdata file containing all of the assessment inputs and data pulls for the current year

14. SAFE Writing 
	a. Copy final figures into Figure folder; need to copy the figures from alternate analyses (retros, sensitivity, profiles, etc.) from those folders and add appropriate figure numbers.
	b. Copy the raw table inputs into the Table folder.
	c. Each numbered table spreadsheet should have the formatted table, the outputs from the R script can then be pasted directly into these tables.
	d. In the final SAFE doc, paste the new formatted table or just the values into the table already in the SAFE; formatting can be a beast, so latter may be easier.
	e. Replace figures with updated figures; only figure not produced by this code are those from the 'other sablefish assessments' in Figure 1 (need to email individual authors to get these).
	f. Start writing, mostly just replacing values and adding any new trends (though usually not much changes in a year!); focus is probably on the model bridging if any major updates.
	g. Make sure to communicate early if expecting others to contribute (e.g., Appendices, getting updates from other regions that manage sablefish, and the Ecosystem Risk Table section from ESP and ESR authors).

Post Assessment:
	1. Pull ages for Jon Short (Survey usually in Nov., Fishery in Feb...once AKFIN updated with final data)
		i. Just run the 'fishagesampleYR.r' code and provide the resulting spreadsheet to Jon....it is mostly automated.
		ii. Put in an associated ageing request in the AGP portal (https://afsc-apps-internal.nmfs.local/al/agps2/agps2_home.php) and let Jon know.
	2. Provide the SARA files and SIS to Kalei Shotwell.
		i. These are auto-produced by the assessment code tem.tpl, but may need some slight tweaking based on guidance in the email request sent out in Nov.
	3. Make .ppt for JPT and director's briefing (can use templates from previous year).
	4. Update github with final models, etc.
	5. Have a giant beer.