# CCHF
Cryosphere hydrology modeling package written in Matlab. This package requires that the 'shared_subroutines_Matlab' repository also be in the search path.

I’ve packaged all the necessary functions to implement the CCHF for Wolverine glacier, which is a benchmark glacier maintained by the USGS (there are glacier stake and streamgage measurements from the 1950’s through present at this site).

The first step is to run the main script (“CCHF_main_*.m”), which is located at “./CCHF_demo_*/CCHF”. Once you run the main script, you’ll then be prompted to locate inputs, in this order:
1.	Select the digital elevation model (DEM), which is located at “./CCHF_demo_*/Wolverine_inputs/wolverine_dem.asc”
2.	Select the flow direction grid (FDR), which is located at “./CCHF_demo_*/Wolverine_inputs/wolverine_fdr_manual.asc”
3.	Select the glacier presence grid (i.e. Boolean values), which is located at “./CCHF_demo_*/Wolverine_inputs/wolverine_rgi5.asc”
4.	Select the folder with the precipitation time-series, which is located at “./CCHF_demo_*/Wolverine_inputs/CFSR_daily/pre”
5.	Select the folder with the mean temperature time-series, which is located at “./CCHF_demo_*/Wolverine_inputs/CFSR_daily/tasmean”
6.	Select the folder with the max temperature time-series, which is located at “./CCHF_demo_*/Wolverine_inputs/CFSR_daily/tasmax”
7.	Select the folder with the minimum temperature time-series, which is located at “./CCHF_demo_*/Wolverine_inputs/CFSR_daily/tasmin”
8.	Select the observation data (used for calibration and validation), which is located at “./CCHF_demo_*/Wolverine_inputs/SETI_gage-flow-stake-sca.txt”
After making the above selections, the program will run on its own. Several actions occur in the main script, including calculating watershed geometric relationships, creating structure arrays of a specific format for handling the necessary inputs, intermediaries, and outputs. Lines 587 – 605 implement the model run(s) (i.e. calibration, validation, or default). I have the model setup to perform a calibration run, using a hybrid optimization routine I’ve developed (uses Monte Carlo simulations, particle swarm optimization, and linear sensitivity). Also, note that during calibration, CCHF runs in parallel by default.


A “calibration” parameter set run for Wolverine takes about a 24 hours using a ten year period with a daily time step because calibration runs the model several thousand times. A “default” parameter set run takes about ten minutes. A “validation” run also takes about ten minutes, but can only be used once calibration has completed. I’ve included a parameter set created during calibration at the path “./CCHF_demo_*/Wolverine_inputs/gulkana_LST_parameters.txt”, which can be used in a validation run. You will be prompted to locate this file approximately 30 seconds after the prompts to locate the inputs listed above.


Several plots are generated after any model run. All of these figures and several other files are written to a folder unique to the run within “./CCHF_demo_*/Wolverine_inputs”. Plots include time-series of requested output at points of interest (default is domain-averaged values and observation points). During calibration, plots also include optimized parameter values.



 
Required Modules (must have one of each):
*	'energy' (calculates heatflux):
  *	'simple'
  *	'Pelli' (Pellicciotti version of ETI)
  *	'Hock' (Hock version of ETI)
  *	'SETI'
*	'mass' (determines link between energy and melt):
  *	'tstep'
  *	'hstep'
  *	'tinertia'
  *	'hinertia'
  *	'enbal'
*	'trans' (atm transmissivity) (used in shortwave radiation):
  *	'Coops'
  *	'DeWalle'
*	'albedo' (not used if 'cryo-simple') (used to determine absorption of shortwave radiation):
  *	'Brock'
  *	'Pelli'
*	'PET' (Potential evapotranspiration; actual ET limited by water availability):
  *	'Hammon'
*	'runoff':
  *	'bucket' (groundwater bucket model)
  *	'direct' (surface only)
*	'time' (refers to travel time between cells):
  *	'Johnstone'
  *	'Liston' (allows variable travel time depending on landcover of cell)
*	'flow':
  *	'Muskingum'
  *	'lumped'
  *	'Liston'
Optional modules (not needed to run)
*	%'glacier' (if term included, model will allow glacier sliding and calculate changes in ice thickness)
*	%'holding' (if term included, snowpack has a liquid water holding capacity of 3%)

Naming
*	Cryosphere (‘sCryo’)
  *	‘heat’ – net heat flux
*	Land-layer (‘sLand’)

*	Atmosphere (‘sAtm’)
  *	‘tas’
  *	‘tasmax’
  *	‘tasmin’
  *	‘pr’
  *	‘rain’
  *	‘prsn’
=======
•	'energy' (calculates heatflux): 
o	'simple'
o	'Pelli' (Pellicciotti version of ETI)
o	'Hock' (Hock version of ETI)
o	'SETI'
•	'mass' (determines link between energy and melt): 
o	'tstep'
o	'hstep'
o	'tinertia'
o	'hinertia'
o	'enbal'
•	'trans' (atm transmissivity) (used in shortwave radiation): 
o	'Coops'
o	'DeWalle'
•	'albedo' (not used if 'cryo-simple') (used to determine absorption of shortwave radiation): 
o	'Brock'
o	'Pelli'
•	'PET' (Potential evapotranspiration; actual ET limited by water availability): 
o	'Hammon'
•	'runoff': 
o	'bucket' (groundwater bucket model) 
o	'direct' (surface only)
•	'time' (refers to travel time between cells): 
o	'Johnstone'
o	'Liston' (allows variable travel time depending on landcover of cell)
•	'flow': 
o	'Muskingum'
o	'lumped'
o	'Liston'
Optional modules (not needed to run)
•	%'glacier' (if term included, model will allow glacier sliding and calculate changes in ice thickness)
•	%'holding' (if term included, snowpack has a liquid water holding capacity of 3%)

Naming
•	Cryosphere (‘sCryo’)
o	‘heat’ – net heat flux
•	Land-layer (‘sLand’)
o	
•	Atmosphere (‘sAtm’)
o	‘tas’
o	‘tasmax’
o	‘tasmin’
o	‘pr’
o	‘rain’
o	‘prsn’
