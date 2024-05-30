# SBCEX-2022-merchant-ship
This repository contains all code and information required for a deep sediment heterogeneity analysis of the New England Mud Patch using very low-frequency features from merchant ships. The data from this analysis come from the Seabed Characterization Experiment 2022 from three vertical lines arrays from Marine Physical Laboratory (MPL) Scripps Institute of Oceanography (SIO) and Applied Research Laboratories, The University of Texas at Austin. 

Below is a explanation of what is contained in each folder of the repository. 

Figures Codes
    -	JASA_PLOTS.m ------ creates the figures from this analysis that appear in the JASA “Deep sediment heterogeneity inferred using very low-frequency features from merchant ships”  
        This folder also contains the powerpoint and .jpg figures from the JASA. 

Inversion Codes
    -	SBCEX_featurebased_inversion.m ------ uses a maximum entropy inversion using the beta -> infinity singularities of very low-frequency merchant ship spectrograms to estimate the sound speed and thickness of         the deeper layers of the New England Mud Patch. User must specify the ship, related parameters, and the singularities found from the measured data. Output is a modeled spectrogram using optimal seabed              parameters that returned the lowest cost function in the inversion.      
        It also returns a distribution file containing the cost for each iteration of the inversion. For this code to work, the MATLAB ORCA repository is required. 
        This code requires the following codes: s00spectrum.m, SBCEX22_SSP, and wateratten.m
    -	marginal_probability_distributions.m ------ This code used the distribution file output of the maximum entropy inversion to plot the marginal probability distribution of the four inferred deep seabed               parameters.
        Requires: stats.m

Mapping Ship Tracks 
    -	plot_ship_track.m ------ Plots a map of the NEMP with the 24 interested ships. The folder Ship coordinates spreadsheets contains the coordinates of the shipsbased on the AIS data. Other files in this               folder are requires for this script to run. 

Merchant Ship Spreadsheets
    -	The three .xlsx files with information about the 24 ships. SHIPS_SBCEX2022_PROTEUS and VLA1_VLA2.xlsx contain the list of all ships that passes within 15 km of the arrays with the distance, SOG, ship ID,           and time stamps. SBCEX_ALL_INTERESTED_SHIPS_FINAL.xlsx contains the inferred results for the deep seabed layers for all interested ships. 

Processing AIS Data
    -	It contains the two raw .csv files downloaded from AIS database which can be obtained from https://marinecadastre.gov/accessais/. The two scripts use the latitude and longitude coordinates of the three             VLAs to create the .xlsx files with the list of ships that had passed within 15 km of the arrays. 

Processing array file types
    -	This folder contains codes to read the data from the arrays and convert the file types to one that are useable with the spectrogram plotting codes. 

SOO .mat files 
    -	The .mat files contain the time, frequency, and power spectral levels for each of the ships. These files are created with the spectrogram_plotting.m code. 

SOO dist files
    -	The marginal probability distribution plots for the four inferred parameters as well as the .mat dist files for each of the ships after the maximum entropy inversion and including information about the             ship and the optimal and peaked parameters.

Spectrogram Plotting Codes 
    -	Spectrogram_plotting.m ------ Plots the spectrograms from any number of channels of a particular ship and save the SOO .mat files. 
    -	Plot_measured_and_modeled_soo_specs.m ------ Uses the information from SBCEX_ALL_INTERESTED_SHIPS_FINAL.xlsx to plot both the measured and modeled spectrograms for each ship from only the specified channel         in the spreadsheet.

Spectrograms 
    -	The measured and modeled spectrogram for all ships. 

SSP 
    -	The code to make the SSP 2022 .mat file required for the orca call in the inversion code. It also contains the raw CTD files. 

The order of file use from this repository: 
1.	Processing array file types.
2.	Processing AIS data to find interested ships.
3.	Spectrogram plotting codes (spectrogram_plotting.m) to find the ships with clear SNR. The SOO .mat files are saved to a folder.
4.	Analyze all ship spectrograms and record the beta=infinity frequencies for each ship
5.	Maximum Entropy Inversion with 4,000 Monte Carlo iterations (SBCEX_featurebased_inversion.m). The initial dist file for each ship is saved to a folder. 
6.	Marginal probability distributions (marginal_probability_distributions.m). The dist files for each ship are saved again with the added information from the marginal probability distribution. 
7.	Mapping ship tracks. Create the .csv files with the latitude and longitude coordinates for each of the ships and then plot the tracks on a map of the NEMP.
8.	JASA_PLOTS.m to provide an analysis of the deep layer sound speeds. 



