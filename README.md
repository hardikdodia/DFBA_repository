# DFBA_repository
Install the toolbox and packages given below before running the program.
1. CPLEX Optimization Studio by IBM
2. ARESLab Toolbox
3. COBRA Toolbox

In the MATLAB, Open folder titled "DFBA_Inputs"
Open "FBA_Model_Inputs.m"
Add path to the ARESLab toolbox.

Set few parameters when prompted. The recommended parameters for our dataset are

Please enter the time-step for MARS: 0.01
Please enter the time-step you want to extract from MARS: 0.01
Please enter the final number of data points desired: 50

Please select the directory to save PreliminaryDataFBA.mat
Here, give the path to "DFBA_Maincode" folder

Would you like to read the the Excel Data file again? (1/0): 1

Select the excel data file...
Here, give the path to "DFBA Model Inputs" file

Is Experimental Dissolved Oxygen data available? (1/0): 1

Is True Absolute Metabolites Mf data available? (1/0) : 1

Is Converted Absolute Metabolites Mf data available? (1/0): 1
This imples the Mf concentration calculated from standard curve

This program will also create presentations of plots of all the metabolites concentrations and fluxes.

After the programs runs, move to the folder, "DFBA_Maincode". You should be able to see a file named "PreliminaryDataFBA.mat" generated.
Open "Main_FBA_script.m"
If you plan to run flux balance analysis (FBA), set FVAcheckwt = 0 and run the program.
If you plan to run flux variability analysis (FVA), set FVAcheckwt = 1 and run the program.

The output flux results will be stored in "DFBA_results.xlsx" file.

After running this program, run "FBA_Model_Results.m" to obtain the plots of all the output fluxes.
