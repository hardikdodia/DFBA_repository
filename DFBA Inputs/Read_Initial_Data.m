%% Read Initial Data from Excel File %%

%% Browse PC to select Excel File Destination Folder %%
excelfilename = uigetfile('.xlsx');

%% Specify if Mf & DO data is available and read accordingly %%
Expt_Data_DO_Available = input("Is Experimental Dissolved Oxygen data available? (1/0) \n");

% Specify if True Absolute Metabolites Mf data is available and read accordingly %
Mf_Data_TAbs_Available = input("Is True Absolute Metabolites Mf data available? (1/0) \n");

% Specify if Converted Absolute Metabolites Mf data is available and read accordingly %
Mf_Data_CAbs_Available = input("Is Converted Absolute Metabolites Mf data available? (1/0) \n");

disp("Reading Excel File");

if Expt_Data_DO_Available == 1
    disp("Reading DO Data");
    Expt_Data_Culture_DO_1 = readtable(excelfilename,'Sheet','Expt Data Culture DO 1','VariableNamingRule','preserve');
    Expt_Data_Culture_DO_2 = readtable(excelfilename,'Sheet','Expt Data Culture DO 2','VariableNamingRule','preserve');
end

if Mf_Data_TAbs_Available == 1
    disp("Reading Mf Data for True Absolute Metabolites");
    Mf_Data_TAbs = readtable(excelfilename,'Sheet','MfKd TAbs','VariableNamingRule','preserve');
end

if Mf_Data_CAbs_Available == 1
    disp("Reading Mf Data for Converted Absolute Metabolites");
    Mf_Data_CAbs = readtable(excelfilename,'Sheet','MfKd CAbs','VariableNamingRule','preserve');
end

%% Metabolite Experimental Data %%
% Standard Curve Data %
disp("Reading Standard Curve Data");
Std_Curve_Data = readtable(excelfilename,'Sheet','Std Curve','VariableNamingRule','preserve');

% Experimental Data - Culture - True Absolute Metabolites %
disp("Reading Experimental Data - Culture - True Absolute Metabolites");
Expt_Data_Culture_TAbs_1 = readtable(excelfilename,'Sheet','Expt Data Culture TAbs 1','VariableNamingRule','preserve');
Expt_Data_Culture_TAbs_2 = readtable(excelfilename,'Sheet','Expt Data Culture TAbs 2','VariableNamingRule','preserve');

% Experimental Data - Culture - Converted Absolute Metabolites %
disp("Reading Experimental Data - Culture - Converted Absolute Metabolites");
Expt_Data_Culture_CAbs_1 = readtable(excelfilename,'Sheet','Expt Data Culture CAbs 1','VariableNamingRule','preserve');
Expt_Data_Culture_CAbs_2 = readtable(excelfilename,'Sheet','Expt Data Culture CAbs 2','VariableNamingRule','preserve');

% Experimental Data - Culture - Relative Metabolites %
disp("Reading Experimental Data - Culture - Relative Metabolites");
Expt_Data_Culture_Rel_1 = readtable(excelfilename,'Sheet','Expt Data Culture Rel 1','VariableNamingRule','preserve');
Expt_Data_Culture_Rel_2 = readtable(excelfilename,'Sheet','Expt Data Culture Rel 2','VariableNamingRule','preserve');

% Experimental Data - NCC - True Absolute Metabolites %
disp("Reading Experimental Data - NCC - True Absolute Metabolites ");
Expt_Data_NCC_TAbs_1 = readtable(excelfilename,'Sheet','Expt Data NCC TAbs 1','VariableNamingRule','preserve');
Expt_Data_NCC_TAbs_2 = readtable(excelfilename,'Sheet','Expt Data NCC TAbs 2','VariableNamingRule','preserve');

% Experimental Data - NCC - Converted Absolute Metabolites %
disp("Reading Experimental Data - NCC - Converted Absolute Metabolites");
Expt_Data_NCC_CAbs_1 = readtable(excelfilename,'Sheet','Expt Data NCC CAbs 1','VariableNamingRule','preserve');
Expt_Data_NCC_CAbs_2 = readtable(excelfilename,'Sheet','Expt Data NCC CAbs 2','VariableNamingRule','preserve');

% Experimental Data - NCC - Relative Metabolites %
disp("Reading Experimental Data - NCC - Relative Metabolites");
Expt_Data_NCC_Rel_1 = readtable(excelfilename,'Sheet','Expt Data NCC Rel 1','VariableNamingRule','preserve');
Expt_Data_NCC_Rel_2 = readtable(excelfilename,'Sheet','Expt Data NCC Rel 2','VariableNamingRule','preserve');

% Ind_Flux & Ind_Fluxr %
disp("Reading IndFlux, IndFluxr Data");
IndFlux_Data = readtable(excelfilename,'Sheet','Ind_Flux','VariableNamingRule','preserve');
IndFluxr_Data = readtable(excelfilename,'Sheet','Ind_Fluxr','VariableNamingRule','preserve');

save('Initial_Data.mat');

disp("Reading Excel File Complete");


