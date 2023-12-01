%% FBA Model Inputs %

% This script must be run first to obtain the input matrix PreliminaryDataFBA if user wishes to provide different data.
% It requires an excel sheet as input with all relevant experimental data
% along with installation of Adaptive regression Splines toolbox for MATLAB
% (Jekabsons, G. (2013)).

%%

clear
clc

%% User Inputs %
MARS_dt = input("Please enter the time-step for MARS: \n");
Extr_dt = input("Please enter the time-step you want to extract from MARS: \n");
Final_Data_Points = input("Please enter the final number of data points desired: \n");

disp("Please select the directory to save PreliminaryDataFBA.mat")
filepath = uigetdir("","Please select the directory to save PreliminaryDataFBA.mat");
filename = fullfile(filepath,"PreliminaryDataFBA.mat");

%% Import Initial Data %
Read_Excel_Data = input("Would you like to read the the Excel Data file again? (1/0) \n");
if Read_Excel_Data == 1
    disp("Select the excel data file...")
    Read_Initial_Data;
else
    disp("Reading existing .MAT files for initial data from this folder...")
    load('Initial_Data.mat');
    disp("Initial Data Loaded...");
end

%% Global Variables %

% OD to gDCW/L Conversion Factor - 'ODCF' %
ODCF = 0.37; 

F = 0.8; % mL/h
V0 = 100; % mL
Feed_Time = 6; % h

Names_TAbs = Expt_Data_Culture_TAbs_1.Properties.VariableNames(2:end);
Names_CAbs = Expt_Data_Culture_CAbs_1.Properties.VariableNames(2:end);
Names_Rel = Expt_Data_Culture_Rel_1.Properties.VariableNames(2:end);

MolWt_Gly = 92;
MolWt_eYFP = 26976.57;

%% Calculating Absolute Values for Converted Absolute Metabolites for all replicates %
disp("Calculating Absolute Values for Converted Absolute Metabolites for all replicates...");

Expt_Culture_Conc_CAbs_1 = FindAbsConc(Expt_Data_Culture_CAbs_1{:,2:end},Std_Curve_Data{:,2:end}); % mM
Expt_NCC_Conc_CAbs_1 = FindAbsConc(Expt_Data_NCC_CAbs_1{:,2:end},Std_Curve_Data{:,2:end}); % mM

Expt_Culture_Conc_CAbs_2 = FindAbsConc(Expt_Data_Culture_CAbs_2{:,2:end},Std_Curve_Data{:,2:end}); % mM
Expt_NCC_Conc_CAbs_2 = FindAbsConc(Expt_Data_NCC_CAbs_2{:,2:end},Std_Curve_Data{:,2:end}); % mM

%% Calculating Mean and Standard Error values for all replicates %
disp("Calculating Mean and Standard Error values for all replicates...");

% Mean %
Expt_Culture_Conc_TAbs = mean(cat(3,Expt_Data_Culture_TAbs_1{:,2:end},Expt_Data_Culture_TAbs_2{:,2:end}),3); 
Expt_NCC_Conc_TAbs = mean(cat(3,Expt_Data_NCC_TAbs_1{:,2:end},Expt_Data_NCC_TAbs_2{:,2:end}),3); 

Expt_Culture_Conc_CAbs = mean(cat(3,Expt_Culture_Conc_CAbs_1,Expt_Culture_Conc_CAbs_2),3); % mM
Expt_NCC_Conc_CAbs = mean(cat(3,Expt_NCC_Conc_CAbs_1,Expt_NCC_Conc_CAbs_2),3); % mM

Expt_Culture_Conc_Rel = mean(cat(3,Expt_Data_Culture_Rel_1{:,2:end},Expt_Data_Culture_Rel_2{:,2:end}),3);
Expt_NCC_Conc_Rel = mean(cat(3,Expt_Data_NCC_Rel_1{:,2:end},Expt_Data_NCC_Rel_2{:,2:end}),3);

% Standard Error %
Expt_Culture_Conc_TAbs_StdErr = std(cat(3,Expt_Data_Culture_TAbs_1{:,2:end},Expt_Data_Culture_TAbs_2{:,2:end}),1,3)/sqrt(2);
Expt_NCC_Conc_TAbs_StdErr = std(cat(3,Expt_Data_NCC_TAbs_1{:,2:end},Expt_Data_NCC_TAbs_2{:,2:end}),1,3)/sqrt(2);

Expt_Culture_Conc_CAbs_StdErr = std(cat(3,Expt_Culture_Conc_CAbs_1,Expt_Culture_Conc_CAbs_2),1,3)/sqrt(2); % mM
Expt_NCC_Conc_CAbs_StdErr = std(cat(3,Expt_NCC_Conc_CAbs_1,Expt_NCC_Conc_CAbs_2),1,3)/sqrt(2); % mM

Expt_Culture_Conc_Rel_StdErr = std(cat(3,Expt_Data_Culture_Rel_1{:,2:end},Expt_Data_Culture_Rel_2{:,2:end}),1,3)/sqrt(2);
Expt_NCC_Conc_Rel_StdErr = std(cat(3,Expt_Data_NCC_Rel_1{:,2:end},Expt_Data_NCC_Rel_2{:,2:end}),1,3)/sqrt(2);

%% Initializing Experimental Time points for various datasets %
disp("Initializing Experimental Time points for various datasets...");

Expt_Culture_Time_TAbs = Expt_Data_Culture_TAbs_2.(1);
Expt_Culture_Time_TAbs(end) = 48;

Expt_NCC_Time_TAbs = Expt_Data_NCC_TAbs_1.(1);

Expt_Culture_Time_CAbs = Expt_Data_Culture_CAbs_1.(1);
Expt_NCC_Time_CAbs = Expt_Data_NCC_CAbs_1.(1);

Expt_Culture_Time_Rel = Expt_Data_Culture_Rel_1.(1);
Expt_NCC_Time_Rel = Expt_Data_NCC_Rel_1.(1);

%% Mf and Kd Determination %
disp("Determining Mf (feed concentrations) and Kd (degradation constants)....");

MfKd0 = rand(2,1);

Mf_LB = 0;
Mf_UB = inf;

Kd_LB = 0; %Bounds changed 
Kd_UB = inf;

% True Absolute Metabolites %
% M0_TAbs = Expt_NCC_Conc_TAbs(1,:);
% 
% if Mf_Data_TAbs_Available == 1
%     Mf_TAbs = Mf_Data_TAbs{1,2:end};
%     for i = 1:size(Expt_NCC_Conc_TAbs,2)
%         Kd_TAbs(i) = lsqcurvefit(@(p,t) FindKd(p, t, M0_TAbs(i), V0, F, Mf_TAbs(i), Feed_Time), MfKd0(2), Expt_NCC_Time_TAbs, Expt_NCC_Conc_TAbs(:,i), Kd_LB, Kd_UB); % h-1
%     end
% elseif Mf_Data_TAbs_Available == 0 
%     for i = 1:size(Expt_NCC_Conc_TAbs,2)
%         MfKd_TAbs(:,i) = lsqcurvefit(@(p,t) FindMfKd(p, t, M0_TAbs(i), V0, F, Feed_Time), MfKd0, Expt_NCC_Time_TAbs, Expt_NCC_Conc_TAbs(:,i), [Mf_LB Kd_LB], [Mf_UB Kd_UB]); % h-1
%         Mf_TAbs(i) = MfKd_TAbs(1,i);
%         Kd_TAbs(i) = MfKd_TAbs(2,i);
%     end
% end

Mf_TAbs = Mf_Data_TAbs{1,2:end};
Kd_TAbs = zeros(size(Mf_TAbs));
MfKd_Data_TAbs = array2table(["Mf" Mf_TAbs; "Kd" Kd_TAbs],"VariableNames",["MfKd" Names_TAbs]);
writetable(MfKd_Data_TAbs,excelfilename,"Sheet","MfKd TAbs Results");

% Converted Absolute Metabolites %
M0_CAbs = Expt_NCC_Conc_CAbs(1,:);

if Mf_Data_CAbs_Available == 1
    Mf_CAbs = Mf_Data_CAbs{1,2:end};
    for i = 1:size(Expt_NCC_Conc_CAbs,2)
        Kd_CAbs(i) = lsqcurvefit(@(p,t) FindKd(p, t, M0_CAbs(i), V0, F, Mf_CAbs(i), Feed_Time), MfKd0(2), Expt_NCC_Time_CAbs, Expt_NCC_Conc_CAbs(:,i), Kd_LB, Kd_UB); % h-1
    end
elseif Mf_Data_CAbs_Available == 0 
    for i = 1:size(Expt_NCC_Conc_CAbs,2)
        MfKd_CAbs(:,i) = lsqcurvefit(@(p,t) FindMfKd(p, t, M0_CAbs(i), V0, F, Feed_Time), MfKd0, Expt_NCC_Time_CAbs, Expt_NCC_Conc_CAbs(:,i), [Mf_LB Kd_LB], [Mf_UB Kd_UB]); % h-1
        Mf_CAbs(i) = MfKd_CAbs(1,i);
        Kd_CAbs(i) = MfKd_CAbs(2,i);
    end
end

MfKd_Data_CAbs = array2table(["Mf" Mf_CAbs; "Kd" Kd_CAbs],"VariableNames",["MfKd" Names_CAbs]);
writetable(MfKd_Data_CAbs,excelfilename,"Sheet","MfKd CAbs Results");

% Relative Metabolites %
M0_Rel = Expt_NCC_Conc_Rel(1,:);

for i = 1:size(Expt_NCC_Conc_Rel,2)
    MfKd_Rel(:,i) = lsqcurvefit(@(p,t) FindMfKd(p, t, M0_Rel(i), V0, F, Feed_Time), MfKd0, Expt_NCC_Time_Rel, Expt_NCC_Conc_Rel(:,i), [Mf_LB Kd_LB], [Mf_UB Kd_UB]); % h-1
    Mf_Rel(i) = MfKd_Rel(1,i);
    Kd_Rel(i) = MfKd_Rel(2,i);
end

MfKd_Data_Rel = array2table(["Mf" Mf_Rel; "Kd" Kd_Rel],"VariableNames",["MfKd" Names_Rel]);
writetable(MfKd_Data_Rel,excelfilename,"Sheet","MfKd Rel Results");

% % Testing MfKd Fit %
Test_Time = 0:MARS_dt:Expt_Culture_Time_TAbs(end);
for i = 1:length(Mf_CAbs)
    Test_CAbs(:,i) = TestMfKd([Mf_CAbs(i), Kd_CAbs(i)],Test_Time,M0_CAbs(i),V0,F,Feed_Time);
end
for i = 1:length(Mf_Rel)
    Test_Rel(:,i) = TestMfKd([Mf_Rel(i), Kd_Rel(i)],Test_Time,M0_Rel(i),V0,F,Feed_Time);
end

%% MARS Interpolation for Experimental Culture Data %
disp("Performing MARS Interpolation for Experimental Culture Data...");

[MARS_OD,MARS_OD_LB,MARS_OD_UB] = MARS(MARS_dt,Expt_Culture_Time_TAbs,0,48,Expt_Data_Culture_TAbs_1.OD,Expt_Data_Culture_TAbs_2.OD,"OD",2); % OD600
[MARS_eYFP,MARS_eYFP_LB,MARS_eYFP_UB] = MARS(MARS_dt,Expt_Culture_Time_TAbs,12,48,Expt_Data_Culture_TAbs_1.eYFP,Expt_Data_Culture_TAbs_2.eYFP,"eYFP",2); % g/L
[MARS_TAbs,MARS_TAbs_LB,MARS_TAbs_UB] = MARS(MARS_dt,Expt_Culture_Time_TAbs,0,48,Expt_Data_Culture_TAbs_1{:,4:end},Expt_Data_Culture_TAbs_2{:,4:end},Names_TAbs(3:end),2); % g/L
[MARS_CAbs,MARS_CAbs_LB,MARS_CAbs_UB] = MARS(MARS_dt,Expt_Culture_Time_CAbs,0,48,Expt_Culture_Conc_CAbs_1,Expt_Culture_Conc_CAbs_2,Names_CAbs,2); %mM
[MARS_Rel,MARS_Rel_LB,MARS_Rel_UB] = MARS(MARS_dt,Expt_Culture_Time_Rel,0,48,Expt_Data_Culture_Rel_1{:,2:end},Expt_Data_Culture_Rel_2{:,2:end},Names_Rel,2); 

MARS_Data_Time = [0:MARS_dt:48].';

MARS_eYFP = [zeros((length(MARS_Data_Time) - length(MARS_eYFP)),1); MARS_eYFP];
MARS_eYFP_LB = [zeros((length(MARS_Data_Time) - length(MARS_eYFP_LB)),1); MARS_eYFP_LB];
MARS_eYFP_UB = [zeros((length(MARS_Data_Time) - length(MARS_eYFP_UB)),1); MARS_eYFP_UB];

MARS_Data_OD = array2table([MARS_Data_Time MARS_OD MARS_OD_LB MARS_OD_UB],'VariableNames',["Time" Names_TAbs(1) Names_TAbs(1)+" LB" Names_TAbs(1)+" UB"]);
writetable(MARS_Data_OD,excelfilename,"Sheet","MARS OD");

MARS_Data_eYFP = array2table([MARS_Data_Time MARS_eYFP MARS_eYFP_LB MARS_eYFP_UB],'VariableNames',["Time" Names_TAbs(2) Names_TAbs(2)+" LB" Names_TAbs(2)+" UB"]);
writetable(MARS_Data_eYFP,excelfilename,"Sheet","MARS eYFP");

MARS_Data_TAbs = array2table([MARS_Data_Time MARS_TAbs],'VariableNames',["Time" Names_TAbs(3:end)]);
MARS_Data_TAbs_LB = array2table([MARS_Data_Time MARS_TAbs_LB],'VariableNames',["Time" Names_TAbs(3:end)]);
MARS_Data_TAbs_UB = array2table([MARS_Data_Time MARS_TAbs_UB],'VariableNames',["Time" Names_TAbs(3:end)]);
writetable(MARS_Data_TAbs,excelfilename,"Sheet","MARS TAbs");
writetable(MARS_Data_TAbs_LB,excelfilename,"Sheet","MARS TAbs LB");
writetable(MARS_Data_TAbs_UB,excelfilename,"Sheet","MARS TAbs UB");

MARS_Data_CAbs = array2table([MARS_Data_Time MARS_CAbs],'VariableNames',["Time" Names_CAbs]);
MARS_Data_CAbs_LB = array2table([MARS_Data_Time MARS_CAbs_LB],'VariableNames',["Time" Names_CAbs]);
MARS_Data_CAbs_UB = array2table([MARS_Data_Time MARS_CAbs_UB],'VariableNames',["Time" Names_CAbs]);
writetable(MARS_Data_CAbs,excelfilename,"Sheet","MARS CAbs");
writetable(MARS_Data_CAbs_LB,excelfilename,"Sheet","MARS CAbs LB");
writetable(MARS_Data_CAbs_UB,excelfilename,"Sheet","MARS CAbs UB");

MARS_Data_Rel = array2table([MARS_Data_Time MARS_Rel],'VariableNames',["Time" Names_Rel]);
MARS_Data_Rel_LB = array2table([MARS_Data_Time MARS_Rel_LB],'VariableNames',["Time" Names_Rel]);
MARS_Data_Rel_UB = array2table([MARS_Data_Time MARS_Rel_UB],'VariableNames',["Time" Names_Rel]);
writetable(MARS_Data_Rel,excelfilename,"Sheet","MARS Rel");
writetable(MARS_Data_Rel_LB,excelfilename,"Sheet","MARS Rel LB");
writetable(MARS_Data_Rel_UB,excelfilename,"Sheet","MARS Rel UB");

%% 1st Filtration - Extracting Time Points from MARS Data based on time-step given by User %
disp("1st Filtration - Extracting Time Points from MARS Data based on time-step given by User");

[OD_End,OD_End_Idx] = max(MARS_OD);
MARS_Time = [0:MARS_dt:MARS_Data_Time(OD_End_Idx)].';
Extr_Time = [0:Extr_dt:MARS_Time(end)].';

% Extracted Time Index - 'ETI' %
ETI = [1:(Extr_dt/MARS_dt):length(MARS_Time)].';

% Extracted Feed Time Index - 'EFTI' %
EFTI = find(Extr_Time == Feed_Time);

Extr_OD = MARS_OD(ETI,:); % OD600
Extr_eYFP = (MARS_eYFP(ETI,:)*1000)./MolWt_eYFP; % mM
Extr_TAbs = (MARS_TAbs(ETI,:)*1000)./MolWt_Gly; % mM
Extr_CAbs = MARS_CAbs(ETI,:); % mM
Extr_Rel = MARS_Rel(ETI,:);

%% Calculating Flux Mean Value at all Extracted Time Points %
disp("Calculating Flux Mean Value at all Extracted Time Points...");

% Finding dM/dt for all metabolites at all Extracted Time Points %
Extr_dMdt_OD = FindTimeDer(Extr_OD,Extr_Time);
Extr_dMdt_eYFP = FindTimeDer(Extr_eYFP,Extr_Time);
Extr_dMdt_TAbs = FindTimeDer(Extr_TAbs,Extr_Time);
Extr_dMdt_CAbs = FindTimeDer(Extr_CAbs,Extr_Time);
Extr_dMdt_Rel = FindTimeDer(Extr_Rel,Extr_Time);

% Finding Flux for all metabolites at all Extracted Time Points (Type: 1-OD; 2-TAbs/CAbs/Rel) % 
Extr_Flux_OD = FindFlux(Extr_OD,ODCF,Extr_OD,Extr_dMdt_OD,EFTI,1,[],[],CultureVolume(V0,Feed_Time,Extr_Time,F),F); % h-1
Extr_Flux_eYFP = FindFlux(Extr_OD,ODCF,Extr_eYFP,Extr_dMdt_eYFP,EFTI,2,0,0,CultureVolume(V0,Feed_Time,Extr_Time,F),F); % mmol/gDCW.h
Extr_Flux_TAbs = FindFlux(Extr_OD,ODCF,Extr_TAbs,Extr_dMdt_TAbs,EFTI,2,Mf_TAbs(3:end),Kd_TAbs(3:end),CultureVolume(V0,Feed_Time,Extr_Time,F),F); % mmol/gDCW.h
Extr_Flux_CAbs = FindFlux(Extr_OD,ODCF,Extr_CAbs,Extr_dMdt_CAbs,EFTI,2,Mf_CAbs,Kd_CAbs,CultureVolume(V0,Feed_Time,Extr_Time,F),F); % mmol/gDCW.h
Extr_Flux_Rel = FindFlux(Extr_OD,ODCF,Extr_Rel,Extr_dMdt_Rel,EFTI,2,Mf_Rel,Kd_Rel,CultureVolume(V0,Feed_Time,Extr_Time,F),F);

%% Calculating Flux Standard Deviation at all Extracted Time Points %
disp("Calculating Flux Standard Deviation at all Extracted Time Points...");

% % Finding All Metabolites Upper and Lower Bounds at all Extracted Time Points %
Extr_OD_LB = MARS_OD_LB(ETI,:); % OD600
Extr_OD_UB = MARS_OD_UB(ETI,:); % OD600

Extr_eYFP_LB = (MARS_eYFP_LB(ETI,:)*1000)./MolWt_eYFP; % mM
Extr_eYFP_UB = (MARS_eYFP_UB(ETI,:)*1000)./MolWt_eYFP; % mM

Extr_TAbs_LB = (MARS_TAbs_LB(ETI,:)*1000)./MolWt_Gly; % mM
Extr_TAbs_UB = (MARS_TAbs_UB(ETI,:)*1000)./MolWt_Gly; % mM

Extr_CAbs_LB = MARS_CAbs_LB(ETI,:); % mM
Extr_CAbs_UB = MARS_CAbs_UB(ETI,:); % mM

Extr_Rel_LB = MARS_Rel_LB(ETI,:);
Extr_Rel_UB = MARS_Rel_UB(ETI,:);
% 
% % Finding All Metabolites Concentration Standard Deviation at all Extracted Time Points %
% Extr_OD_StdDev = 0;
% Extr_TAbs_StdDev = (Extr_TAbs_UB - Extr_TAbs_LB)./2;
% Extr_CAbs_StdDev = (Extr_CAbs_UB - Extr_CAbs_LB)./2;
% Extr_Rel_StdDev = (Extr_Rel_UB - Extr_Rel_LB)./2;

% % Finding All Metabolites dq/dM at all Extracted Time Points %
% Extr_dqdt_OD = FindTimeDer(Extr_Flux_OD,Extr_Time);
% Extr_dqdt_TAbs = FindTimeDer(Extr_Flux_TAbs,Extr_Time);
% Extr_dqdt_CAbs = FindTimeDer(Extr_Flux_CAbs,Extr_Time);
% Extr_dqdt_Rel = FindTimeDer(Extr_Flux_Rel,Extr_Time);
% 
% Extr_dqdM_OD = Extr_dqdt_OD./Extr_dMdt_OD;
% Extr_dqdM_TAbs = Extr_dqdt_TAbs./Extr_dMdt_TAbs;
% Extr_dqdM_CAbs = Extr_dqdt_CAbs./Extr_dMdt_CAbs;
% Extr_dqdM_Rel = Extr_dqdt_Rel./Extr_dMdt_Rel;

% Finding All Metabolites Flux Standard Deviation at all Extracted Time Points %
Extr_Flux_OD_StdDev = 0;
Extr_Flux_eYFP_StdDev = 0.1*abs(Extr_Flux_eYFP);
Extr_Flux_TAbs_StdDev = 0.1*abs(Extr_Flux_TAbs); % abs(Extr_dqdM_TAbs.*Extr_TAbs_StdDev);
Extr_Flux_CAbs_StdDev = 0.1*abs(Extr_Flux_CAbs); % abs(Extr_dqdM_CAbs.*Extr_CAbs_StdDev);
Extr_Flux_Rel_StdDev = 0.1*abs(Extr_Flux_Rel); % abs(Extr_dqdM_Rel.*Extr_Rel_StdDev);

% % Fixing Infinity Value Flux Upper and Lower Bounds for All Metabolites at all Extracted Time Points %
% Extr_Flux_OD_StdDev = FixInf(Extr_Flux_OD_StdDev,Extr_Time,V0,Feed_Time,ODCF,Extr_OD,0,Extr_OD_StdDev,F);
% Extr_Flux_TAbs_StdDev = FixInf(Extr_Flux_TAbs_StdDev,Extr_Time,V0,Feed_Time,ODCF,Extr_OD,0,Extr_TAbs_StdDev,F);
% Extr_Flux_CAbs_StdDev = FixInf(Extr_Flux_CAbs_StdDev,Extr_Time,V0,Feed_Time,ODCF,Extr_OD,Kd_CAbs,Extr_CAbs_StdDev,F);
% Extr_Flux_Rel_StdDev = FixInf(Extr_Flux_Rel_StdDev,Extr_Time,V0,Feed_Time,ODCF,Extr_OD,Kd_Rel,Extr_Rel_StdDev,F);

% Finding Flux Upper and Lower Bounds for All Metabolites at all Extracted Time Points %
Extr_Flux_OD_LB = Extr_Flux_OD + Extr_Flux_OD_StdDev;
Extr_Flux_OD_UB = Extr_Flux_OD - Extr_Flux_OD_StdDev;

Extr_Flux_eYFP_LB = Extr_Flux_eYFP + Extr_Flux_eYFP_StdDev;
Extr_Flux_eYFP_UB = Extr_Flux_eYFP - Extr_Flux_eYFP_StdDev;

Extr_Flux_TAbs_LB = Extr_Flux_TAbs - Extr_Flux_TAbs_StdDev;
Extr_Flux_TAbs_UB = Extr_Flux_TAbs + Extr_Flux_TAbs_StdDev;

Extr_Flux_CAbs_LB = Extr_Flux_CAbs - Extr_Flux_CAbs_StdDev;
Extr_Flux_CAbs_UB = Extr_Flux_CAbs + Extr_Flux_CAbs_StdDev;

Extr_Flux_Rel_LB = Extr_Flux_Rel - Extr_Flux_Rel_StdDev;
Extr_Flux_Rel_UB = Extr_Flux_Rel + Extr_Flux_Rel_StdDev;

%% 2nd Filtration - Creating Equally Spaced OD Intervals %
disp("2nd Filtration - Creating Equally Spaced OD Intervals...");

OD_Start = 0.05;
OD_Step = (OD_End - OD_Start)/Final_Data_Points;
OD = [OD_Start:OD_Step:OD_End].';

% Check if User given filtration is feasible %
if length(OD) > length(Extr_Time)
    disp("Please Decrease Time-Step or Final Number of Data Points");
    return;
end

% Final Time Index - 'FTI' %
[Time,FTI] = FindClosestTime(OD,Extr_OD,Extr_Time); 
Time = Time.';
FTI = FTI.';

FTI(1) = find(double(Extr_Flux_OD>0 & Extr_Flux_OD<1.5)>0,1,'first');
Time(1) = Extr_Time(FTI(1));

FTI(1) = 101;
Time(1) = 1;



OD = Extr_OD(FTI);
OD_LB = Extr_OD_LB(FTI);
OD_UB = Extr_OD_UB(FTI);

eYFP = Extr_eYFP(FTI);
eYFP_LB = Extr_eYFP_LB(FTI);
eYFP_UB = Extr_eYFP_UB(FTI);

TAbs = Extr_TAbs(FTI,:);
TAbs_LB = Extr_TAbs_LB(FTI,:);
TAbs_UB = Extr_TAbs_UB(FTI,:);

CAbs = Extr_CAbs(FTI,:);
CAbs_LB = Extr_CAbs_LB(FTI,:);
CAbs_UB = Extr_CAbs_UB(FTI,:);

Rel = Extr_Rel(FTI,:);
Rel_LB = Extr_Rel_LB(FTI,:);
Rel_UB = Extr_Rel_UB(FTI,:);

Flux_OD = Extr_Flux_OD(FTI);
Flux_OD_LB = Extr_Flux_OD_LB(FTI);
Flux_OD_UB = Extr_Flux_OD_UB(FTI);

Flux_eYFP = Extr_Flux_eYFP(FTI);
Flux_eYFP_LB = Extr_Flux_eYFP_LB(FTI);
Flux_eYFP_UB = Extr_Flux_eYFP_UB(FTI);

Flux_TAbs = Extr_Flux_TAbs(FTI,:);
Flux_TAbs_LB = Extr_Flux_TAbs_LB(FTI,:);
Flux_TAbs_UB = Extr_Flux_TAbs_UB(FTI,:);

Flux_CAbs = Extr_Flux_CAbs(FTI,:);
Flux_CAbs_LB = Extr_Flux_CAbs_LB(FTI,:);
Flux_CAbs_UB = Extr_Flux_CAbs_UB(FTI,:);

Flux_Rel = Extr_Flux_Rel(FTI,:);
Flux_Rel_LB = Extr_Flux_Rel_LB(FTI,:);
Flux_Rel_UB = Extr_Flux_Rel_UB(FTI,:);

%% Final Flux Matrices % Tables %
disp("Creating Final Flux Matrices and Tables...");

Final_Flux_Abs_Names = ["OD Interval" "Time" Names_TAbs Names_CAbs];
Final_Flux_Abs = [OD Time Flux_OD Flux_eYFP Flux_TAbs Flux_CAbs];
Final_Flux_Abs_LB = [OD Time Flux_OD_LB Flux_eYFP_LB Flux_TAbs_LB Flux_CAbs_LB];
Final_Flux_Abs_UB = [OD Time Flux_OD_UB Flux_eYFP_UB Flux_TAbs_UB Flux_CAbs_UB];

Final_Flux_Rel_Names = ["OD Interval" "Time" Names_Rel];
Final_Flux_Rel = [OD Time Flux_Rel];
Final_Flux_Rel_LB = [OD Time Flux_Rel_LB];
Final_Flux_Rel_UB = [OD Time Flux_Rel_UB];

%% Calculating Oxygen Flux if DO data is available %

if Expt_Data_DO_Available == 1
    disp("Calculating Oxygen Flux if DO data is available....");
    MolWt_DO = 32;
    kla = 133;
    C_Oxy = (0.0075/MolWt_DO)*1000; % mM
    
    Expt_Culture_Conc_DO = mean(cat(2,Expt_Data_Culture_DO_1{:,2:end},Expt_Data_Culture_DO_2{:,2:end}),2); 
    Expt_Culture_Conc_DO_StdErr = std(cat(2,Expt_Data_Culture_DO_1{:,2:end},Expt_Data_Culture_DO_2{:,2:end}),1,2)/sqrt(2);
    
    Expt_Culture_Time_DO = Expt_Data_Culture_DO_1.(1);
    
    [Extr_Time_DO,ETI_DO] = FindClosestTime(Extr_Time,Expt_Culture_Time_DO,Expt_Culture_Time_DO);
    Extr_Time_DO = Extr_Time_DO.';
    ETI_DO = ETI_DO.';

    Extr_DO = Expt_Culture_Conc_DO(ETI_DO);

    Extr_dMdt_DO = FindTimeDer(Extr_DO,Extr_Time);

    Extr_Flux_DO = FindFluxDO(Extr_OD,ODCF,Extr_DO,Extr_dMdt_DO,EFTI,kla,C_Oxy,CultureVolume(V0,Feed_Time,Extr_Time,F),F); % mmol/gDCW.h

    Extr_DO_LB = Extr_DO - Expt_Culture_Conc_DO_StdErr(ETI,:); 
    Extr_DO_UB = Extr_DO + Expt_Culture_Conc_DO_StdErr(ETI,:); 
% 
%     Extr_DO_StdDev = (Extr_DO_UB - Extr_DO_LB)./2;
% 
%     Extr_dqdt_DO = FindTimeDer(Extr_Flux_DO,Extr_Time);
% 
%     Extr_dqdM_DO = Extr_dqdt_DO./Extr_dMdt_DO;

    Extr_Flux_DO_StdDev = 0.1*abs(Extr_Flux_DO); % abs(Extr_dqdM_DO.*Extr_DO_StdDev);
%     Extr_Flux_DO_StdDev = FixInf(Extr_Flux_DO_StdDev,Extr_Time,V0,Feed_Time,ODCF,Extr_OD,0,Extr_DO_StdDev,F);

    Extr_Flux_DO_LB = Extr_Flux_DO - Extr_Flux_DO_StdDev;
    Extr_Flux_DO_UB = Extr_Flux_DO + Extr_Flux_DO_StdDev;

    DO = Extr_DO(FTI);
    DO_LB = Extr_DO_LB(FTI);
    DO_UB = Extr_DO_UB(FTI);

    Flux_DO = Extr_Flux_DO(FTI);
    Flux_DO_LB = Extr_Flux_DO_LB(FTI);
    Flux_DO_UB = Extr_Flux_DO_UB(FTI);

    Final_Flux_Abs_Names = ["OD Interval" "Time" Names_TAbs "DO Flux" Names_CAbs];
    Final_Flux_Abs = [OD Time Flux_OD Flux_eYFP Flux_TAbs Flux_DO Flux_CAbs];
    Final_Flux_Abs_LB = [OD Time Flux_OD_LB Flux_eYFP_LB Flux_TAbs_LB Flux_DO_LB Flux_CAbs_LB];
    Final_Flux_Abs_UB = [OD Time Flux_OD_UB Flux_eYFP_UB Flux_TAbs_UB Flux_DO_UB Flux_CAbs_UB];
end

%% Removing Extra Metabolites %
disp("Removing Extra Metabolites that are not the part of model...");
Model_Present_Abs = find(IndFlux_Data{1,2:end} ~= 0);
Model_Present_Rel = find(IndFluxr_Data{1,2:end} ~= 0);

Final_Flux_Data_Abs = array2table(Final_Flux_Abs(:,Model_Present_Abs+2),'VariableNames',Final_Flux_Abs_Names(:,Model_Present_Abs+2));
Final_Flux_Data_Abs_LB = array2table(Final_Flux_Abs_LB(:,Model_Present_Abs+2),'VariableNames',Final_Flux_Abs_Names(:,Model_Present_Abs+2));
Final_Flux_Data_Abs_UB = array2table(Final_Flux_Abs_UB(:,Model_Present_Abs+2),'VariableNames',Final_Flux_Abs_Names(:,Model_Present_Abs+2));
writetable(Final_Flux_Data_Abs,excelfilename,'Sheet',"Final Flux Abs");
writetable(Final_Flux_Data_Abs_LB,excelfilename,'Sheet',"Final Flux Abs LB");
writetable(Final_Flux_Data_Abs_UB,excelfilename,'Sheet',"Final Flux Abs UB");

Final_Flux_Data_Rel = array2table(Final_Flux_Rel(:,Model_Present_Rel+2),'VariableNames',Final_Flux_Rel_Names(:,Model_Present_Rel+2));
Final_Flux_Data_Rel_LB = array2table(Final_Flux_Rel_LB(:,Model_Present_Rel+2),'VariableNames',Final_Flux_Rel_Names(:,Model_Present_Rel+2));
Final_Flux_Data_Rel_UB = array2table(Final_Flux_Rel_UB(:,Model_Present_Rel+2),'VariableNames',Final_Flux_Rel_Names(:,Model_Present_Rel+2));
writetable(Final_Flux_Data_Rel,excelfilename,'Sheet',"Final Flux Rel");
writetable(Final_Flux_Data_Rel_LB,excelfilename,'Sheet',"Final Flux Rel LB");
writetable(Final_Flux_Data_Rel_UB,excelfilename,'Sheet',"Final Flux Rel UB");

FC1 = MARS_Rel(1,:);
FC2 = max(MARS_Rel,[],1);
FC = FC2./FC1;
FC = FC(Model_Present_Rel);

MFC = Mf_Rel./FC1;
MFC = MFC(Model_Present_Rel);

V = zeros(size(Time));
FI = find(Time < Feed_Time, 1, 'last' );
V(1:FI,1) = V0;
V(FI+1:end,1) = V0 + F*(Time(FI+1:end)-Feed_Time);
F_V = zeros(size(Final_Flux_Rel,1),1);
F_V(1:FI,1) = 0; 
F_V(FI+1:end,1) =F./V(FI+1:end,1);

%% Creating Preliminary Data.mat File %
disp("Creating Preliminary Data.mat File for FBA calculations...");

load("iECD_1391.mat");
modelpFBA = iECD_1391;

timesFBA = Time';
ODist = OD';

meas_flux = table2array(Final_Flux_Data_Abs);
fluxm_lb = table2array(Final_Flux_Data_Abs_LB);
fluxm_ub = table2array(Final_Flux_Data_Abs_UB);
ind_flux = IndFlux_Data{1,Model_Present_Abs+1};

meas_fluxr = table2array(Final_Flux_Data_Rel);
ind_fluxr = IndFluxr_Data{1,Model_Present_Rel+1};
C_rel = Rel(:,Model_Present_Rel);
exchrxn = double(findExcRxns(modelpFBA));

prefluxr;

save(filename,'exchrxn','fluxm_lb','fluxm_ub','ind_flux','ind_fluxr','meas_flux','meas_fluxr','modelpFBA','ODist','pre_fluxr','timesFBA','FC','MFC','F_V','C_rel');
save("FBA_Model_Input_Data.mat");

% %% Making Presentations of all Plots %
% disp("Creating Presentations of plots of Metabolite concentrations and Fluxes...");
% 
% % X-Label Type = 1-Time(h); 2-OD(600) %
% % Y-Label Type = 1-Concentration; 2-Normalized Area Ratios; 3-Abs Flux; 4-Rel Flux %
% % Expt_Data = 1-Available; 0-Not Available %
% 
% % % MfKd Fit %
% PPT_MfKd_CAbs = MakePPT("PPT_MfKd_CAbs",Names_CAbs,1,1,1,Test_Time,Test_CAbs,zeros(size(Test_CAbs)),zeros(size(Test_CAbs)),Expt_NCC_Time_CAbs,zeros(size(Expt_NCC_Time_CAbs)),Expt_NCC_Conc_CAbs,Expt_NCC_Conc_CAbs_StdErr);
% PPT_MfKd_Rel = MakePPT("PPT_MfKd_Rel",Names_Rel,1,2,1,Test_Time,Test_Rel,zeros(size(Test_Rel)),zeros(size(Test_Rel)),Expt_NCC_Time_Rel,zeros(size(Expt_NCC_Time_Rel)),Expt_NCC_Conc_Rel,Expt_NCC_Conc_Rel_StdErr);
% 
% % % MARS %
% PPT_MARS_OD = MakePPT("PPT_MARS_OD","OD",1,1,1,MARS_Data_Time,MARS_OD,MARS_OD_LB,MARS_OD_UB,Expt_Culture_Time_TAbs,zeros(size(Expt_Culture_Time_TAbs)),Expt_Culture_Conc_TAbs(:,1),Expt_Culture_Conc_TAbs_StdErr(:,1));
% PPT_MARS_eYFP = MakePPT("PPT_MARS_eYFP","eYFP",1,1,1,MARS_Data_Time,MARS_eYFP,MARS_eYFP_LB,MARS_eYFP_UB,Expt_Culture_Time_TAbs,zeros(size(Expt_Culture_Time_TAbs)),Expt_Culture_Conc_TAbs(:,2),Expt_Culture_Conc_TAbs_StdErr(:,2));
% PPT_MARS_TAbs = MakePPT("PPT_MARS_TAbs",Names_TAbs(3:end),1,1,1,MARS_Data_Time,MARS_TAbs,MARS_TAbs_LB,MARS_TAbs_UB,Expt_Culture_Time_TAbs,zeros(size(Expt_Culture_Time_TAbs)),Expt_Culture_Conc_TAbs(:,3:end),Expt_Culture_Conc_TAbs_StdErr(:,3:end));
% PPT_MARS_CAbs = MakePPT("PPT_MARS_CAbs",Names_CAbs,1,1,1,MARS_Data_Time,MARS_CAbs,MARS_CAbs_LB,MARS_CAbs_UB,Expt_Culture_Time_CAbs,[],Expt_Culture_Conc_CAbs,Expt_Culture_Conc_CAbs_StdErr);
% PPT_MARS_Rel = MakePPT("PPT_MARS_Rel",Names_Rel,1,2,1,MARS_Data_Time,MARS_Rel,MARS_Rel_LB,MARS_Rel_UB,Expt_Culture_Time_Rel,[],Expt_Culture_Conc_Rel,Expt_Culture_Conc_Rel_StdErr);
% 
% % % Metabolite vs OD %
% PPT_TAbs_eYFP = MakePPT("PPT_eYFP_OD","eYFP",2,1,1,OD,eYFP,eYFP_LB,eYFP_UB,Expt_Culture_Conc_TAbs(:,1),Expt_Culture_Conc_TAbs_StdErr(:,1),Expt_Culture_Conc_TAbs(:,2),Expt_Culture_Conc_TAbs_StdErr(:,2));
% PPT_TAbs_OD = MakePPT("PPT_TAbs_OD",Names_TAbs(3:end),2,1,1,OD,TAbs,TAbs_LB,TAbs_UB,Expt_Culture_Conc_TAbs(:,1),Expt_Culture_Conc_TAbs_StdErr(:,1),Expt_Culture_Conc_TAbs(:,3:end),Expt_Culture_Conc_TAbs_StdErr(:,3:end));
% PPT_CAbs_OD = MakePPT("PPT_CAbs_OD",Names_CAbs,2,1,1,OD,CAbs,CAbs_LB,CAbs_UB,Expt_Culture_Conc_TAbs([1,3,5,10,14,15,end],1),Expt_Culture_Conc_TAbs_StdErr([1,3,5,10,14,15,end],1),Expt_Culture_Conc_CAbs,Expt_Culture_Conc_CAbs_StdErr);
% PPT_Rel_OD = MakePPT("PPT_Rel_OD",Names_Rel,2,2,1,OD,Rel,Rel_LB,Rel_UB,Expt_Culture_Conc_TAbs([1,3,5,10,14,15,end],1),Expt_Culture_Conc_TAbs_StdErr([1,3,5,10,14,15,end],1),Expt_Culture_Conc_Rel,Expt_Culture_Conc_Rel_StdErr);
% 
% % % Fluxes vs OD%
% PPT_Flux_OD_OD = MakePPT("PPT_Flux_OD_OD","OD Flux",2,3,0,OD,Flux_OD,Flux_OD_LB,Flux_OD_UB,0,0,0,0);
% PPT_Flux_eYFP_OD = MakePPT("PPT_Flux_eYFP_OD","eYFP Flux",2,3,0,OD,Flux_eYFP,Flux_eYFP_LB,Flux_eYFP_UB,0,0,0,0);
% PPT_Flux_TAbs_OD = MakePPT("PPT_Flux_TAbs_OD",Names_TAbs(3:end)+"Flux",2,3,0,OD,Flux_TAbs,Flux_TAbs_LB,Flux_TAbs_UB,0,0,0,0);
% PPT_Flux_CAbs_OD = MakePPT("PPT_Flux_CAbs_OD",Names_CAbs+"Flux",2,3,0,OD,Flux_CAbs,Flux_CAbs_LB,Flux_CAbs_UB,0,0,0,0);
% PPT_Flux_Rel_OD = MakePPT("PPT_Flux_Rel_OD",Names_Rel+"Flux",2,3,0,OD,Flux_Rel,Flux_Rel_LB,Flux_Rel_UB,0,0,0,0);
