clear
clc

MARS_dt = 0.01;
Expt_Time = Expt_Culture_Time_TAbs;
Expt_Data_1 = Expt_Data_Culture_TAbs_1.OD;
Expt_Data_2 = Expt_Data_Culture_TAbs_2.OD;
Names_Met = ["OD"];
N_Replicates = 2;
Combinations = MARS_Combinations_Data.TAbs;
[MARS_M,MARS_M_LB,MARS_M_UB] = MARS(MARS_dt, Expt_Time, Expt_Data_1, Expt_Data_2, Names_Met, N_Replicates, Combinations.');
