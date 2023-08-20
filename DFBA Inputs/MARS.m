function [MARS_M,MARS_M_LB,MARS_M_UB] = MARS(MARS_dt, Expt_Time, t_start, t_end, Expt_Data_1, Expt_Data_2, Names_Met, N_Replicates)

[~,t_start_idx] = FindClosestTime(t_start,Expt_Time,Expt_Time);
[~,t_end_idx] = FindClosestTime(t_end,Expt_Time,Expt_Time);

Combinations = MakeCombinations((t_end_idx - t_start_idx) + 1);

MARS_Input_Data = array2table(zeros(size(Expt_Data_1,2),2+(2*((t_end_idx - t_start_idx) + 1))));
MARS_Input_Data.(2) = Names_Met.';
MARS_Input_Data{:,3:2:end} = transpose(Expt_Data_1(t_start_idx:t_end_idx,:));
MARS_Input_Data{:,4:2:end} = transpose(Expt_Data_2(t_start_idx:t_end_idx,:));

Expt_Time = Expt_Time(t_start_idx:t_end_idx);
Expt_Time_Points = (t_end_idx - t_start_idx) + 1;

MARS_Time_Start = Expt_Time(1);
MARS_Time_End = Expt_Time(end);
MARS_Time_Points = ((MARS_Time_End - MARS_Time_Start)/MARS_dt)+1;

MARS_results_mean_overall = [];
MARS_results_lb_overall = [];
MARS_results_ub_overall = [];

%% Bootstrapping
parfor m = 1:size(MARS_Input_Data,1)
    
met_values = [] ;
MARS_results = zeros(MARS_Time_Points,1000) ; % time plotted from 0 to 48 , fitted 1000 times

for i = 1:Expt_Time_Points

    met_values_temp = zeros(1,N_Replicates) ;
    met_values_temp = table2array(MARS_Input_Data(m,(2*i)+1:(2*i)+2)) ; % Scope to generalize
    met_values = [met_values ; met_values_temp] ;

end   

time_metvalues = [Expt_Time met_values] ; % time_metvalues - first column contains time values , second onwards contains corresponding levels for each replicate 
metvalues_mean = mean(time_metvalues(:,2:end),2) ;

%% shuffling within row elements , followed by random row extraction

time_metvalues_model_temp = [] ;


for j=1:1000


    b= Combinations(:, randperm(size(Combinations, 2))); % combinations matrix is a set of combinations through which elements can be extracted from each time point so that total data used is nearly 90 percent

    Y  = datasample(b,1,1,'Replace',true);

    metval =[] ;

    for i=1:Expt_Time_Points
    
        metval(i,1) = median(datasample(met_values(i,:),Y(i),'Replace',false)) ;

    end

    time_metvalues_model = [Expt_Time metval] ;

    time_metvalues_model_temp = [time_metvalues_model_temp metval] ;

%% Model building and implementation for culture vs time

    params = aresparams2('useMinSpan',1,'useEndSpan',1,'cubic',true,'c',0) ;

    [model,time, resultsEval] = aresbuild(Expt_Time,time_metvalues_model(:,2),params) ;

    MARS_results(:,j) =  arespredict(model,[MARS_Time_Start:MARS_dt:MARS_Time_End]');

end

time_metvalues_model_mean = mean(time_metvalues_model_temp,2) ;
MARS_results_mean = mean(MARS_results,2) ;
MARS_results_std = std(MARS_results,0,2) ;
MARS_results_lb = MARS_results_mean - MARS_results_std;
MARS_results_ub = MARS_results_mean + MARS_results_std;

%% Below variable stores the model results of all metabolites row wise
MARS_results_mean_overall = [MARS_results_mean_overall MARS_results_mean]; 
MARS_results_lb_overall = [MARS_results_lb_overall MARS_results_lb];
MARS_results_ub_overall = [MARS_results_ub_overall MARS_results_ub];
end

MARS_M = MARS_results_mean_overall;
MARS_M_LB = MARS_results_lb_overall;
MARS_M_UB = MARS_results_ub_overall;
end


