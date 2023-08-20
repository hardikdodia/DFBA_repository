%% Kd Validation %%

% run FBA_Model_Inputs.m manually until Mf-Kd Determination Section

CAbs_min = min(Expt_NCC_Conc_CAbs_1,Expt_NCC_Conc_CAbs_2);
CAbs_max = max(Expt_NCC_Conc_CAbs_1,Expt_NCC_Conc_CAbs_2);

Rel_min = min(Expt_Data_NCC_Rel_1{:,2:end},Expt_Data_NCC_Rel_2{:,2:end});
Rel_max = max(Expt_Data_NCC_Rel_1{:,2:end},Expt_Data_NCC_Rel_2{:,2:end});

%% When Kd is predicted %%
Test_CAbs_Expt = Test_CAbs([1,61,121,241,361,421,481],:);
Test_Rel_Expt = Test_Rel([1,61,241,481],:);

Score_CAbs = sum((Test_CAbs_Expt < CAbs_max & Test_CAbs_Expt > CAbs_min),1);
Score_Rel = sum((Test_Rel_Expt < Rel_max & Test_Rel_Expt > Rel_min),1);

%% When Kd = 0 %%
Kd_CAbs = zeros(size(Mf_CAbs));
Kd_Rel = zeros(size(Mf_Rel));

Test_Time = 0:MARS_dt:Expt_Culture_Time_TAbs(end);
for i = 1:length(Mf_CAbs)
    Test_CAbs2(:,i) = TestMfKd([Mf_CAbs(i), Kd_CAbs(i)],Test_Time,M0_CAbs(i),V0,F,Feed_Time);
end
for i = 1:length(Mf_Rel)
    Test_Rel2(:,i) = TestMfKd([Mf_Rel(i), Kd_Rel(i)],Test_Time,M0_Rel(i),V0,F,Feed_Time);
end

Test_CAbs2_Expt = Test_CAbs2([1,61,121,241,361,421,481],:);
Test_Rel2_Expt = Test_Rel2([1,61,241,481],:);

Score_CAbs2 = sum((Test_CAbs2_Expt < CAbs_max & Test_CAbs2_Expt > CAbs_min),1);
Score_Rel2 = sum((Test_Rel2_Expt < Rel_max & Test_Rel2_Expt > Rel_min),1);

Score_CAbs_Total = sum(Score_CAbs)/(14*7);
Score_Rel_Total = sum(Score_Rel)/(232*4);
Score_CAbs2_Total = sum(Score_CAbs2)/(14*7);
Score_Rel2_Total = sum(Score_Rel2)/(232*4);

SSE_CAbs1 = 




