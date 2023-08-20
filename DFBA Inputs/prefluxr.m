Extr_dMdt_Rel_NCC = FindTimeDer(Test_Rel(ETI,:),Extr_Time);

ncc_dMdt_Rel = Extr_dMdt_Rel_NCC(FTI,:);

culture_dMdt_Rel = Extr_dMdt_Rel(FTI,:);

pre_fluxr = (double(culture_dMdt_Rel <ncc_dMdt_Rel))' ;

pre_fluxr = pre_fluxr(Model_Present_Rel,:);

save("prefluxr.mat","pre_fluxr");



