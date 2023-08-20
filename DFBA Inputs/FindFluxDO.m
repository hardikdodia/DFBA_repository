function q = FindFluxDO(OD,ODCF,DO,dMdt_DO,Feed_Time_Index,kla,C_Oxy,Volume,F)
    q = zeros(size(DO));
    V = Volume;
    
    Batch = [1:Feed_Time_Index];
    FedBatch = [Feed_Time_Index+1:size(DO,1)];

    q(Batch) = (dMdt_DO(Batch) - (kla*(100 - DO(Batch))))./((100/C_Oxy)*(ODCF*OD(Batch)));
    q(FedBatch) = (dMdt_DO(FedBatch) - (kla*(100 - DO(FedBatch))) + ((F*DO(FedBatch))./V(FedBatch)))./((100/C_Oxy)*(ODCF*OD(FedBatch)));
end




