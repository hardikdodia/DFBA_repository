function q = FindFlux(OD,ODCF,M,dMdt,Feed_Time_Index,Type,p1,p2,Volume,F)
    q = zeros(size(M));
    V = Volume;
    
    Batch = [1:Feed_Time_Index];
    FedBatch = [Feed_Time_Index+1:size(M,1)];

    if Type == 1
        q(Batch) = dMdt(Batch)./OD(Batch);
        q(FedBatch) = (dMdt(FedBatch)./OD(FedBatch)) + (F./V(FedBatch));
   
    elseif Type == 2
        for i = 1:size(M,2)
            q(Batch,i) = (dMdt(Batch,i) + (p2(i)*(M(Batch,i))))./(ODCF*OD(Batch));% changed sign
            q(FedBatch,i) = (dMdt(FedBatch,i) + (p2(i)*(M(FedBatch,i))) - ((F*(p1(i) - M(FedBatch,i)))./V(FedBatch)))./(ODCF*OD(FedBatch));% changed sign
        end
    end
end




