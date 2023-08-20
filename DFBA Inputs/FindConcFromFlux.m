function M = FindConcFromFlux(q,X,Time,Feed_Time,dMdt,V,F,Mf,Kd,M0)
    M = zeros(size(q));
    [Closest_Feed_Time, Closest_Feed_Time_Index] = FindClosestTime(Feed_Time,Time,Time);

    batch = [1:Closest_Feed_Time_Index];
    fedBatch = [Closest_Feed_Time_Index+1:size(M,1)];

    if Kd == 0
        M(batch) = M0 + q(batch).*(X(batch)).*(Time(batch));
    else
        M(batch) = (dMdt(batch) - (q(batch).*X(batch)))/Kd;
    end
    M(fedBatch) = (dMdt(fedBatch) - ((F*Mf)./V(fedBatch)) - (q(fedBatch).*X(fedBatch)))./(Kd - (F./V(fedBatch)));

end
