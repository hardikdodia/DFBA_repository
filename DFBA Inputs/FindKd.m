function M = FindKd(p,t,M0,V0,F,Mf,Feed_Time)
    
    Feed_Time_Index = find(t == Feed_Time);

    [T_Batch,M_Batch] = ode45(@(t,x) massbalance_Batch(t,x), t(1:Feed_Time_Index), M0);

    [T_FedBatch,M_FedBatch] = ode45(@(t,x) massbalance_FedBatch(t,x), t(Feed_Time_Index:end), [M_Batch(end),V0]);

    function dM_Batch = massbalance_Batch(t,x)
        dM_Batch = zeros(1,1);
        dM_Batch(1) = p.*x(1);
    end

    function dM_FedBatch = massbalance_FedBatch(t,x) 
        dM_FedBatch = zeros(2,1);
        dM_FedBatch(1) = ((F./x(2)).*(Mf - x(1))) + (p.*x(1));
        dM_FedBatch(2) = F;
    end
    
    M = [M_Batch(1); M_Batch(end); M_FedBatch(2:end,1)];
end

