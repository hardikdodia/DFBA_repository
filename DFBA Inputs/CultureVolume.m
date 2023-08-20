function V = Culture_Volume(V0,T0,T,F)
    V = zeros(size(T));
    Feed_Time_Index = find(T == T0);
    V(1:Feed_Time_Index,1) = V0;
    V(Feed_Time_Index+1:end,1) = V0 + F*(T(Feed_Time_Index+1:end)-T0);
end