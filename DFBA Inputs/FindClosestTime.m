function [Time,Index] = FindClosestTime(OD,MARS_OD,MARS_Time) % all column matrices

    for i = 1:length(OD)
        [Error(i),Index(i)] = min(abs(OD(i) - MARS_OD));
        Time(i) = MARS_Time(Index(i));
    end
    
    Time(1) = 0;
    Time(end) = MARS_Time(end);
end
