function TimeDer = FindTimeDer(Variable,Time)
    dt = Time(2) - Time(1);
    TimeDer = zeros(size(Variable));
    TimeDer(1,:) = ((-3*Variable(1,:)) + (4*Variable(2,:)) - (1*Variable(3,:)))./(2*dt);
    TimeDer([2:end-1],:) = (Variable([3:end],:) - Variable([1:end-2],:))./(2*dt);
    TimeDer(end,:) = ((1*Variable(end-2,:)) - (4*Variable(end-1,:)) + (3*Variable(end,:)))./(2*dt);    
end

