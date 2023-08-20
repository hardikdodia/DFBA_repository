function AbsConc = FindAbsConc(Expt_Data_Conc_CAbs, Std_Curve_Data)
    AbsConc = zeros(size(Expt_Data_Conc_CAbs));
    AbsConc = (((Expt_Data_Conc_CAbs - Std_Curve_Data(2,:))./Std_Curve_Data(1,:))./Std_Curve_Data(3,:))*1000; % mM
    [Idx1,Idx2] = find(AbsConc < 0);
    AbsConc(Idx1,Idx2) = 0;
end
