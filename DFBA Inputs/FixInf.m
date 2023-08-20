function f = fixinf(Extracted_Flux_StdDev,Extracted_Time,V0,Feed_Time,OD_gDCW_Factor_Input,Extracted_OD,K,Extracted_conc_StdDev,F)


    [r,c] = find(double(isinf(Extracted_Flux_StdDev)));
    
    if ~isempty([r,c])
    
        for i=1:size(r)
            if Extracted_Time(r(i)) > Feed_Time
        
                vol = V0+(F*(Extracted_Time(r(i))-Feed_Time));  % generalize time point
            else 
                vol = V0;
                F =0;
            end
        
            Extracted_Flux_StdDev(r(i),c(i)) = abs((1/(OD_gDCW_Factor_Input*Extracted_OD(r(i))))*((-K(c(i))+((F/vol))*((Extracted_conc_StdDev(r(i),c(i)))))));
        
           
        end
  


    end
    f = Extracted_Flux_StdDev;  
  
end