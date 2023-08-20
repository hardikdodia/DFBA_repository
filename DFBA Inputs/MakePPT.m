function PPT = MakePPT(Title, MetNames, X_Label_Type, Y_Label_Type, Expt_Data, X, Y, Y_LB, Y_UB, X_Expt, X_Expt_StdDev, Y_Expt, Y_Expt_StdDev)

    import mlreportgen.ppt.*
    
    PPT = Presentation(Title + ".pptx");
    
    titleSlide = add(PPT,'Title Slide');
    replace(titleSlide,'Title', Title);
    
    % Assigning Axes Labels %
    if X_Label_Type == 1 
        X_Label = "Time (h)";
    elseif X_Label_Type == 2
        X_Label = "OD 600";
    end

    if Y_Label_Type == 1 
        Y_Label = "Concentration (mM)";
    elseif Y_Label_Type == 2
        Y_Label = "Normalized Area Ratios";
    elseif Y_Label_Type == 3
        Y_Label = "Absolute Flux (mmol/gDCW.h)";
    elseif Y_Label_Type == 4
        Y_Label = "Relative Flux";
    end
    
    if Expt_Data == 1
        if Y_LB == 0
            for i = 1:length(MetNames)
                f(i) = figure('Visible','off');
                plot(X,Y(:,i),'-',Color="black"); % Continuous Data - Mean Value
                title(MetNames(i));
                xlabel(X_Label);
                ylabel(Y_Label);
                hold on
                plot(X_Expt,Y_Expt(:,i),'o',Color="red"); % Experimental Data - Mean Value
                hold on
                errorbar(X_Expt,Y_Expt(:,i),X_Expt_StdDev,'horizontal',"LineStyle","none",Color="red"); % Experimental Data - X_StdDev
                hold on
                errorbar(X_Expt,Y_Expt(:,i),Y_Expt_StdDev(:,i),'vertical',"LineStyle","none",Color="red"); % Experimental Data - Y_StdDev
                hold off
                legend('Mean', 'Experimental Data'); % legend
                saveas(gcf,MetNames(i) + ".jpg");
                pictureSlide = add(PPT,'Title and Picture');
                replace(pictureSlide,'Title',MetNames(i));
                replace(pictureSlide,'Picture',Picture(MetNames(i) + ".jpg"));
            end
        else
            for i = 1:length(MetNames)
                f(i) = figure('Visible','off');
                plot(X,Y(:,i),'-',Color="black"); % Continuous Data - Mean Value
                title(MetNames(i));
                xlabel(X_Label);
                ylabel(Y_Label);
                hold on
                plot(X,Y_LB(:,i),':',Color="blue"); % Continuous Data - Lower Bound
                hold on
                plot(X,Y_UB(:,i),':',Color="blue"); % Continuous Data - Upper Bound
                hold on
                plot(X_Expt,Y_Expt(:,i),'o',Color="red"); % Experimental Data - Mean Value
                hold on
                errorbar(X_Expt,Y_Expt(:,i),X_Expt_StdDev,'horizontal',"LineStyle","none",Color="red"); % Experimental Data - X_StdDev
                hold on
                errorbar(X_Expt,Y_Expt(:,i),Y_Expt_StdDev(:,i),'vertical',"LineStyle","none",Color="red"); % Experimental Data - Y_StdDev
                hold off
                legend('Mean', 'Lower Bound', 'Upper Bound', 'Experimental Data'); % legend
                saveas(gcf,MetNames(i) + ".jpg");
                pictureSlide = add(PPT,'Title and Picture');
                replace(pictureSlide,'Title',MetNames(i));
                replace(pictureSlide,'Picture',Picture(MetNames(i) + ".jpg"));
            end
        end
    elseif Expt_Data == 0
        if Y_LB == 0
            for i = 1:length(MetNames)
                f(i) = figure('Visible','off');
                plot(X,Y(:,i),'-',Color="black"); % Continuous Data - Mean Value
                title(MetNames(i));
                xlabel(X_Label);
                ylabel(Y_Label);
                hold off
                legend('Mean'); % legend
                saveas(gcf,MetNames(i) + ".jpg");
                pictureSlide = add(PPT,'Title and Picture');
                replace(pictureSlide,'Title',MetNames(i));
                replace(pictureSlide,'Picture',Picture(MetNames(i) + ".jpg"));
            end
        else
            for i = 1:length(MetNames)
                f(i) = figure('Visible','off');
                plot(X,Y(:,i),'-',Color="black"); % Continuous Data - Mean Value
                title(MetNames(i));
                xlabel(X_Label);
                ylabel(Y_Label);
                hold on
                plot(X,Y_LB(:,i),':',Color="blue"); % Continuous Data - Lower Bound
                hold on
                plot(X,Y_UB(:,i),':',Color="blue"); % Continuous Data - Upper Bound
                hold off
                legend('Mean', 'Lower Bound', 'Upper Bound'); % legend
                saveas(gcf,MetNames(i) + ".jpg");
                pictureSlide = add(PPT,'Title and Picture');
                replace(pictureSlide,'Title',MetNames(i));
                replace(pictureSlide,'Picture',Picture(MetNames(i) + ".jpg"));
            end
        end
    end

    close(PPT);

    delete *.jpg;

end
