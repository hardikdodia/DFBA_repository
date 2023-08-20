%% FBA Model Results %
clear
clc

%% Import Input and Results Data %%
disp("Select Model Input Data .MAT file");
FBA_Model_Input_Data_Path = uigetfile;
load(FBA_Model_Input_Data_Path);

disp("Select FBA Results Excel file");
FBA_Results_Path = uigetfile;
FBA_Results_Data = readtable(FBA_Results_Path);

Absolute_PPT = input("Would you like to plot Absolute Metabolites? (1/0) \n");
Relative_PPT = input("Would you like to plot Relative Metabolites? (1/0) \n");
Extra_Exch_Rxns_PPT = input("Would you like to plot Extra Exchange Metabolites? (1/0) \n");
All_Remaining_Rxns_PPT = input("Would you like to plot All Remaining Metabolites? (1/0) \n");

X_Axis_Type = input("Would you like the X-Axis to be Time(1) or OD(0) ? (1/0) \n");

Img_Save = input("Would you like to save the images? (1/0) \n");

%% Ind_Flux and Ind_Fluxr according to Results Excel Sheet %%
Ind_Flux = ind_flux + 3;
Ind_Fluxr = ind_fluxr + 3;

% Finding extra exch rxns %
ind_flux_AllExch = find(exchrxn ~= 0);
ind_flux_Exch = setdiff(ind_flux_AllExch,[ind_flux,ind_fluxr]);
Ind_Flux_Exch = ind_flux_Exch + 3;

% Finding all remaining rxns %
ind_flux_All = find(exchrxn == 0);
Ind_Flux_All = ind_flux_All + 3;

%% Flux Names %
FluxNames_Abs = FBA_Results_Data.Var1(Ind_Flux);
FluxNames_Rel = FBA_Results_Data.Var1(Ind_Fluxr);
FluxNames_Extra = FBA_Results_Data.Var1(Ind_Flux_Exch);
FluxNames_All = FBA_Results_Data.Var1(Ind_Flux_All);

%% Plots %

% X-Axis %
if X_Axis_Type == 0
    X_Input = OD;
    X_Output = FBA_Results_Data{2,2:end};
    X_Label_Text = "OD_{600}";
elseif X_Axis_Type == 1
    X_Input = Time;
    X_Output = FBA_Results_Data{1,2:end};
    X_Label_Text = "Time (h)";
end

% Y-Axis %
Y_Input_Abs = table2array(Final_Flux_Data_Abs);
Y_Input_Abs_LB = table2array(Final_Flux_Data_Abs_LB);
Y_Input_Abs_UB = table2array(Final_Flux_Data_Abs_UB);

Y_Input_Rel = table2array(Final_Flux_Data_Rel);
Y_Input_Rel_LB = table2array(Final_Flux_Data_Rel_LB);
Y_Input_Rel_UB = table2array(Final_Flux_Data_Rel_UB);

Y_Output_Abs = transpose(FBA_Results_Data{Ind_Flux,2:end});

Y_Output_Rel = transpose(FBA_Results_Data{Ind_Fluxr,2:end});

Y_Output_Exch = transpose(FBA_Results_Data{Ind_Flux_Exch,2:end});

Y_Output_All = transpose(FBA_Results_Data{Ind_Flux_All,2:end});

% Colors %
Color_Output = 'red';
Color_Input = 'blue';

%% Absolute Metabolites %%
if Absolute_PPT == 1
    
    mkdir PPT_Abs;
    addpath PPT_Abs;
    cd PPT_Abs; 

    import mlreportgen.ppt.*
    PPT_Abs = Presentation("FBA Model Results - Absolute.pptx");
    titleSlide = add(PPT_Abs,'Title Slide');
    replace(titleSlide,'Title', "FBA Model Results - Absolute");
    
    disp("Plotting Absolute Metabolites");
    for i = 1:length(Ind_Flux)
        f = figure('Visible','off');
        plot(X_Input,Y_Input_Abs(:,i),'LineStyle','-','Color',Color_Input);
        title(FluxNames_Abs(i));
        ylabel("Absolute Flux (mmol/gDCW.h)");
        xlabel(X_Label_Text);
        if Y_Input_Abs(:,i)==0
            ylim([-0.00000001 0.00000001])
        else
            ylim([-max(max(abs(Y_Input_Abs(:,i)),[],'all'),max(max(abs(Y_Input_Abs_LB(:,i)),[],'all'),max(abs(Y_Input_Abs_UB(:,i)),[],'all'))) max(max(abs(Y_Input_Abs(:,i)),[],'all'),max(max(abs(Y_Input_Abs_LB(:,i)),[],'all'),max(abs(Y_Input_Abs_UB(:,i)),[],'all')))]);
        end
        hold on
        plot(X_Input,Y_Input_Abs_LB(:,i),'LineStyle',':','Color','black');
        plot(X_Input,Y_Input_Abs_UB(:,i),'LineStyle',':','Color','black');
        plot(X_Output,Y_Output_Abs(:,i),'LineStyle','-','Color',Color_Output);
        legend("Input","Input LB","Input UB","Output",'Location','best');
        hold off
        saveas(gcf,FluxNames_Abs(i) + ".jpg");
        pictureSlide = add(PPT_Abs,'Title and Picture');
        replace(pictureSlide,'Title',FluxNames_Abs(i));
        replace(pictureSlide,'Picture',Picture(FluxNames_Abs(i) + ".jpg"));
        disp(i);
    end
    
    close(PPT_Abs);
    if Img_Save == 0
        delete *.jpg;
    end
    cd ../;
end



%% Relative Metabolites %%
if Relative_PPT == 1
    
    mkdir PPT_Rel;
    addpath PPT_Rel;
    cd PPT_Rel; 

    import mlreportgen.ppt.*
    PPT_Rel = Presentation("FBA Model Results - Relative.pptx");
    titleSlide = add(PPT_Rel,'Title Slide');
    replace(titleSlide,'Title', "FBA Model Results - Relative");
    
    disp("Plotting Relative Metabolites");
    for i = 1:length(Ind_Fluxr)
        f = figure('Visible','off');
        yyaxis left
        plot(X_Input,Y_Input_Rel(:,i),'LineStyle','-','Color',Color_Input);
        title(FluxNames_Rel(i));
        ylabel("Input Relative Flux");
        xlabel(X_Label_Text);
        ylim([-max(abs(Y_Input_Rel(:,i)),[],'all') max(abs(Y_Input_Rel(:,i)),[],'all')]);
        hold on
        plot(X_Input,Y_Input_Rel_LB(:,i),'LineStyle',':','Color','black');
        plot(X_Input,Y_Input_Rel_UB(:,i),'LineStyle',':','Color','black');
        yyaxis right
        plot(X_Output,Y_Output_Rel(:,i),'LineStyle','-','Color',Color_Output);
        ylabel("Output Absolute Flux (mmol/gDCW.h)");
        ylim([-(max(abs(Y_Output_Rel(:,i)),[],'all')+0.00001) (max(abs(Y_Output_Rel(:,i)),[],'all')+0.00001)]);
        legend("Input","Input LB","Input UB","Output",'Location','best');
        hold off
        saveas(gcf,FluxNames_Rel(i) + ".jpg");
        pictureSlide = add(PPT_Rel,'Title and Picture');
        replace(pictureSlide,'Title',FluxNames_Rel(i));
        replace(pictureSlide,'Picture',Picture(FluxNames_Rel(i) + ".jpg"));
        disp(i);
    end
    
    close(PPT_Rel);
    if Img_Save == 0
        delete *.jpg;
    end
    cd ../;
end



%% Exchange Metabolites %%
if Extra_Exch_Rxns_PPT == 1
    
    mkdir PPT_Exch;
    addpath PPT_Exch;
    cd PPT_Exch;

    import mlreportgen.ppt.*
    PPT_Exch = Presentation("FBA Model Results - Exchange.pptx");
    titleSlide = add(PPT_Exch,'Title Slide');
    replace(titleSlide,'Title', "FBA Model Results - Exchange");
    
    disp("Plotting Exchange Metabolites");
    for i = 1:length(Ind_Flux_Exch)
        f = figure('Visible','off');
        plot(X_Output,Y_Output_Exch(:,i),'LineStyle','-','Color',Color_Output);
        title(FluxNames_Extra(i));
        ylabel("Output Absolute Flux (mmol/gDCW.h)");
        xlabel(X_Label_Text);
        saveas(gcf,FluxNames_Extra(i) + ".jpg");
        pictureSlide = add(PPT_Exch,'Title and Picture');
        replace(pictureSlide,'Title',FluxNames_Extra(i));
        replace(pictureSlide,'Picture',Picture(FluxNames_Extra(i) + ".jpg"));
        disp(i);
    end

    close(PPT_Exch);
    if Img_Save == 0
        delete *.jpg;
    end
    cd ../;
end



%% All Metabolites %%
if All_Remaining_Rxns_PPT == 1
    
    mkdir PPT_All;
    addpath PPT_All;
    cd PPT_All 

    Max_Itr = 1000;
    c = 1;
    Q = floor(length(Ind_Flux_All)/Max_Itr);
    R = rem(length(Ind_Flux_All),Max_Itr);
    
    for i = 1:Q
        import mlreportgen.ppt.*
        PPT_All(i) = Presentation("FBA Model Results - All " + i + ".pptx");
        titleSlide = add(PPT_All,'Title Slide');
        replace(titleSlide,'Title', "FBA Model Results - All" + i);
    
        disp("Plotting All Remaining Metabolites" + i);
        for j = c:c+Max_Itr
            f = figure('Visible','off');
            plot(X_Output,Y_Output_All(:,j),'LineStyle','-','Color',Color_Output);
            title(FluxNames_All(j));
            ylabel("Output Absolute Flux (mmol/gDCW.h)");
            xlabel(X_Label_Text);
            saveas(gcf,FluxNames_All(j) + ".jpg");
            pictureSlide = add(PPT_All,'Title and Picture');
            replace(pictureSlide,'Title',FluxNames_All(j));
            replace(pictureSlide,'Picture',Picture(FluxNames_All(j) + ".jpg"));
            disp(j);
        end
        
        close(PPT_All(i));
        if Img_Save == 0
            delete *.jpg;
        end

        c = c + 1000;
    end

    import mlreportgen.ppt.*
    PPT_All(Q+1) = Presentation("FBA Model Results - All " + (Q+1) + ".pptx");
    titleSlide = add(PPT_All,'Title Slide');
    replace(titleSlide,'Title', "FBA Model Results - All" + (Q+1));
    
    disp("Plotting All Remaining Metabolites" + (Q+1));
    for j = ((Q*Max_Itr)+1):length(Ind_Flux_All)
        f = figure('Visible','off');
        plot(X_Output,Y_Output_All(:,j),'LineStyle','-','Color',Color_Output);
        title(FluxNames_All(j));
        ylabel("Output Absolute Flux (mmol/gDCW.h)");
        xlabel(X_Label_Text);
        saveas(gcf,FluxNames_All(j) + ".jpg");
        pictureSlide = add(PPT_All,'Title and Picture');
        replace(pictureSlide,'Title',FluxNames_All(j));
        replace(pictureSlide,'Picture',Picture(FluxNames_All(j) + ".jpg"));
        disp(j);
    end
    
    close(PPT_All(Q+1));
    if Img_Save == 0
        delete *.jpg;
    end
    cd ../;
end



