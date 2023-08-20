model.rxnNames = strrep(model.rxnNames,char(58),"_"); % colon
model.rxnNames = strrep(model.rxnNames,char(59),"_"); % semi-colon
model.rxnNames = strrep(model.rxnNames,char(92),"_"); % back slash
model.rxnNames = strrep(model.rxnNames,char(47),"_"); % forward slash
model.rxnNames = strrep(model.rxnNames,char(39),"_"); % inverted comma
model.rxnNames = strrep(model.rxnNames,char(34),"_"); % double inverted comma
model.rxnNames = strrep(model.rxnNames,char(60),"_"); % less than
model.rxnNames = strrep(model.rxnNames,char(62),"_"); % greater than
model.rxnNames = strrep(model.rxnNames,char(124),"_"); % straight slash

writematrix(FBA_results,'DFBA_results.xlsx','Sheet','Flux_results_','Range','B4')
writematrix(model.rxnNames,'DFBA_results.xlsx','Sheet','Flux_results_','Range','A4')

headings = {'Time','OD'}';

writecell(headings,'DFBA_results.xlsx','Sheet','Flux_results_')

OD_time = [timesFBA ; ODist];

writematrix(OD_time,'DFBA_results.xlsx','Sheet','Flux_results_','Range','B1')