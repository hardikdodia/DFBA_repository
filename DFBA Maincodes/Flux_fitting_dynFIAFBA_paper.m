function [solutionsWT ] = Flux_fitting_dynFIAFBA_paper(modelWT, model,fluxm,ind_fluxm,fluxm_lb,fluxm_ub,fluxr,ind_fluxr,obj_tol, FVAcheckwt,ODist,timesFBA,pre_fluxr,FC,MFC,F_V,C_rel)

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% changeCobraSolver('glpk');
% set constraints
N_REAC = size(modelWT.S,2);
N_MET = size(modelWT.S,1);

N_MEAS = size(fluxm,2);
N_MEAS_REL = size(fluxr,2);


Nm = size(fluxm,1) % number of mutants

%% round inputs to reduce complexity
fluxr = round(fluxr*1e4)/1e4;
%% Aeq
Aeqs = sparse([],[],[],N_MET+N_MEAS,2*N_REAC + 3*N_MEAS, 2*nnz(modelWT.S~=0) + N_MEAS*3);
Aeqs(1:N_MET,1:N_REAC) =  modelWT.S;
% measured fluxes
Aeqs(N_MET+1:N_MET+N_MEAS,2*N_REAC+1:2*N_REAC+N_MEAS) = -eye(N_MEAS, N_MEAS);   %vm
Aeqs(N_MET+1:N_MET+N_MEAS,2*N_REAC+N_MEAS+1:2*N_REAC+2*N_MEAS) = -eye(N_MEAS, N_MEAS);   %e
Aeqs(N_MET+1:N_MET+N_MEAS,ind_fluxm) = eye(N_MEAS, N_MEAS);                 %vwt

beqs = zeros(size(Aeqs,1),1);
% repeating the constraints for each time
Aeq = sparse([]);
for i=1:Nm
    Aeq = [Aeq sparse(size(Aeq,1),size(Aeqs,2));sparse(size(Aeqs,1),size(Aeq,2)) Aeqs];
end
beq = repmat(beqs,Nm,1);
% add proportionality constants (C) , and dummy variables for the
% (absolute) distance (D and |D|)
temp = sparse(size(Aeq,1),N_MEAS_REL+N_MEAS_REL*Nm*2)
Aeq = [Aeq temp];

for i=1:Nm
    for j=1:length(ind_fluxr)
        Aeq(size(Aeq,1) + 1,ind_fluxr(j) + (i-1)*size(Aeqs,2)) = 1;
        Aeq(size(Aeq,1),Nm*size(Aeqs,2)+j) = -fluxr(i,j);
        Aeq(size(Aeq,1),Nm*size(Aeqs,2)+N_MEAS_REL+(i-1)*N_MEAS_REL+j) = -1;
    end
end
beq = [beq; zeros(Nm*N_MEAS_REL,1)];
%% Aineq

Aineqs = sparse(2*N_REAC+2*N_MEAS,size(Aeqs,2));
% absolute vwt
Aineqs(1:N_REAC,1:N_REAC) =  eye(N_REAC,N_REAC);
Aineqs(1:N_REAC,N_REAC+1:2*N_REAC) =  -eye(N_REAC,N_REAC);

Aineqs(N_REAC+1:2*N_REAC,1:N_REAC) =  -eye(N_REAC,N_REAC);
Aineqs(N_REAC+1:2*N_REAC,N_REAC+1:2*N_REAC) =  -eye(N_REAC,N_REAC);

% absolute distance value with measured values |e|
Aineqs(2*N_REAC+1:2*N_REAC+N_MEAS,2*N_REAC+N_MEAS+1:2*N_REAC+2*N_MEAS) =  eye(N_MEAS,N_MEAS);
Aineqs(2*N_REAC+1:2*N_REAC+N_MEAS,2*N_REAC+2*N_MEAS+1:2*N_REAC+3*N_MEAS) =  -eye(N_MEAS,N_MEAS);

Aineqs(2*N_REAC+N_MEAS+1:2*N_REAC+2*N_MEAS,2*N_REAC+N_MEAS+1:2*N_REAC+2*N_MEAS) =  -eye(N_MEAS,N_MEAS);
Aineqs(2*N_REAC+N_MEAS+1:2*N_REAC+2*N_MEAS,2*N_REAC+2*N_MEAS+1:2*N_REAC+3*N_MEAS) =  -eye(N_MEAS,N_MEAS);

bineqs = zeros(size(Aineqs,1),1);
% repeating the constraints for each time
Aineq = sparse([]);
for i=1:Nm
    Aineq = [Aineq sparse(size(Aineq,1),size(Aineqs,2));sparse(size(Aineqs,1),size(Aineq,2)) Aineqs];
end

% FC AND MFC

FC(FC<0) = 1;
MFC(MFC<0) = 1;


% delta t

del_t = diff([0 timesFBA]);


% adding constraints for absolute D
Aineq = [Aineq sparse(size(Aineq,1),size(Aeq,2)-size(Aineq,2))];
temp = eye(Nm*N_MEAS_REL,Nm*N_MEAS_REL);
Aineq = [Aineq; [sparse(N_MEAS_REL*Nm,size(Aeq,2)- 2*N_MEAS_REL*Nm) temp -temp]];
Aineq = [Aineq; [sparse(N_MEAS_REL*Nm,size(Aeq,2)- 2*N_MEAS_REL*Nm) -temp -temp]];
bineq = zeros(size(Aineq,1),1);
% adding temporal constraints on the maximum amount of extracellular metabolite
temp = [sparse(N_MEAS_REL*Nm+N_MEAS_REL, size(Aineq,2))];
count =1;
for i=1:N_MEAS_REL
    temp(count,(0:Nm-1)*size(Aeqs,2)+ind_fluxr(i)) = ODist.*0.370.*diff([0 timesFBA]);
    bineq(length(bineq)+1) = 10000;
    count = count + 1;
    for j=1:Nm
%         temp(count,(0:j-1)*size(Aeqs,2)+ind_fluxr(i)) = -ODist(1:j).*0.370.*diff([0 timesFBA(1:j)]);
        temp(count,(j-1)*size(Aeqs,2)+ind_fluxr(i)) = -ODist(j).*0.370.*del_t(j);
        bineq(length(bineq)+1) = (0.05*pre_fluxr(i,j)*(FC(i)))+(F_V(j)*MFC(i)*0.05*del_t(j));
        count = count + 1;
    end
    i
end

Aineq = [Aineq; temp];

%% prepare Aineq for the next 3 constraints
Aineq = [Aineq; sparse(3,size(Aineq,2))];
bineq = [bineq; zeros(3,1)];
%% constraints
LB = [];
UB = [];

fluxm_lb = (fluxm)-abs(0.1*fluxm);
fluxm_ub = (fluxm)+abs(0.1*fluxm);


count = 0;

for i = 1:Nm
%     LB = [LB; [modelWT.lb ;  0*ones(N_REAC,1); fluxm_lb(i,:)'; -1e4*ones(2*N_MEAS,1)]];
%     UB = [UB; [modelWT.ub ;  1e3*ones(N_REAC,1); fluxm_ub(i,:)'; 1e4*ones(2*N_MEAS,1)]];
    for j=1:length(ind_fluxr)          % Modification made here to be more stringent in ub and lb for fluxes
        if fluxr(i,j)>0
            modelWT.lb(ind_fluxr(j))=0;
            modelWT.ub(ind_fluxr(j))=1000; 
            count = count+1;
        else 
             if fluxr(i,j)<0
                modelWT.ub(ind_fluxr(j))=0;
                modelWT.lb(ind_fluxr(j))=-1000;
             end
        end
    end
%     modelWT.lb(ind_fluxm) = fluxm_lb(i,:);
%     modelWT.ub(ind_fluxm) = fluxm_ub(i,:);
    LB = [LB; [modelWT.lb ;  0*ones(N_REAC,1); fluxm_lb(i,:)'; -1e4*ones(2*N_MEAS,1)]];
    UB = [UB; [modelWT.ub ;  1e3*ones(N_REAC,1); fluxm_ub(i,:)'; 1e4*ones(2*N_MEAS,1)]];

end

c_lb_size = size(LB,1) ;
% imposing constraints on C
% C_LB = zeros(1,N_MEAS_REL);
% for i = 1:N_MEAS_REL
%     if nnz((cumsum(ODist.*0.370.*diff([0 timesFBA]).*fluxr(:,i)'.*fluxm(:,1)'))>0)>0 & pre_fluxr(i)==0
%         [frmax imax] = max((cumsum(ODist.*0.370.*diff([0 timesFBA]).*fluxr(:,i)'.*fluxm(:,1)')));
%         imin = min(setdiff(find((cumsum(ODist.*0.370.*diff([0 timesFBA]).*fluxr(:,i)'.*fluxm(:,1)'))>0),1));
%         if ~isempty(imin)
%             C_LB(i) = 0.05/ sum(ODist(imin:imax).*0.370.*diff([0 timesFBA(imin:imax)]).*fluxr(imin:imax,i)'.*fluxm(imin:imax,1)');
%         end
%     end
% end

C_LB = zeros(1,N_MEAS_REL);
for i = 1:N_MEAS_REL
    if nnz((cumsum(ODist.*0.370.*diff([0 timesFBA]).*fluxr(:,i)'))>0)>0 
        [frmax imax] = max((cumsum(ODist.*0.370.*diff([0 timesFBA]).*fluxr(:,i)')));
        imin = min(setdiff(find((cumsum(ODist.*0.370.*diff([0 timesFBA]).*fluxr(:,i)'))>0),1));
        if ~isempty(imin)
            C_LB(i) = (0.05)/ (sum(ODist(imin:imax).*0.370.*diff([0 timesFBA(imin:imax)]).*fluxr(imin:imax,i)'));
        end
    
    else
        C_LB(i) = (0.05)/(C_rel(1,i)+C_rel(end,i));
    end
end

% LB = [LB; 0*ones(N_MEAS_REL,1); -1e6*ones(2*N_MEAS_REL*Nm,1)];
% LB = [LB; C_LB'; -1e6*ones(2*N_MEAS_REL*Nm,1)];
% UB = [UB; 1e6*ones(N_MEAS_REL,1); 1e6*ones(2*N_MEAS_REL*Nm,1)];


LB = [LB; C_LB'; -1e6*ones(2*N_MEAS_REL*Nm,1)];
UB = [UB; 1e6*ones(N_MEAS_REL,1); 1e6*ones(2*N_MEAS_REL*Nm,1)];



VARTYPE = repmat('C',1,size(Aeq,2));



% fix growth rate
for i=1:Nm

    LB((i-1)*size(Aeqs,2)+22) = fluxm(i,1);
    UB((i-1)*size(Aeqs,2)+22) = fluxm(i,1);
%     LB((i-1)*size(Aeqs,2)+2759) = fluxm_lb(i,2);
%     UB((i-1)*size(Aeqs,2)+2759) = fluxm_ub(i,2);

end


% minimize difference between measured and inferred flux in the wt
Obj=zeros(1,size(Aeq,2));
for i=1:Nm
    Obj((i-1)*size(Aeqs,2)+2*N_REAC+2*N_MEAS+1:(i-1)*size(Aeqs,2)+2*N_REAC+2*N_MEAS+length(ind_fluxm)) = 1;
end
%% solve 1
options = cplexoptimset ('cplex')
options.simplex.tolerances.feasibility=0;
options.feasopt.tolerance=0;
options.emphasis.numerical = 1;
options.lpmethod = 2;
options.diagnostics = 'on';
options.display = 'on';
% options.preprocessing.presolve = 0;
%% minimize L1 norm between predicted and mesured fluxes

tic
[xmax, fval, exitflag, output] = cplexlp (Obj', Aineq, bineq, Aeq , beq, LB, UB, [], options);

toc

solutionsWT.obj1 = fval;
solutionsWT.c1 = xmax(c_lb_size+1:c_lb_size+N_MEAS_REL);
solutionsWT.xmax1 = xmax;

%% minimize distance with relative fluxes
% constrain measured fluxes
Aineq(size(Aineq,1)-2,find(Obj)) = 1;
bineq(length(bineq)-2) = fval +(fval*obj_tol*0);

% minimize D
Obj = 0.*Obj;
Obj((Nm*size(Aeqs,2)+N_MEAS_REL+Nm*N_MEAS_REL+1:Nm*size(Aeqs,2)+N_MEAS_REL+Nm*N_MEAS_REL+Nm*N_MEAS_REL)) =1;

tic
[xmax2, fval2, exitflag2, output2] = cplexlp (Obj', Aineq, bineq, Aeq , beq, LB, UB, [], options);
toc

solutionsWT.xmax2 = xmax2(Obj((Nm*size(Aeqs,2)+N_MEAS_REL+Nm*N_MEAS_REL+1:Nm*size(Aeqs,2)+N_MEAS_REL+Nm*N_MEAS_REL+Nm*N_MEAS_REL)));

solutionsWT.obj2 = fval2;

solutionsWT.c2 = xmax2(c_lb_size+1:c_lb_size+N_MEAS_REL);

c2_temp = xmax2(c_lb_size+1:c_lb_size+N_MEAS_REL);

temp_fluxdr = [];
temp_fluxd = [];
for i=1:Nm
    temp_fluxd = [temp_fluxd xmax2((i-1)*size(Aeqs,2)+(1:N_REAC))];
end
[model.rxnNames(ind_fluxr) num2cell(temp_fluxd(ind_fluxr,:))]
[model.rxnNames(:) num2cell(temp_fluxd)]
%% minimize total sum of absolute fluxes
% constrain measured fluxes
Aineq(size(Aineq,1)-1,find(Obj)) = 1;
bineq(length(bineq)-1) = sum(xmax2(find(Obj))) +(sum(xmax2(find(Obj)))*1e-4);

% minimize C
Obj = 0.*Obj;
for i=1:Nm
    Obj((i-1)*size(Aeqs,2)+N_REAC+1:(i-1)*size(Aeqs,2)+N_REAC+N_REAC) =1;
end

options.lpmethod = 1;
% options.read.scale = -1;
tic
[xmax3, fval3, exitflag3, output3] = cplexlp (Obj', Aineq, bineq, Aeq , beq, LB, UB, [], options);
toc

solutionsWT.obj3 = fval3;

temp_fluxdr = [];
temp_fluxd = [];
for i=1:Nm
    temp_fluxd = [temp_fluxd xmax3((i-1)*size(Aeqs,2)+(1:N_REAC))];
end
[model.rxnNames(ind_fluxr) num2cell(temp_fluxd(ind_fluxr,:))];
[model.rxnNames(:) num2cell(temp_fluxd)];
solutionsWT.minabs = temp_fluxd;


%% perform FVA
% % constrain measured fluxes
if FVAcheckwt==1

Aineq(size(Aineq,1),find(Obj)) = 1;
bineq(length(bineq)) = sum(xmax3(find(Obj))) +(sum(xmax3(find(Obj)))*0.001);

vwt_max = zeros(N_REAC,Nm);
vwt_min = zeros(N_REAC,Nm);

cplx = Cplex()
cplx.Model.obj = Obj';
cplx.Model.lb = xmax3;
cplx.Model.ub = xmax3;
cplx.Model.A = [Aeq; Aineq];
cplx.Model.rhs = [beq; bineq];
cplx.Model.lhs = [beq; -Inf.*ones(size(bineq))];

cplx.Param.lpmethod.Cur =1;

out3 = solve(cplx)

count = 1;

%% Hard constraint growth rate for FVA

fluxm_lb(:,1) = fluxm(:,1);
fluxm_ub(:,1) = fluxm(:,1);

%% FVA for each reaction across time points

for i=1:N_REAC
    tic
    if i ==22 || i ==23|| i ==40
        [~,count] = ismember(i,ind_fluxm);
        vwt_max(i,:) = fluxm(:,count);
        vwt_min(i,:) = fluxm(:,count);
    else


    Objwt=sparse(1,size(Aeq,2));
    Objwt((0:Nm-1)*size(Aeqs,2)+i) = 1;

    if ismember(i,ind_fluxm) 
        [~,count] = ismember(i,ind_fluxm) ;   
        LB((0:Nm-1)*size(Aeqs,2)+i,1) = fluxm_lb(:,count);
        UB((0:Nm-1)*size(Aeqs,2)+i,1) = fluxm_ub(:,count);
    end

    options.emphasis.numerical = 4;
    options.simplex.tolerances.markowitz = 0.9999;
    options.simplex.tolerances.feasibility = 1e-9;
    options.simplex.tolerances.optimality = 1e-9;
    [xmax_max, fval, exitflag, output] = cplexlp (-Objwt', Aineq, bineq, Aeq , beq, LB, UB, xmax3, options);
    
    vwt_max(i,:) = xmax_max((0:Nm-1)*size(Aeqs,2)+i);
    [xmax_min, fval, exitflag, output] = cplexlp (Objwt', Aineq, bineq, Aeq , beq, LB, UB, xmax3, options);
    vwt_min(i,:) = xmax_min((0:Nm-1)*size(Aeqs,2)+i);
    
    i

    end
    toc
end


solutionsWT.FVA_min = vwt_min;
solutionsWT.FVA_max = vwt_max;
end

%% results
FBA_results = solutionsWT.minabs;
results;
save("model_results_param_modified_c_final.mat","solutionsWT");

end