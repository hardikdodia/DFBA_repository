%% Estimate of fluxes

%% load data
% load PreliminaryDataFBA.mat

%% variable description
% model and modelpFBA  = genome-scale model of E. coli from BiGG database
% meas_flux = average estimates of absolute uptake/secretion rates
% fluxm_lb, fluxm_ub = upper/lower bounds on the absolute measurements of uptake/secretion rates
% ind_flux = indexes corresponding to the reactions for which absolute
% measruments of uptake/secretion rates are available
% meas_fluxr =  =  estimates of relative uptake/secretion rates
% ind_fluxr = indexes corresponding to the reactions for which relative
% measurements of uptake/secretion rates are available
% obj_tol = tolerance on the objective function for FVA analysis
% FVAcheckwt = if 1 FVA analysis is performed
% ODist = Optical Density 
% timesFBA = times (hour)
% pre_fluxr = vector of 0/1 indicating exchange fluxes for which corresponding metabolites were detected in the MS analysis of the supernatant

additions;

%% initialize parameters
obj_tol = 0.00;
FVAcheckwt = 0; 

modelpFBA.lb(find(exchrxn)) = -1000;
modelpFBA.ub(find(exchrxn)) = 1000;
model = modelpFBA;
% model.S

%% optimal solution 
% The following matlab function uses IBM ILOG CPLEX Optimization Studio. Please notice that the MATLAB connector from CPLEX is tested for the available versions of CPLEX at development time. 
% As such, each version of CPLEX is only guaranteed to be compatible with a subset of MATLAB releases.

[modelpFBAFVA] = Flux_fitting_dynFIAFBA_paper(modelpFBA, model, meas_flux, ind_flux, fluxm_lb, fluxm_ub, meas_fluxr, ind_fluxr, obj_tol, FVAcheckwt,ODist,timesFBA,pre_fluxr,FC,MFC,F_V,C_rel);
