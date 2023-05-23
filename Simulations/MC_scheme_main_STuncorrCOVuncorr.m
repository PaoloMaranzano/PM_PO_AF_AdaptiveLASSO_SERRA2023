%% Paper: Adaptive LASSO estimation for functional hidden dynamic geostatistical models
%% Journal: Stochastic Environmental Research and Risk Assessment
%% Date: July 2022


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%% Monte Carlo simulation of penalized   %%%%%%%%%% %%
%% %%%%%%%%%%  f-HDGM using adaptive LASSO penalty  %%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% %%%%%%%%%% Main program %%%%%%%%%% %%
clear
clc


%%% Upload real data from ARPA Lombardia network
addpath(genpath('../../../VisitingLUH2021'));
% addpath(genpath('../../../SPASTA2021/Code/DSTEM_software'));
% addpath(genpath('../../../SPASTA2021/Code/Matlab_auxfuns'));
addpath(genpath('../../../DSTEM_Code_General/DSTEM_software'));
addpath(genpath('../../../DSTEM_Code_General/Matlab_auxfuns'));
addpath(genpath('../../../SPASTA2021/Data'));

if ~exist('Ground')
	paper = "CSDA2021";
    run('Reshape_to_HDGM.m');
end
clearvars -except crossval_step data data_fhdgm Ground log_transform poll standardization


%%%%% Settings
% Estimate varcov of the full model
varcov_estim_flag = 1;
% Estimate the log-likelihood of the full model
logL_estim_flag = 1;
% Correlated covariates
cov_corr_flag = 0;
% Spatio-temporal correlation
st_corr_flag = 0;
% Standardize the data before estimation
standardization_flag = 0;
% Use de-standardized coefficients and varcov for PMLE
beta_dest_flag = 0;
% Penalize the intercept or not
beta0_pen_flag = 0;
% Exact or approximated VarCov computation
ExactVarcov = 1;

%%% File number (distributed computing)
File_number = 1;

%%% Number of Monte Carlo simulation for each file
MC_reps = 1;
%%% Number of cumulated Monte Carlo simulation for the curret simulation scheme
MC_cum = 1;

% Save the final output
save_output_flag = 0;


%%%%% Run simulation code
run('MC_scheme_alg_penlik.m');