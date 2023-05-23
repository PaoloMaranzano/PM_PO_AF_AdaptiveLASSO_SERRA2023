%% Paper: Adaptive LASSO estimation for functional hidden dynamic geostatistical models
%% Journal: Stochastic Environmental Research and Risk Assessment
%% Date: July 2022


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%% RMSE and MAE reference values %%%%%%%%%% %%
%% %%%%%%%%%%           Simulated           %%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% %%%%%%%%%% Main program %%%%%%%%%% %%
clear
clc

%%%%% Settings working directory and clear environment
if 1
    cd 'H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim'
    addpath(genpath('../../../VisitingLUH2021'));
    addpath(genpath('../../../DSTEM_Code_General/DSTEM_software'));
    addpath(genpath('../../../DSTEM_Code_General/Matlab_auxfuns'));
    addpath(genpath('../../../SPASTA2021/Data'));
end

if ~exist('Ground')
	paper = "CSDA2021";
    run('Reshape_to_HDGM.m');
end
clearvars -except crossval_step data data_fhdgm Ground log_transform poll standardization



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      General inputs for simulations     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DSTEM_str_ground = Ground;
DSTEM_str_ground.poll = {'NO2'};
n_sites = {15};
n_covs = 3;
nvars = n_covs + 1;
datestamp_begin = '01-01-2017 00:00';
datestamp_end = '31-12-2017 23:00';
model_type = 'f-HDGM';

%%% Setup when B-spline are used
fda_setup.spline_type = 'Bspline';      % Spline type
fda_setup.spline_order = 3;
fda_setup.knots_number = 5;
fda_setup.spline_range = [0 24];        % Spline basis domain
fda_setup.spline_knots = linspace(fda_setup.spline_range(1),fda_setup.spline_range(2),fda_setup.knots_number);
fda_setup.n_basis = (fda_setup.knots_number + 2*fda_setup.spline_order - (fda_setup.spline_order+1))*nvars;
n_basis = fda_setup.n_basis;

%%% Fixed-effects and random-effects components
ones_numb = 4;
DSTEM_str_sim_setup.beta = [[1 1 1 1 0 0 0],repmat(repelem([1 0],[ones_numb n_basis/nvars-ones_numb]),1,n_covs)]';
beta_true = DSTEM_str_sim_setup.beta;
DSTEM_str_sim_setup.nan_rate = [];
DSTEM_str_sim_setup.nan_pattern_par = [];
DSTEM_str_sim_setup.sigma_eps = repelem(0,n_basis/nvars)';





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Simulating reference values      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsims = 100;

%% %%%%% Setting I: uncorrelated covariates and ST uncorrelated observations %%%%% 
%%% Setup
cov_corr_flag = 0;
st_corr_flag = 0;
clear i Decomp Sims_unex
% Correlated/uncorrelated covariates
if cov_corr_flag == 1
    cov_msg = 'COVcorr';
	X_varcov1 = [
    1     0.9  0.70 ;
    0.90    1  0.50 ;
    0.70  0.5   1
    ];
else
    cov_msg = 'COVuncorr';
	X_varcov1 = [
    1    0   0 ;
    0    1   0 ;
    0    0   1
    ];
end
X_varcov = repmat({X_varcov1},1,length(DSTEM_str_ground.poll));
X_names = {'X1','X2','X3'};
% Spatio-temporal correlation
if st_corr_flag == 1
    st_msg = 'STcorr';
    DSTEM_str_sim_setup.theta_z = repelem(km2deg(50),n_basis/nvars);
    DSTEM_str_sim_setup.G = diag(repelem(0.85,n_basis/nvars));
    DSTEM_str_sim_setup.v_z = diag(repelem(1,n_basis/nvars));
else
    st_msg = 'STuncorr';
    DSTEM_str_sim_setup.theta_z = repelem(km2deg(0.00001),n_basis/nvars);
    DSTEM_str_sim_setup.G = diag(repelem(0,n_basis/nvars));
    DSTEM_str_sim_setup.v_z = diag(repelem(0,n_basis/nvars));
end
%%% For loop
for i = 1:nsims
    %%%%% Simulate
    [DSTEM_obj_sim,obj_stem_par,DSTEM_obj_sim_str,sim_ground] = DSTEM_fHDGM_sim(...
        DSTEM_str_ground,DSTEM_str_sim_setup,...
        n_covs,n_sites,X_names,X_varcov,datestamp_begin, datestamp_end, model_type,0,...
        fda_setup);
    %%%%% Decomposing variance and error
    Sims_unex{i} = DSTEM_fHDGM_VarDecomp('ErrMat.csv','REMat.csv',...
        'YFEMat.csv','YFEREMat.csv',DSTEM_obj_sim);
end
%%% Save output
Decomp = vertcat(Sims_unex{:});
writetable(Decomp,'SimRefVals_Sett1.xlsx');
min_MSE = mean(Decomp.sigma2_eps_sim); std_min_MSE = std(Decomp.sigma2_eps_sim)/sqrt(nsims);
min_RMSE = mean(sqrt(Decomp.sigma2_eps_sim)); std_min_RMSE = std(sqrt(Decomp.sigma2_eps_sim))/sqrt(nsims);
min_MAE = mean(Decomp.MAE_eps_sim); std_min_MAE = std(Decomp.MAE_eps_sim)/sqrt(nsims);
max_MSE = mean(Decomp.MSE_eps_mu1); std_max_MSE = std(Decomp.MSE_eps_mu1)/sqrt(nsims);
max_RMSE = mean(Decomp.RMSE_eps_mu1); std_max_RMSE = std(Decomp.RMSE_eps_mu1)/sqrt(nsims);
max_MAE = mean(Decomp.MAE_eps_mu1); std_max_MAE = std(Decomp.MAE_eps_mu1)/sqrt(nsims);
SimRefVals_Sett1 = [
    min_MSE, min_RMSE, min_MAE ;
    std_min_MSE, std_min_RMSE, std_min_MAE ;
    max_MSE, max_RMSE, max_MAE ;
    std_max_MSE, std_max_RMSE, std_max_MAE
    ];


%% %%%%% Setting II: uncorrelated covariates and ST correlated observations %%%%% 
%%% Setup
cov_corr_flag = 0;
st_corr_flag = 1;
clear i Decomp Sims_unex
% Correlated/uncorrelated covariates
if cov_corr_flag == 1
    cov_msg = 'COVcorr';
	X_varcov1 = [
    1     0.9  0.70 ;
    0.90    1  0.50 ;
    0.70  0.5   1
    ];
else
    cov_msg = 'COVuncorr';
	X_varcov1 = [
    1    0   0 ;
    0    1   0 ;
    0    0   1
    ];
end
X_varcov = repmat({X_varcov1},1,length(DSTEM_str_ground.poll));
X_names = {'X1','X2','X3'};
% Spatio-temporal correlation
if st_corr_flag == 1
    st_msg = 'STcorr';
    DSTEM_str_sim_setup.theta_z = repelem(km2deg(50),n_basis/nvars);
    DSTEM_str_sim_setup.G = diag(repelem(0.85,n_basis/nvars));
    DSTEM_str_sim_setup.v_z = diag(repelem(1,n_basis/nvars));
else
    st_msg = 'STuncorr';
    DSTEM_str_sim_setup.theta_z = repelem(km2deg(0.00001),n_basis/nvars);
    DSTEM_str_sim_setup.G = diag(repelem(0,n_basis/nvars));
    DSTEM_str_sim_setup.v_z = diag(repelem(0,n_basis/nvars));
end
%%% For loop
for i = 1:nsims
    %%%%% Simulate
    [DSTEM_obj_sim,obj_stem_par,DSTEM_obj_sim_str,sim_ground] = DSTEM_fHDGM_sim(...
        DSTEM_str_ground,DSTEM_str_sim_setup,...
        n_covs,n_sites,X_names,X_varcov,datestamp_begin, datestamp_end, model_type,0,...
        fda_setup);
    %%%%% Decomposing variance and error
    Sims_unex{i} = DSTEM_fHDGM_VarDecomp('ErrMat.csv','REMat.csv',...
        'YFEMat.csv','YFEREMat.csv',DSTEM_obj_sim);
end
%%% Save output
Decomp = vertcat(Sims_unex{:});
writetable(Decomp,'SimRefVals_Sett2.xlsx');
min_MSE = mean(Decomp.sigma2_eps_sim); std_min_MSE = std(Decomp.sigma2_eps_sim)/sqrt(nsims);
min_RMSE = mean(sqrt(Decomp.sigma2_eps_sim)); std_min_RMSE = std(sqrt(Decomp.sigma2_eps_sim))/sqrt(nsims);
min_MAE = mean(Decomp.MAE_eps_sim); std_min_MAE = std(Decomp.MAE_eps_sim)/sqrt(nsims);
max_MSE = mean(Decomp.MSE_eps_mu1); std_max_MSE = std(Decomp.MSE_eps_mu1)/sqrt(nsims);
max_RMSE = mean(Decomp.RMSE_eps_mu1); std_max_RMSE = std(Decomp.RMSE_eps_mu1)/sqrt(nsims);
max_MAE = mean(Decomp.MAE_eps_mu1); std_max_MAE = std(Decomp.MAE_eps_mu1)/sqrt(nsims);
SimRefVals_Sett2 = [
    min_MSE, min_RMSE, min_MAE ;
    std_min_MSE, std_min_RMSE, std_min_MAE ;
    max_MSE, max_RMSE, max_MAE ;
    std_max_MSE, std_max_RMSE, std_max_MAE
    ];


%% %%%%% Setting III: correlated covariates and ST correlated observations %%%%% 
%%% Setup
cov_corr_flag = 1;
st_corr_flag = 1;
clear i Decomp Sims_unex
% Correlated/uncorrelated covariates
if cov_corr_flag == 1
    cov_msg = 'COVcorr';
	X_varcov1 = [
    1     0.9  0.70 ;
    0.90    1  0.50 ;
    0.70  0.5   1
    ];
else
    cov_msg = 'COVuncorr';
	X_varcov1 = [
    1    0   0 ;
    0    1   0 ;
    0    0   1
    ];
end
X_varcov = repmat({X_varcov1},1,length(DSTEM_str_ground.poll));
X_names = {'X1','X2','X3'};
% Spatio-temporal correlation
if st_corr_flag == 1
    st_msg = 'STcorr';
    DSTEM_str_sim_setup.theta_z = repelem(km2deg(50),n_basis/nvars);
    DSTEM_str_sim_setup.G = diag(repelem(0.85,n_basis/nvars));
    DSTEM_str_sim_setup.v_z = diag(repelem(1,n_basis/nvars));
else
    st_msg = 'STuncorr';
    DSTEM_str_sim_setup.theta_z = repelem(km2deg(0.00001),n_basis/nvars);
    DSTEM_str_sim_setup.G = diag(repelem(0,n_basis/nvars));
    DSTEM_str_sim_setup.v_z = diag(repelem(0,n_basis/nvars));
end
%%% For loop
for i = 1:nsims
    %%%%% Simulate
    [DSTEM_obj_sim,obj_stem_par,DSTEM_obj_sim_str,sim_ground] = DSTEM_fHDGM_sim(...
        DSTEM_str_ground,DSTEM_str_sim_setup,...
        n_covs,n_sites,X_names,X_varcov,datestamp_begin, datestamp_end, model_type,0,...
        fda_setup);
    %%%%% Decomposing variance and error
    Sims_unex{i} = DSTEM_fHDGM_VarDecomp('ErrMat.csv','REMat.csv',...
        'YFEMat.csv','YFEREMat.csv',DSTEM_obj_sim);
end
%%% Save output
Decomp = vertcat(Sims_unex{:});
writetable(Decomp,'SimRefVals_Sett3.xlsx');
min_MSE = mean(Decomp.sigma2_eps_sim); std_min_MSE = std(Decomp.sigma2_eps_sim)/sqrt(nsims);
min_RMSE = mean(sqrt(Decomp.sigma2_eps_sim)); std_min_RMSE = std(sqrt(Decomp.sigma2_eps_sim))/sqrt(nsims);
min_MAE = mean(Decomp.MAE_eps_sim); std_min_MAE = std(Decomp.MAE_eps_sim)/sqrt(nsims);
max_MSE = mean(Decomp.MSE_eps_mu1); std_max_MSE = std(Decomp.MSE_eps_mu1)/sqrt(nsims);
max_RMSE = mean(Decomp.RMSE_eps_mu1); std_max_RMSE = std(Decomp.RMSE_eps_mu1)/sqrt(nsims);
max_MAE = mean(Decomp.MAE_eps_mu1); std_max_MAE = std(Decomp.MAE_eps_mu1)/sqrt(nsims);
SimRefVals_Sett3 = [
    min_MSE, min_RMSE, min_MAE ;
    std_min_MSE, std_min_RMSE, std_min_MAE ;
    max_MSE, max_RMSE, max_MAE ;
    std_max_MSE, std_max_RMSE, std_max_MAE
    ];


%% %%%%% Combining outputs %%%%% 
RefVals = [SimRefVals_Sett1 ; SimRefVals_Sett2 ; SimRefVals_Sett3];
v1 = ["Setting I" , "Setting II" , "Setting III"];
c1 = categorical([repelem(v1,4)]');
v2 = ["min","max","min","max","min","max"];
c2 = categorical(repelem(v2,2)');
v3 = ["avg","sem","avg","sem","avg","sem"];
c3 = categorical([v3 , v3]');
RefVals = [array2table(c1), array2table(c2), array2table(c3), array2table(RefVals)];
RefVals.Properties.VariableNames = {'Setting','min/max','stat','MSE','RMSE','MAE'};
% RefVals = [array2table(c1), array2table(RefVals)];
% RefVals.Properties.VariableNames = {'Setting','MSE','RMSE','MAE'};
RefVals
writetable(RefVals,'SimRefVals.xlsx');


