%% Paper: Adaptive LASSO estimation for functional hidden dynamic geostatistical models
%% Journal: Stochastic Environmental Research and Risk Assessment
%% Date: March 2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Simulations results analysis %% %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Settings working directory and clear environment
if 1
    cd 'H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim'
    addpath(genpath('../../../VisitingLUH2021'));
    addpath(genpath('../../../DSTEM_Code_General/DSTEM_software'));
    addpath(genpath('../../../DSTEM_Code_General/Matlab_auxfuns'));
    clear
end

sim_scheme = {"STcorrCOVcorr","STcorrCOVuncorr","STuncorrCOVuncorr"};

for s = 1:length(sim_scheme)
    %%%%% Upload functional setup
    load('fda_setup.mat')

    %%%%% Upload simulation results
    if ismember(sim_scheme{s},'STuncorrCOVuncorr')
        load('Results_simulations\r1_sims\MCrep_1_STuncorrCOVuncorr_std0_betadest0_beta0unpen_1_500.mat')
        sup_text_tit = 'Setting I: spatio-temporal uncorrelation and uncorrelated covariates';
        save_plot_path = 'STuncorrCOVuncorr/';
        save_plot_name = 'STuncorrCOVuncorr';
    end
    
    if ismember(sim_scheme{s},'STcorrCOVuncorr')
        load('Results_simulations\r1_sims\MCrep_1_STcorrCOVuncorr_std0_betadest0_beta0unpen_1_500.mat')
        sup_text_tit = 'Setting II: spatio-temporal correlation and uncorrelated covariates';
        save_plot_path = 'STcorrCOVuncorr/';
        save_plot_name = 'STcorrCOVuncorr';
    end
    
    if ismember(sim_scheme{s},'STcorrCOVcorr')
        load('Results_simulations\r1_sims\MCrep_1_STcorrCOVcorr_std0_betadest0_beta0unpen_1_500.mat')
        sup_text_tit = 'Setting III: spatio-temporal correlation and correlated covariates';
        save_plot_path = 'STcorrCOVcorr/';
        save_plot_name = 'STcorrCOVcorr';
    end
       
    %% %%%%% Step 2. Computing summary statistics
    [CV_perf_metrics] = DSTEM_HDGM_CV_PenLik_SimsOptRes(MC_CV_perf_met,...
        beta_true,{'Intercept','X_1','X_2','X_3'},1,1,save_plot_path);
    
    %% %%%%% Step 3. Plot results
    %%% Plots A, B, C, D and F
    DSTEM_HDGM_CV_PenLik_OptRes_Plot(...
        CV_perf_metrics,beta_true,fda_setup,input_fda,...
        {'Intercept','X1','X2','X3'},...
        {'A','B','C','D','F'},...
        1,0,0,...
        [],[],...
        sup_text_tit,...
        1,save_plot_path,save_plot_name);
    {'A','B','C','D','F'},...
    %%% Plot E (Boxplot)
    crits = {'beta_min_MAE','beta_1se_MAE','beta_min_RMSE','beta_1se_RMSE'};
    box_tit = {
        'Distribution of the estimated coefficients at \lambda=\lambda^{*}_{MAE}',...
        'Distribution of the estimated coefficients at \lambda=\lambda^{*}_{1SE MAE}',...
        'Distribution of the estimated coefficients at \lambda=\lambda^{*}_{RMSE}',...
        'Distribution of the estimated coefficients at \lambda=\lambda^{*}_{1SE RMSE}'};
    for i = 1:length(crits)
        DSTEM_HDGM_CV_PenLik_OptRes_Plot(...
            CV_perf_metrics,beta_true,fda_setup,input_fda,...
            {'Intercept','X1','X2','X3'},...
            {'E'},...
            [],[],[],...
            crits{i},box_tit{i},...
            sup_text_tit,...
            1,save_plot_path,save_plot_name);
    end
    close all;
end