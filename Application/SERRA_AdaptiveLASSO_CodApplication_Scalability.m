%% Paper: Adaptive LASSO estimation for functional hidden dynamic geostatistical models
%% Journal: Stochastic Environmental Research and Risk Assessment
%% Date: July 2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Application of PenLik to functional HDGM model: Scalability analysis %% %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% Change folder
cd 'H:/.shortcut-targets-by-id/1hRMzSqAcE5AcmrsOu5k3nw1rmUY1NHyO/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application'
% cd('C:/Users/paulm/Google Drive/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application')

% Auxiliary functions and data
addpath(genpath('../../../../VisitingLUH2021'));
addpath(genpath('../../../../SPASTA2021/Code/DSTEM_software'));
addpath(genpath('../../../../SPASTA2021/Code/Matlab_auxfuns'));
addpath(genpath('../../../../SPASTA2021/Data'));

% Load original data
load('Application_data_ENV.mat')

% Export data in long format for functional boxplot
% Long_exp = Ground.data_long(:,{'Date_day','y','m','d','h','IDStat','NO2',...
%     'Pressure','RelHumid','Temperature','Rainfall','WindU','WindV'});
% writetable(Long_exp,'LongData.csv')

% Models combinations
comb = readtable("Models_combinations.xlsx");

comb = comb(1,:);


for mod = 1:size(comb,1)

    % Optimal mod = 1

    % Load application results
    SpPart_k = comb(mod,:).Spat_part;       % 1, 2, 3, 4 or 5
    nbasis = comb(mod,:).Basis_nbr;         % 5 or 7 or 9
    VarCov_type = comb(mod,:).VarCov{:};    % Exact or Approx
    load(['Application_ENV_SpPart' num2str(SpPart_k) '_VarCov' VarCov_type '_nbasis' num2str(nbasis) '.mat'])

    % Destination folder
    save_plot_path = ['VarCov' VarCov_type '_nbasis' num2str(nbasis) '_SpPart' num2str(SpPart_k) '/'];
    if ~exist(save_plot_path)
        mkdir(save_plot_path)
    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %% %%%%% Functional box-plot of NO2 %% %%
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    X_name = 'Hour';
    Y_name = 'NO_2';
    sup_tit_text = 'Functional box-plot of NO_2 concentrations';
    sub_tit_text = 'Hourly data from 1st March 2020 to 31th May 2020';
    X_lab = 'Hour of the day';
    Y_lab = 'NO_2 [\mu g/m^3 ]';
    save_plot_name = 'Application';
    data_long_boxplot = Ground.data_long(:,{'Date_day','y','m','d','h','IDStat','NO2'});
    data_long_boxplot.Properties.VariableNames{'NO2'} = 'NO_2';
    data_long_boxplot.Properties.VariableNames{'h'} = 'Hour';
    [stat_grp] = Functional_BoxPlot(data_long_boxplot,X_name,Y_name,...
        sup_tit_text,sub_tit_text,X_lab,Y_lab);
    saveas(gcf,[save_plot_path save_plot_name '_Y.png'])
    close all




    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %% %%%%% Functional box-plots of covariates %% %%
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    MC_CV_perf_met{1,1} = CV_perf_metrics;
    beta_true = [];
    fda_setup_true = [];
    vars_names = obj_stem_model.stem_data.stem_varset_p.X_beta_name{1,1}';
    sup_text_tit = 'Application - AQ in Lombardy during COVID-19 pandemic';
    input_fda.n_basis = n_basis;

    [optimal_vals] = DSTEM_HDGM_CV_PenLik_SimsOptRes(MC_CV_perf_met,...
        beta_true,vars_names,1,1,save_plot_path);
    %  {'A','B','C','D','F'},...
    DSTEM_HDGM_CV_PenLik_OptRes_Plot(...
        optimal_vals,[],[],input_fda,...
        vars_names,...
        {'A','B','C','D','F'},...
        1,0,0,...
        [],[],...
        sup_text_tit,...
        0,save_plot_path,save_plot_name);
    close all




    %% %%%%% Spatio-temporal parameters and variance plots
    %%% Computing basis splines
    tobs = input_fda.spline_range(1):0.01:input_fda.spline_range(2);

    sigma_eps = obj_stem_model.stem_par.sigma_eps;
    theta = obj_stem_model.stem_par.theta_z';
    G = diag(obj_stem_model.stem_par.G);
    vz = diag(obj_stem_model.stem_par.v_z);

    if input_fda.spline_type == "Fourier"
        basis_obj = create_fourier_basis(input_fda.spline_range, input_fda.n_basis);
        basis_obj_mat = eval_basis(tobs, basis_obj);
        basis_mat = full(basis_obj_mat);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_sigma_eps] = simulate_fourier_basis(...
            input_fda,sigma_eps,[],[]);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_theta] = simulate_fourier_basis(...
            input_fda,theta,[],[]);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_G] = simulate_fourier_basis(...
            input_fda,G,[],[]);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_vz] = simulate_fourier_basis(...
            input_fda,vz,[],[]);
    else
        basis_obj = create_bspline_basis(input_fda.spline_range, input_fda.n_basis,...
            input_fda.spline_order, input_fda.spline_knots);
        basis_obj_mat = eval_basis(tobs, basis_obj);
        basis_mat = full(basis_obj_mat);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_sigma_eps] = simulate_bspline_basis(...
            input_fda,sigma_eps,[],[]);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_theta] = simulate_bspline_basis(...
            input_fda,theta,[],[]);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_G] = simulate_bspline_basis(...
            input_fda,G,[],[]);
        [basis_obj,basis_obj_mat,basis_mat,yhat_MLE_vz] = simulate_bspline_basis(...
            input_fda,vz,[],[]);
    end

    if 1
        % Significance level
        alpha = 0.05;
        %%% Computing functional coefficients
        sigma_eps_h = exp(yhat_MLE_sigma_eps);
        eps_idx = (nbasis*length(vars_names) + 0*nbasis + 1) : (nbasis*length(vars_names) + 0*nbasis + nbasis);
        Varcov_eps = obj_stem_model.stem_par.varcov(eps_idx,eps_idx);
        Varcov = diag(basis_mat*Varcov_eps*basis_mat');
        Varcov1 = (exp(Varcov)-ones(length(Varcov),1)).*((sigma_eps_h).^2).*(exp(Varcov));
        sigma_eps_se = sqrt(Varcov1);
        sigma_eps_h_u = sigma_eps_h + norminv(1 - alpha/2)*sigma_eps_se;
        sigma_eps_h_l = sigma_eps_h - norminv(1 - alpha/2)*sigma_eps_se;
        % Theta
        theta_idx = (nbasis*length(vars_names) + 1*nbasis + 1) : (nbasis*length(vars_names) + 1*nbasis + nbasis);
        Varcov_theta = obj_stem_model.stem_par.varcov(theta_idx,theta_idx);
        Varcov = diag(basis_mat*Varcov_theta*basis_mat');
        Varcov1 = Varcov;
        theta_se = sqrt(Varcov1);
        theta_h_u = yhat_MLE_theta + norminv(1 - alpha/2)*theta_se;
        theta_h_l = yhat_MLE_theta - norminv(1 - alpha/2)*theta_se;
        % G
        G_idx = (nbasis*length(vars_names) + 2*nbasis + 1) : (nbasis*length(vars_names) + 2*nbasis + nbasis);
        Varcov_G = obj_stem_model.stem_par.varcov(G_idx,G_idx);
        Varcov = diag(basis_mat*Varcov_G*basis_mat');
        Varcov1 = Varcov;
        G_se = sqrt(Varcov1);
        G_h_u = yhat_MLE_G + norminv(1 - alpha/2)*G_se;
        G_h_l = yhat_MLE_G - norminv(1 - alpha/2)*G_se;
        % V_z
        vz_idx = (nbasis*length(vars_names) + 3*nbasis + 1) : (nbasis*length(vars_names) + 3*nbasis + nbasis);
        Varcov_vz = obj_stem_model.stem_par.varcov(vz_idx,vz_idx);
        Varcov = diag(basis_mat*Varcov_vz*basis_mat');
        Varcov1 = Varcov;
        vz_se = sqrt(Varcov1);
        vz_h_u = yhat_MLE_vz + norminv(1 - alpha/2)*vz_se;
        vz_h_l = yhat_MLE_vz - norminv(1 - alpha/2)*vz_se;
        %%% Plotting functional coefficients
        f = figure('Name','plot G');
        f.WindowState = 'maximized';
        sup_t1 = sgtitle([sup_text_tit],...
            'Color','black','fontweight','bold','fontsize',20);
        subplot(1,1,1)
        plot(tobs,sigma_eps_h,'Color','black','linewidth',2)
        patch([fliplr(tobs) tobs],[fliplr(sigma_eps_h_u') sigma_eps_h_l'],'r','FaceAlpha',0.35,'EdgeAlpha',0)
        title(['\sigma^2_\epsilon(h)'])
        % subplot(2,2,2)
        % plot(tobs,deg2km(yhat_MLE_theta),'Color','black','linewidth',2)
        % patch([fliplr(tobs) tobs],[fliplr(deg2km(theta_h_u')) deg2km(theta_h_l')],'r','FaceAlpha',0.35,'EdgeAlpha',0)
        % title(['\theta_Z(h) in km'])
        % subplot(2,2,3)
        % plot(tobs,yhat_MLE_G,'Color','black','linewidth',2)
        % patch([fliplr(tobs) tobs],[fliplr(G_h_u') G_h_l'],'r','FaceAlpha',0.35,'EdgeAlpha',0)
        % title(['G(h)'])
        % subplot(2,2,4)
        % plot(tobs,yhat_MLE_vz,'Color','black','linewidth',2)
        % patch([fliplr(tobs) tobs],[fliplr(vz_h_u') vz_h_l'],'r','FaceAlpha',0.35,'EdgeAlpha',0)
        % title(['V_Z(h) in km^2'])
        [ax1,h1]=suplabel('Time (hours)','x');
        [ax2,h2]=suplabel('Functional coefficients over the 24 hours','y');
        saveas(gcf,[save_plot_path save_plot_name '_G.png'])
        %%% Table
        [ G , sqrt(diag(Varcov_G)) , theta , sqrt(diag(Varcov_theta)) , vz , sqrt(diag(Varcov_vz)) ]
        close all;
    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %% %%%%% Functional coefficients at optimal values %% %%
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    tab_coefs = [optimal_vals.Avg_after_opt.Full_sample.beta_LASSO_LQA.beta_min_MAE , ...
        optimal_vals.Avg_after_opt.Full_sample.beta_LASSO_LQA.beta_1se_MAE, ...
        optimal_vals.Avg_after_opt.Full_sample.beta_LASSO_LQA.beta_min_RMSE, ...
        optimal_vals.Avg_after_opt.Full_sample.beta_LASSO_LQA.beta_1se_RMSE];
    tab_coefs = array2table(tab_coefs);
    betas = optimal_vals.Coef_names';
    betas = array2table(betas);
    spline_nbasis_fourier = size(betas,1) / length(vars_names);
    var = repelem(vars_names,spline_nbasis_fourier)';
    var = array2table(var);
    beta_MLE = array2table(obj_stem_model.stem_EM_result.stem_par.beta);
    tab_coefs = [var,betas,beta_MLE,tab_coefs];
    tab_coefs.Properties.VariableNames = {'Variable','Beta',...
        'beta_MLE',...
        'beta_min_MAE','beta_1se_min_MAE',...
        'beta_min_RMSE','beta_1se_min_RMSE'};
    writetable(tab_coefs,[save_plot_path 'Coefficients.csv']);



    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    %% %%%%% Functional box-plots of the covariates %% %%
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    Xregs = {'Pressure','RelHumid','Temperature','Rainfall','WindU','WindV'};
    units = {'hPa','%','°','mm','m/s','m/s'};
    f = figure('Name','plot H');
    f.WindowState = 'maximized';
    sup_t1 = sgtitle('Functional box-plots',...
        'Color','black','fontweight','bold','fontsize',20);
    for i = 1:length(Xregs)
        subplot(3,2,i)
        X_name = 'Hour';
        Y_name = Xregs{i};
        sup_tit_text = Xregs{i};
        sub_tit_text = 'Hourly data from 1st March 2020 to 31th May 2020';
        X_lab = 'Hour of the day';
        Y_lab = units{i};
        data_long_boxplot = Ground.data_long(:,{'Date_day','y','m','d','h','IDStat',Xregs{i}});
        data_long_boxplot.Properties.VariableNames{Xregs{i}} = Xregs{i};
        data_long_boxplot.Properties.VariableNames{'h'} = 'Hour';
        [stat_grp] = Functional_BoxPlot(data_long_boxplot,X_name,Y_name,...
            sup_tit_text,sub_tit_text,X_lab,Y_lab);
    end
    saveas(gcf,[save_plot_path save_plot_name '_H.png'])
    close all;

    comb_out = comb(mod,:);

    %%%%%%%%%% Memory
    det = dir(['Application_ENV_SpPart' num2str(SpPart_k) '_VarCov' VarCov_type '_nbasis' num2str(nbasis) '.mat']);

    %%%%%%%%%% Timing
    if ~exist('End_time','var')
        End_time = datetime(det.date,"InputFormat","dd-MMM-uuuu HH:mm:ss","Locale","it_IT");
    end

    out_vals = [minutes(VarCov_time_begin - Begin_time), ...         %%% Time DSTEM
        minutes(PenLik_time_begin - VarCov_time_begin), ...     %%% Time VarCov
        minutes(End_time - PenLik_time_begin), ...              %%% Time PenLik
        minutes(End_time - Begin_time),...                      %%% Total time
        det.bytes/(1024^2), ...
        optimal_vals.Avg_before_opt.Optimal_vals.min_RMSE , ...
        optimal_vals.Avg_before_opt.Optimal_vals.lambda_min_RMSE, ...
        optimal_vals.Avg_before_opt.Optimal_vals.oneSE_min_RMSE , ...
        optimal_vals.Avg_before_opt.Optimal_vals.lambda_oneSE_min_RMSE, ...
        optimal_vals.Avg_before_opt.Optimal_vals.min_MAE , ...
        optimal_vals.Avg_before_opt.Optimal_vals.lambda_min_MAE, ...
        optimal_vals.Avg_before_opt.Optimal_vals.oneSE_min_MAE,...
        optimal_vals.Avg_before_opt.Optimal_vals.lambda_oneSE_min_MAE];
    comb_out = [comb_out , array2table(out_vals)];

    output_tab(mod,:) = comb_out;
end

output_tab.Properties.VariableNames = {'VarCov','Spatial partitioning', 'Basis number (p=x10+4)', ...
    'Computation time - DSTEM (Minutes)','Computation time - VarCov (Minutes)', ...
    'Computation time - PenLik (Minutes)', 'Computation time - Total (Minutes)', ...
    'Storage dimension (MB)', ...
    'min RMSE', 'Lambda min RMSE', '1-SE min RMSE', 'Lambda 1-SE min RMSE',...
    'min MAE', 'Lambda min MAE', '1-SE min MAE', 'Lambda 1-SE min MAE'};

writetable(output_tab,'Models_summary_auto.xlsx');


