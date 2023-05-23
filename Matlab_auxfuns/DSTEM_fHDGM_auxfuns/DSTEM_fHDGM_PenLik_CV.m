%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%% DSTEM_fHDGM_penlik_KfoldCV_par %%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%%% Computes the CV (using a CVpart object) for a functional HDGM with penalization
%%% of the regression coefficients

function [CV_perf_metrics] = DSTEM_fHDGM_PenLik_CV(cv_part,DSTEM_model,...
    InputData_fHDGM,lambda_seq,beta_dest_flag,Fisher_opts,ActiveSets,ASi_idx)

if 0
    DSTEM_model = obj_stem_model;
    cv_part = cv_part;
    ActiveSets = [];
    ASi_idx = [];
    InputData_fHDGM = input_data;
    Fisher_opts = opts;
    lambda_seq = lambda_seq;
end

if isempty(ActiveSets) == 1
    Variable = DSTEM_model.stem_data.X_beta_name{1,1}';
    ActiveSets = array2table(Variable);
    ActiveSets.Response_idx = repmat(1,size(ActiveSets,1),1);
    ActiveSets.ActiveSets1 = repmat(1,size(ActiveSets,1),1);
    ASi_idx = 1;
end

ASi_cols = find(contains(ActiveSets.Properties.VariableNames,'ActiveSets'));
nameAS = ActiveSets.Properties.VariableNames(contains(ActiveSets.Properties.VariableNames,'ActiveSets'));
nvarY = 1;
namevarY = {'Y'};
if contains(cv_part.cv_type,'LeaveOneStatOut')
    varY_pos = horzcat(cv_part.idx_st_poll{:});
else
    varY_pos = ones(size(cv_part.ARPA_stats_reg,1),nvarY);
end
MSE = cell(nvarY,cv_part.K);
RMSE = cell(nvarY,cv_part.K);
MAE = cell(nvarY,cv_part.K);
R2 = cell(nvarY,cv_part.K);
beta_LASSO_LQA = cell(1,cv_part.K);
MSE_full = cell(nvarY,1);
RMSE_full = cell(nvarY,1);
MAE_full = cell(nvarY,1);
R2_full = cell(nvarY,1);
beta_LASSO_LQA_full = cell(1,1);
fold_failed = zeros(1,cv_part.K);

parfor run = 1:(cv_part.K + 1)
    if run ~= cv_part.K + 1
        sprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',...
            'Fold %d of %d',...
            ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%'],run,cv_part.K)
    else
        sprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',...
            'Full sample re-estimation ',...
            ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%'])
    end

    %%%%% Extracting CV objects
    Y_obs = struct();
    Y_obs = cv_part.Y;
    if run ~= cv_part.K + 1
        test_coords = cv_part.test_coords{1,run};
    end
    
    %%%%% STEM settings
    % Extracting input data from the simulation setting
    input_data_temp = InputData_fHDGM;
    obj_stem_fda = stem_fda(input_data_temp.stem_fda);
    % Changing the fda settings to the new estimation setting
    input_data_temp.stem_fda = obj_stem_fda;
    % CVpartition data
    if run ~= cv_part.K + 1
        input_data_temp.data_table.Y = cv_part.Y_fun_cv{1,run}{1,1};
    else
        input_data_temp.data_table.Y = cv_part.Y_fun{1,1};
    end
    % Creation of stem_data object
    obj_temp_stem_data = stem_data(input_data_temp);
    % Creation of stem_par object
    obj_temp_stem_par = DSTEM_model.stem_par;
    % Creation of stem_model object
    obj_temp_stem_model = stem_model(obj_temp_stem_data, obj_temp_stem_par);
    
    %%%%% Data transform
    if DSTEM_model.stem_data.log_transformed == 1
        obj_temp_stem_model.stem_data.log_transform;
    end
    if DSTEM_model.stem_data.standardized == 1
        obj_temp_stem_model.stem_data.standardize;
    end
    
    %%%%% Model estimation
    %%% Starting values
    obj_temp_stem_par.beta = DSTEM_model.stem_EM_result.stem_par.beta;
    obj_temp_stem_par.theta_z = DSTEM_model.stem_EM_result.stem_par.theta_z;
    obj_temp_stem_par.v_z = DSTEM_model.stem_EM_result.stem_par.v_z;
    obj_temp_stem_par.sigma_eta = DSTEM_model.stem_EM_result.stem_par.sigma_eta;
    obj_temp_stem_par.G = DSTEM_model.stem_EM_result.stem_par.G;
    obj_temp_stem_par.sigma_eps = DSTEM_model.stem_EM_result.stem_par.sigma_eps;
    obj_temp_stem_model.set_initial_values(obj_temp_stem_par);
    %%% Parameters estimation
    obj_temp_stem_EM_options = stem_EM_options();
    obj_temp_stem_EM_options.exit_tol_loglike = DSTEM_model.stem_EM_result.exit_tol_loglike;
    obj_temp_stem_EM_options.exit_tol_par = DSTEM_model.stem_EM_result.exit_tol_par;
    obj_temp_stem_EM_options.max_iterations = DSTEM_model.stem_EM_result.iterations;
    %%% CV for loop and performances storage
    y_oos = struct();
    y_oos.MSE = cell(nvarY,1);
    y_oos.RMSE = cell(nvarY,1);
    y_oos.MAE = cell(nvarY,1);
    y_oos.R2 = cell(nvarY,1);
    for p = find(varY_pos(run,:))
        y_oos.MSE{1,p} = zeros(length(lambda_seq),1);
        y_oos.RMSE{1,1} = zeros(length(lambda_seq),1);
        y_oos.MAE{1,1} = zeros(length(lambda_seq),1);
        y_oos.R2{1,1} = zeros(length(lambda_seq),1);
    end
    
    try
        obj_temp_stem_model.EM_estimate(obj_temp_stem_EM_options);
        obj_temp_stem_model.set_logL;
        %%% Var-cov matrix of estimated pars for each k-th fold
        obj_temp_stem_model.set_varcov;
        if beta_dest_flag == 1
            [MLE] = DSTEM_extraction(DSTEM_model,[],[]);
            [MLE_temp] = DSTEM_extraction(obj_temp_stem_model,[],[]);
            MLE_logL = obj_temp_stem_model.stem_EM_result.logL;
            %%% Penalized ML estimates (Adaptive LASSO)
            beta_LQA = zeros(length(obj_temp_stem_model.stem_EM_result.stem_par.beta),length(lambda_seq));
            % 			for lam = 1:length(lambda_seq)
            for lam = [1,3,10,20,51]
                sprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',...
                    'Fold %d : Lambda %d of %d',...
                    ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%'],run,lam,length(lambda_seq))
                obj_temp_stem_model.stem_EM_result.stem_par.beta = MLE_temp.Reg_pars.Beta_destd_tab.Y.Coef;
                obj_temp_stem_model.stem_EM_result.stem_par.varcov(1:length(obj_temp_stem_par.beta),1:length(obj_temp_stem_par.beta)) = table2array(MLE_temp.Reg_pars.VarCov_reg_destd.Y);
                %                 obj_temp_stem_model.stem_EM_result.logL = MLE_logL;
                [beta_adaptive] = DSTEM_Fisher_PenLik(obj_temp_stem_model,lambda_seq(lam),Fisher_opts);
                beta_LQA(:,lam) = beta_adaptive;
                [repelem(lam,28)' repelem(lambda_seq(lam),28)' beta_adaptive obj_temp_stem_model.stem_par.beta obj_temp_stem_model.stem_EM_result.stem_par.beta]
            end
        else
            %%% Penalized ML estimates (Adaptive LASSO)
            beta_LQA = zeros(length(DSTEM_model.stem_EM_result.stem_par.beta),length(lambda_seq));
            for lam = 1:length(lambda_seq)
                sprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',...
                    'Fold %d : Lambda %d of %d',...
                    ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%'],run,lam,length(lambda_seq))
                obj_temp_stem_model.stem_par.beta = DSTEM_model.stem_EM_result.stem_par.beta;
                [beta_adaptive] = DSTEM_Fisher_PenLik(obj_temp_stem_model,lambda_seq(lam),Fisher_opts);
                beta_LQA(:,lam) = beta_adaptive;
            end
        end
        %%% Iterate for each lambda in lambda_seq
        if run ~= cv_part.K + 1
            for lam = 1:length(lambda_seq)
                %%% Calculating fitted values yhat(s,t)_lambda,k|betahat_lambda,k
                [output_y] = DSTEM_fitted_given_beta(...
                    obj_temp_stem_model,obj_temp_stem_EM_options,beta_LQA(:,lam));
                yhat_mat_temp = output_y.y_hat_back_mat;
                %%% Store CV performance indices for each run for a given lambda
                for p = find(varY_pos(run,:))
                    for i = 1:size(cv_part.test_coords{1,run}{1,p},1)
                        y_oos.(namevarY{p})(i,1) = cv_part.test_coords{1,run}{1,p}(i,1);
                        y_oos.(namevarY{p})(i,2) = cv_part.test_coords{1,run}{1,p}(i,2);
                        y_oos.(namevarY{p})(i,3) = Y_obs{1,p}(cv_part.test_coords{1,run}{1,p}(i,1),...
                            cv_part.test_coords{1,run}{1,p}(i,2));
                        y_oos.(namevarY{p})(i,4) = yhat_mat_temp(cv_part.test_coords{1,run}{1,p}(i,1),...
                            cv_part.test_coords{1,run}{1,p}(i,2));
                    end
                    y_oos.(namevarY{p}) = array2table(y_oos.(namevarY{p}));
                    y_oos.(namevarY{p}).Properties.VariableNames = {'row','col','y_obs','y_hat'};
                    y_oos.MSE{p,1}(lam) = nanmean((y_oos.(namevarY{p}).y_hat - y_oos.(namevarY{p}).y_obs).^2);
                    y_oos.RMSE{p,1}(lam) = sqrt(y_oos.MSE{p,1}(lam));
                    y_oos.MAE{p,1}(lam) = nanmean(abs(y_oos.(namevarY{p}).y_hat - y_oos.(namevarY{p}).y_obs));
                    y_oos.R2{p,1}(lam) = 1 - y_oos.MSE{p,1}(lam) ./ nanvar(y_oos.(namevarY{p}).y_obs);
                    y_oos.(namevarY{p}) = [];
                end
            end
        else
            for lam = 1:length(lambda_seq)
                %%% Calculating fitted values yhat(s,t)_lambda,k|betahat_lambda,k
                [output_y] = DSTEM_fitted_given_beta(...
                    obj_temp_stem_model,obj_temp_stem_EM_options,beta_LQA(:,lam));
                yhat_mat_temp = output_y.y_hat_back_mat;
                for p = find(varY_pos(run,:))
                    y_oos.MSE{p,1}(lam) = nanmean((Y_obs{1,p} - yhat_mat_temp).^2,'all');
                    y_oos.RMSE{p,1}(lam) = sqrt(y_oos.MSE{p,1}(lam));
                    y_oos.MAE{p,1}(lam) = nanmean(abs(Y_obs{1,p} - yhat_mat_temp),'all');
                    y_oos.R2{p,1}(lam) = 1 - y_oos.MSE{p,1}(lam) ./ nanvar(Y_obs{1,p}(:));
                end
            end
        end
        sprintf('Fold %d complete', run)
        AIC(run) = obj_temp_stem_model.stem_EM_result.AIC;
        LogL(run) = obj_temp_stem_model.stem_EM_result.logL;
        MSE(:,run) = y_oos.MSE;
        RMSE(:,run) = y_oos.RMSE;
        MAE(:,run) = y_oos.MAE;
        R2(:,run) = y_oos.R2;
        beta_LASSO_LQA(:,run) = {beta_LQA};
    catch
        sprintf('Fold %d failed', run)
        fold_failed(run) = 1;
        AIC(run) = NaN;
        LogL(run) = NaN;
        MSE(:,run) = y_oos.MSE;
        RMSE(:,run) = y_oos.RMSE;
        MAE(:,run) = y_oos.MAE;
        R2(:,run) = y_oos.R2;
        continue
    end
end

%%% Store performance metrics
% General information
CV_perf_metrics.Lambda_seq = lambda_seq;
CV_perf_metrics.Fold_failed = fold_failed;
CV_perf_metrics.NameVarY = namevarY;
CV_perf_metrics.NameAS = nameAS;
% K-fold CV metrics
CV_perf_metrics.AIC = AIC(1:cv_part.K);
CV_perf_metrics.LogL = LogL(1:cv_part.K);
CV_perf_metrics.MSE = MSE(:,1:cv_part.K);
CV_perf_metrics.RMSE = RMSE(:,1:cv_part.K);
CV_perf_metrics.MAE = MAE(:,1:cv_part.K);
CV_perf_metrics.R2 = R2(:,1:cv_part.K);
CV_perf_metrics.beta_LASSO_LQA = beta_LASSO_LQA(:,1:cv_part.K);
% Full sample metrics
CV_perf_metrics.Full_sample.AIC = AIC(cv_part.K + 1);
CV_perf_metrics.Full_sample.LogL = LogL(cv_part.K + 1);
CV_perf_metrics.Full_sample.MSE = MSE(:,cv_part.K + 1);
CV_perf_metrics.Full_sample.RMSE = RMSE(:,cv_part.K + 1);
CV_perf_metrics.Full_sample.R2 = R2(:,cv_part.K + 1);
CV_perf_metrics.Full_sample.beta_LASSO_LQA = beta_LASSO_LQA(:,cv_part.K + 1);

return;
end