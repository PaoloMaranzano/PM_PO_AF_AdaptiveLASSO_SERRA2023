%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%% DSTEM_fHDGM_KfoldCV_par %%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%%% Computes the random K-fold CV for a generic functional HDGM

function [CV_perf_metrics] = DSTEM_fHDGM_CV(cv_part,DSTEM_model,...
    InputData_fHDGM,ActiveSets,ASi_idx)

if 0
    DSTEM_model = stem_model_LASSO;
    cv_part = cv_part;
    ActiveSets = [];
    InputData_fHDGM = input_data;
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
MSE_full = cell(nvarY,1);
RMSE_full = cell(nvarY,1);
MAE_full = cell(nvarY,1);
R2_full = cell(nvarY,1);
fold_failed = zeros(1,cv_part.K);

parfor run = 1:(cv_part.K + 1)
    sprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',...
    'Fold %d of %d (except the full sample fold)',...
    ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%'],run,cv_part.K)
    
    %%%%% Extracting CV objects
    Y_obs = struct();
    Y_obs = cv_part.Y;
    % [~,~,Y_obs] = DSTEM_extract_Y_to_funY(DSTEM_model,DSTEM_obj_sim_str.data_long.gId,0);
    test_coords = cv_part.test_coords{1,run};
    
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
        input_data_temp.data_table.Y = cv_part.Y{1,1};
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
    obj_temp_stem_par.beta = DSTEM_model.stem_EM_result.stem_par.beta(find(ActiveSets{:,ASi_cols(ASi_idx)}));
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
    try
        obj_temp_stem_model.EM_estimate(obj_temp_stem_EM_options);
        obj_temp_stem_model.set_logL;
        [~,~,yhat_mat_temp] = DSTEM_extract_Y_to_funY(obj_temp_stem_model,input_data_temp.data_long.gId,1);
        % [~,~,y_mat_temp] = DSTEM_extract_Y_to_funY(obj_temp_stem_model,DSTEM_obj_sim_str.data_long.gId,0);
        %%%%% Store CV performance indices for each run for a given Active Set
        if run ~= cv_part.K + 1
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
                y_oos.MSE{p,1} = nanmean((y_oos.(namevarY{p}).y_hat - y_oos.(namevarY{p}).y_obs).^2);
                y_oos.RMSE{p,1} = sqrt(y_oos.MSE{p,1});
                y_oos.MAE{p,1} = nanmean(abs(y_oos.(namevarY{p}).y_hat - y_oos.(namevarY{p}).y_obs));
                y_oos.R2{p,1} = 1 - y_oos.MSE{p,1} ./ nanvar(y_oos.(namevarY{p}).y_obs);
            end
        else
            y_oos.MSE{p,1}(lam) = nanmean((Y_obs{1,p} - yhat_mat_temp).^2);
            y_oos.RMSE{p,1}(lam) = sqrt(y_oos.MSE{p,1}(lam));
            y_oos.MAE{p,1}(lam) = nanmean(abs(Y_obs{1,p} - yhat_mat_temp));
            y_oos.R2{p,1}(lam) = 1 - y_oos.MSE{p,1}(lam) ./ nanvar(y_oos.(namevarY{p}).y_obs);
        end
        sprintf('Fold %d complete', run)
        if run ~= cv_part.K + 1
            AIC(run) = obj_temp_stem_model.stem_EM_result.AIC;
            LogL(run) = obj_temp_stem_model.stem_EM_result.logL;
            MSE(:,run) = y_oos.MSE;
            RMSE(:,run) = y_oos.RMSE;
            MAE(:,run) = y_oos.MAE;
            R2(:,run) = y_oos.R2;
        else
            AIC_full = obj_temp_stem_model.stem_EM_result.AIC;
            LogL_full = obj_temp_stem_model.stem_EM_result.logL;
            MSE_full = y_oos.MSE;
            RMSE_full = y_oos.RMSE;
            MAE_full = y_oos.MAE;
            R2_full = y_oos.R2;
        end
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
CV_perf_metrics.NameVarY = namevarY;
CV_perf_metrics.NameAS = nameAS;
% K-fold CV metrics
CV_perf_metrics.AIC = AIC;
CV_perf_metrics.LogL = LogL;
CV_perf_metrics.MSE = MSE;
CV_perf_metrics.RMSE = RMSE;
CV_perf_metrics.MAE = MAE;
CV_perf_metrics.R2 = R2;
% Full sample metrics
CV_perf_metrics.Full_sample.AIC = AIC_full;
CV_perf_metrics.Full_sample.LogL = LogL_full;
CV_perf_metrics.Full_sample.MSE = MSE_full;
CV_perf_metrics.Full_sample.RMSE = RMSE_full;
CV_perf_metrics.Full_sample.R2 = R2_full;

return;
end