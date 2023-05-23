%% Paper: Adaptive LASSO estimation for functional hidden dynamic geostatistical models
%% Journal: Stochastic Environmental Research and Risk Assessment
%% Date: July 2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Application of PenLik to functional HDGM model: Estimation %% %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
Begin_time = datetime(datestr(now))

%% %%%%% Set working directory
% cd('C:/Users/paulm/Google Drive/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application')

%% %%%%% Auxiliary functions
addpath(genpath('../../../../VisitingLUH2021'));
% addpath(genpath('../../../SPASTA2021/Code/DSTEM_software'));
% addpath(genpath('../../../SPASTA2021/Code/Matlab_auxfuns'));
addpath(genpath('../../../DSTEM_Code_General/DSTEM_software'));
addpath(genpath('../../../DSTEM_Code_General/Matlab_auxfuns'));

%% %%%%% Generation / Load data
DataBulding = 0;
if DataBulding == 1
    Data_generation_time_begin = datetime(datestr(now));
    run('Application_AQ_COVID_DataManagement.m');
else
    Data_loading_time_begin = datetime(datestr(now));
    load('Application_data_ENV.mat');
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%      Settings      %%
%%%%%%%%%%%%%%%%%%%%%%%%
Setting_fHDGM_time_begin = datetime(datestr(now))
% Spline type
spline_type = 'Fourier';
spline_nbasis_fourier = 5;
% Estimate varcov of the full model
varcov_estim_flag = 1;
% Estimate the log-likelihood of the full model
logL_estim_flag = 1;
% Standardize the data before estimation
standardization_flag = 1;
% Use de-standardized coefficients and varcov for PMLE
beta_dest_flag = 0;
% Penalize the intercept or not
beta0_pen_flag = 0;
% Exact or approximated VarCov computation
ExactVarcov = 0;
delta_varcov = 0.01;
% Spatial partitioning
SpatPart = 1;
k = 4;
trials = 100;
lambda = 5000;
% Save the final output
save_output_flag = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Model estimation      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Setup inputs
% FDA inputs
n_covs = length(Ground.vars_names);
nvars = n_covs + 1;
input_fda.spline_type = spline_type;        % Fourier (periodic) or B-splines basis
input_fda.spline_range = [0 24];            % Spline basis domain
if contains(spline_type,'Bspline')
    input_fda.spline_order = 3;
    input_fda.knots_number = 5;
    input_fda.spline_knots = linspace(input_fda.spline_range(1),input_fda.spline_range(2),input_fda.knots_number);
    n_basis = (input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1))*nvars;
elseif contains(spline_type,'Fourier')
    input_fda.spline_nbasis_z = spline_nbasis_fourier;
    input_fda.spline_nbasis_beta = spline_nbasis_fourier;
    input_fda.spline_nbasis_sigma = spline_nbasis_fourier;
    n_basis = input_fda.spline_nbasis_beta*nvars;
end
obj_stem_fda = stem_fda(input_fda);
input_fda.n_basis = n_basis;
input_fda.n_vars = nvars;
% stem_object creation
%%% Model type
obj_stem_modeltype = stem_modeltype('f-HDGM');
input_data.stem_modeltype = obj_stem_modeltype;
input_data.data_table = Ground.data_fhdgm;
input_data.data_long = Ground.data_long;
% Changing the fda settings to the new estimation setting
input_data.stem_fda = obj_stem_fda;
% Creation of stem_data object
obj_stem_data = stem_data(input_data);
% Creation of stem_par object
obj_stem_par_constraints = stem_par_constraints();
obj_stem_par_constraints.time_diagonal = 0;
obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
% Creation of stem_model object
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);

%%% Data transform
if standardization_flag == 1
    obj_stem_model.stem_data.standardize;
end

%%% Starting values of EM algorithm
obj_stem_par.beta = obj_stem_model.get_beta0();
obj_stem_par.theta_z = repelem(0.5,n_basis/nvars);
obj_stem_par.v_z = diag(repelem(1,n_basis/nvars));
obj_stem_par.sigma_eta = repelem(1,n_basis/nvars)';
obj_stem_par.G = diag(repelem(0.85,n_basis/nvars));
obj_stem_par.sigma_eps = repelem(1,n_basis/nvars)';
obj_stem_model.set_initial_values(obj_stem_par);

% Parameters estimation
Estimation_time_begin = datetime(datestr(now))
obj_stem_EM_options = stem_EM_options();
obj_stem_EM_options.exit_tol_loglike = 0.0001;
obj_stem_EM_options.exit_tol_par = 0.0001;
obj_stem_EM_options.max_iterations = 300;

% Spatial partitioning options
if SpatPart == 1
partitions = obj_stem_data.kmeans_partitioning(k, trials, lambda);
obj_stem_EM_options.partitions = partitions;
obj_stem_EM_options.workers = 2;
end

% Estimation
obj_stem_model.EM_estimate(obj_stem_EM_options);

% Var-cov matrix of estimated pars
if varcov_estim_flag == 1
    VarCov_time_begin = datetime('now')
	if ExactVarcov == 1
	    % Exact
		obj_stem_model.set_varcov;
	else
		% Approximated
		delta_varcov = delta_varcov;
		obj_stem_model.set_varcov(delta_varcov);
	end
end

% Log-likelihood of the data
if logL_estim_flag == 1
    LogLik_time_begin = datetime(datestr(now))
    obj_stem_model.set_logL;
    logL_estim = obj_stem_model.stem_EM_result.logL
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Model selection using penalized likelihood      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fisher algorithm settings
beta_MLE = obj_stem_model.stem_EM_result.stem_par.beta;
if beta0_pen_flag == 1
    opts.penaltywt = ones(length(beta_MLE),1);
else
    opts.penaltywt = [zeros(n_basis/nvars,1) ; ones(n_basis/nvars*n_covs,1)];
end
penalty_fun = "adaptive";
if ismember(penalty_fun,'adaptive')
    % For adaptive LASSO
    opts.penalty = @p_adaptive;
    opts.gamma = 1;
    opts.adaptivewt = beta_MLE;
end
if ismember(penalty_fun,'lasso')
    % For LASSO
    opts.penalty = @p_lasso;
end
if ismember(penalty_fun,'ridge')
    % For ridge
    opts.penalty = @p_ridge;
end
% Local Quadratica Approximation for the log-likelihood
opts.LQA_flag = 1;
% Re-estimation of VarCov for each CV K-fold
opts.CV_VarCov_reestim_flag = 0;
% Optimizer options
opts.opts_optim = optimoptions('fminunc',...
    'MaxFunctionEvaluations',10^+20,...
    'OptimalityTolerance', 10^-20,...
    'StepTolerance',10^-20,...
    'MaxIterations',1000000);

%%% Random K-fold CV partitioning
Kfold_partitioning_time_begin = datetime(datestr(now));
[cv_part] = CVpart_random_Kfold(Ground,10);

%%% CV penalized likelihood
PenLik_time_begin = datetime(datestr(now))
lambda_seq = Geom_Seq(0.00000001,0.5,100,1);
[CV_perf_metrics] = DSTEM_fHDGM_PenLik_CV(cv_part,obj_stem_model,...
    input_data,lambda_seq,beta_dest_flag,opts,[],[]);


End_time = datetime(datestr(now))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Output storage      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_output_flag == 1
    Save_output_time_begin = datetime(datestr(now));
    save(['Application_ENV_' num2str(day(today),'%02.f') , num2str(month(today),'%02.f') '.mat'])
end