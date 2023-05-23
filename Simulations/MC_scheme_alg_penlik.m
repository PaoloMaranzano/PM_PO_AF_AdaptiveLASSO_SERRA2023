%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%% Monte Carlo simulation of penalized   %%%%%%%%%% %%
%% %%%%%%%%%%  f-HDGM using adaptive LSSSO penalty  %%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% %%%%%%%%%% Simulation and feature selection algorithm (penalized likelihood) %%%%%%%%%% %%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Input data for simulations     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DSTEM_str_ground = Ground;
DSTEM_str_ground.poll = {'NO2'};
n_sites = {15};
n_covs = 3;
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
DSTEM_str_sim_setup.sigma_eps = repelem(0,n_basis/nvars)';		% Mod on 09/02/2023
% Spatio-temporal correlation
if st_corr_flag == 1
    st_msg = 'STcorr';
	DSTEM_str_sim_setup.v_z = diag(repelem(1,n_basis/nvars));
	DSTEM_str_sim_setup.theta_z = repelem(km2deg(50),n_basis/nvars);
	DSTEM_str_sim_setup.G = diag(repelem(0.85,n_basis/nvars));
else
    st_msg = 'STuncorr';
	DSTEM_str_sim_setup.v_z = diag(repelem(0,n_basis/nvars));
	DSTEM_str_sim_setup.theta_z = repelem(km2deg(0.00001),n_basis/nvars);
	DSTEM_str_sim_setup.G = diag(repelem(0,n_basis/nvars));
end


for rep = [(File_number-1)*MC_reps + 1:File_number*MC_reps]
    
	%%%%%%%%%%%%%%%%%%%%%%
    %%      Message     %%
    %%%%%%%%%%%%%%%%%%%%%%
	sprintf(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',...
            'Running simulation %d of %d',...
            ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%'],rep,File_number*MC_reps)
	
	
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Simulation     %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
	rng(rep + MC_cum)
    [DSTEM_obj_sim,obj_stem_par,DSTEM_obj_sim_str] = DSTEM_fHDGM_sim(...
        DSTEM_str_ground,DSTEM_str_sim_setup,...
        n_covs,n_sites,X_names,X_varcov,datestamp_begin, datestamp_end, model_type,0,...
        fda_setup);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Model estimation      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Setup inputs
    % FDA inputs
    input_fda.spline_type = 'Bspline';      % Spline type
    input_fda.spline_order = 3;
    input_fda.knots_number = 5;
    input_fda.spline_range = [0 24];        % Spline basis domain
    input_fda.spline_knots = linspace(input_fda.spline_range(1),input_fda.spline_range(2),input_fda.knots_number);
    obj_stem_fda = stem_fda(input_fda);
    n_basis = (input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1))*nvars;
    % Extracting input data from the simulation setting
    input_data = DSTEM_obj_sim_str.input_data;
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
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Data transform      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if standardization_flag == 1
        obj_stem_model.stem_data.standardize;
    end
    
    %%% Starting values of EM algorithm
    obj_stem_par.beta = obj_stem_model.get_beta0();
	obj_stem_par.sigma_eps = repelem(0,n_basis/nvars)';
	if st_corr_flag == 1
		obj_stem_par.G = diag(repelem(0.85,n_basis/nvars));
		obj_stem_par.theta_z = repelem(0.45,n_basis/nvars);
		obj_stem_par.v_z = diag(repelem(1,n_basis/nvars));
	else
	    obj_stem_par.G = diag(repelem(0.01,n_basis/nvars));
		obj_stem_par.theta_z = repelem(0.01,n_basis/nvars);
		obj_stem_par.v_z = diag(repelem(0.01,n_basis/nvars));
	end
	obj_stem_model.set_initial_values(obj_stem_par);
    
    %%% Parameters estimation
    date_estimation_begin = datetime('now')
    obj_stem_EM_options = stem_EM_options();
    obj_stem_EM_options.exit_tol_loglike = 0.0001;
    obj_stem_EM_options.exit_tol_par = 0.0001;
    obj_stem_EM_options.max_iterations = 200;
    obj_stem_model.EM_estimate(obj_stem_EM_options);
    
    %%% Var-cov matrix of estimated pars
    if varcov_estim_flag == 1
        date_begin_varcov = datetime('now')
		if ExactVarcov == 1
		    % Exact
			obj_stem_model.set_varcov;
		else
			% Approximated
			delta_varcov = 0.001;
			obj_stem_model.set_varcov(delta_varcov);
		end
    end
    
    %%% Log-likelihood of the data
    if logL_estim_flag == 1
        date_begin_logL = datetime('now')
        obj_stem_model.set_logL;
        logL_estim = obj_stem_model.stem_EM_result.logL
    end
	    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Model selection using penalized likelihood      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fisher algorithm settings
    % opts.betascale = 1;
    % opts.thresh = 1e-4;
    % opts.lambdaminratio = 1e-4;
    % opts.debug = 1;
    % opts.maxnewvars = 5;
    % devchange = Improvement of the penalized log-likelihood
    % opts.devchange  = 1e-4;
    % trustrgn = Omega weight of Levenberg-Marquardt
    % opts.trustrgn = 0.1;
    % trustiter = maximum number of iterations of Levember-Marquardt (internal loop)
    % opts.trustiter = 5000;
    % coreiter = maximum number of iterations of Fisher (external loop)
    % opts.coreiter = 5000;
	if beta_dest_flag == 1
        [MLE] = DSTEM_extraction(obj_stem_model,[],[]);
        beta_MLE = MLE.Reg_pars.Beta_destd_tab.Y.Coef;
    else
        beta_MLE = obj_stem_model.stem_EM_result.stem_par.beta;
    end
    if beta0_pen_flag == 1
        opts.penaltywt = ones(length(beta_MLE),1);
    else
        opts.penaltywt = [zeros(n_basis/nvars,1) ; ones(n_basis/nvars*n_covs,1)];
    end
    % betathresh = scale factor of the improvement of the penalized beta coefficients in LM
    % opts.betathresh = 1e-4;
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
    [cv_part] = CVpart_random_Kfold(DSTEM_obj_sim_str,10);
    
    %%% CV penalized likelihood
    % lambda_seq = [Geom_Seq(0.00001,0.001,50,1),Geom_Seq(0.0011,0.5,30,0)];
    % lambda_seq = [Geom_Seq(10^-6,10^-3.5,100,1)];
    % lambda_seq = [Geom_Seq(0.1,10,100,1)];
    lambda_seq = Geom_Seq(0.00000001,0.5,100,1);
    [CV_perf_metrics] = DSTEM_fHDGM_PenLik_CV(cv_part,obj_stem_model,...
        input_data,lambda_seq,beta_dest_flag,opts,[],[]);
    
    %%% Save MC rep CV performance metrics
    idx_rep = rep - (File_number-1)*MC_reps;
    % CV performance metrics
	MC_CV_perf_met{idx_rep,1} = CV_perf_metrics;
	% ML estimates extraction
	[MLE_fullsample,Perform_metrics_fullsample,Setup_fullsample] = DSTEM_extraction(obj_stem_model,[],[]);
	Est_mod_MLE.MLE = MLE_fullsample;
	Est_mod_MLE.Perform_metrics = Perform_metrics_fullsample;
	Est_mod_MLE.Setup = Setup_fullsample;
	Est_mod_MLE.Beta_true = DSTEM_str_sim_setup.beta;
	FullSample_MLE{idx_rep,1} = Est_mod_MLE;
end

% Standardized data output name
if standardization_flag == 1
    std_msg = 'std1';
else
    std_msg = 'std0';
end
% De-standardized coefficients and VarCov
if beta_dest_flag == 1
    betadestd_msg = 'betadest1';
else
    betadestd_msg = 'betadest0';
end
% Penalize the intercept
if beta0_pen_flag == 1
    beta0pen_msg = 'beta0pen';
else
    beta0pen_msg = 'beta0unpen';
end

%%% Save output
if save_output_flag == 1
		save(['MCrep_' num2str(File_number) '_'  st_msg cov_msg '_' std_msg '_' ...
            betadestd_msg '_' beta0pen_msg '_' ...
			num2str(day(today),'%02.f') , num2str(month(today),'%02.f') '.mat'],...
			'MC_CV_perf_met','fda_setup','input_fda','beta_true')
end



