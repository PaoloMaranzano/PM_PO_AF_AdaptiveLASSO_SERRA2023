%% Paper: Adaptive LASSO estimation for functional hidden dynamic geostatistical models
%% Journal: Stochastic Environmental Research and Risk Assessment
%% Date: July 2022


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%% Decomposition of the random effects variance %%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% %%%%%%%%%% Main program %%%%%%%%%% %%
clear
clc


%%% Upload real data from ARPA Lombardia network
addpath(genpath('../../../VisitingLUH2021'));
addpath(genpath('../../../DSTEM_Code_General/DSTEM_software'));
addpath(genpath('../../../DSTEM_Code_General/Matlab_auxfuns'));
addpath(genpath('../../../SPASTA2021/Data'));

if ~exist('Ground')
	paper = "CSDA2021";
    run('Reshape_to_HDGM.m');
end
clearvars -except crossval_step data data_fhdgm Ground log_transform poll standardization



%%      Input data for simulations     %%
DSTEM_str_ground = Ground;
DSTEM_str_ground.poll = {'NO2'};
n_sites = {15};
n_covs = 3;
% Correlated covariates
cov_corr_flag = 0;
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



%% Decomposition of the variance of Y %%%%%%%%%%
% Setup
cov_corr_flag = 0;
st_corr_flag = 1;
ones_numb = 4;
DSTEM_str_sim_setup.beta = [[1 1 1 1 0 0 0],repmat(repelem([1 0],[ones_numb n_basis/nvars-ones_numb]),1,n_covs)]';
beta_true = DSTEM_str_sim_setup.beta;
DSTEM_str_sim_setup.nan_rate = [];
DSTEM_str_sim_setup.nan_pattern_par = [];
DSTEM_str_sim_setup.sigma_eps = repelem(0,n_basis/nvars)';
% For loop
nsims = 1;
seq_theta_z = Geom_Seq(0.0001,150,50,0);
% seq_theta_z = 50;
% seq_g = [0.5:0.05:0.95, 0.98];
seq_g = 0.85;
% seq_vz = 0:0.20:3;
seq_vz = 1;
% Combination
[theta_z,g,vz] = meshgrid(seq_theta_z,seq_g,seq_vz);
mod_comb = [array2table(vz(:),'VariableNames',{'v_z'}), ...
    array2table(g(:),'VariableNames',{'g'}), ...
    array2table(theta_z(:),'VariableNames',{'theta_z'})];
clear i m Sims_dec Decomp
for i = 1:nsims
    % Parameters
    for m = 1:size(mod_comb,1)
        DSTEM_str_sim_setup.theta_z = repelem(km2deg(mod_comb.theta_z(m)),n_basis/nvars);
        DSTEM_str_sim_setup.G = diag(repelem(mod_comb.g(m),n_basis/nvars));
        DSTEM_str_sim_setup.v_z = diag(repelem(mod_comb.v_z(m),n_basis/nvars));
        % Simulate
        [DSTEM_obj_sim,obj_stem_par,DSTEM_obj_sim_str,sim_ground] = DSTEM_fHDGM_sim(...
            DSTEM_str_ground,DSTEM_str_sim_setup,...
            n_covs,n_sites,X_names,X_varcov,datestamp_begin, datestamp_end, model_type,0,...
            fda_setup);
        % Decompose the variance
        Dec_temp = DSTEM_fHDGM_VarDecomp('ErrMat.csv','REMat.csv',DSTEM_obj_sim);
        Decomp(m,:) = [mod_comb(m,:) , Dec_temp];
    end
    Sims_dec{i} = Decomp;
end
%%% Effect of ST pars
1*1/(1-0.85^2)
% Theta_z
Theta_tab = vertcat(Sims_dec{:});
writetable(Theta_tab,'Sims_Theta.xlsx');
% G
G_tab = vertcat(Sims_dec{:});
writetable(G_tab,'Sims_G.xlsx');
% v_z
Vz_tab = vertcat(Sims_dec{:});
writetable(Vz_tab,'Sims_Vz.xlsx');
