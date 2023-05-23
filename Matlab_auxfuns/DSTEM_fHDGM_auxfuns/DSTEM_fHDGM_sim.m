%% %%%%%%%%%%%%%%%%%%%%% %%
%% %%%%% fHDGM_sim %%%%% %%
%% %%%%%%%%%%%%%%%%%%%%% %%

%%% Simulates a functional HDGM model (f-HDGM)

function [DSTEM_obj_sim,DSTEM_obj_sim_stem_par,...
    DSTEM_str_ground_sim, DSTEM_ground] = DSTEM_fHDGM_sim(DSTEM_str_ground,...
    DSTEM_str_sim_setup,n_covs,n_sites,X_names,X_varcov,...
    datestamp_begin, datestamp_end, model_type,spatialCV_step, fda_setup)

d1 = datetime(datestamp_begin,'InputFormat','dd-MM-yyyy HH:mm');
d2 = datetime(datestamp_end,'InputFormat','dd-MM-yyyy HH:mm');
time_stamps = d1:hours(1):d2;
T = length(time_stamps);
    
%%% Computation begin
date_begin = datetime('now')


%% %%%%% HDGM data building
poll = DSTEM_str_ground.poll;
n = cell(length(poll),1);
sites = cell(length(poll),1);
obj_stem_gridlist_p = stem_gridlist();
for p = 1:length(poll)
    %%% Sampling sites
    idx_st_poll = ismember(DSTEM_str_ground.ARPA_stats_reg.New_cod_stz,...
        DSTEM_str_ground.(['Coords_' DSTEM_str_ground.poll{p}]).IDStat);
    SOI = find(idx_st_poll);
    sites{p} = datasample(SOI,n_sites{1,p},'Replace',false);
    sites_idx{p} = find(ismember(SOI,sites{p}));    
    %%% Dependent variables
    ground.Y{p} = DSTEM_str_ground.([poll{p}])(sites_idx{p},1:T);
    ground.Y_name{p} = poll{p};
    n{p} = n_sites{1,p};
    %%% Loading covariates for the selected pollutants at each monitoring stations
    ground.X_beta_name{p} = ['Constant' X_names];
    ground.X_beta{p} = zeros(n{p},n_covs+1,T);
    for tt = 1:T
        ground.X_beta{p}(:,:,tt) = [repmat(1,n{p},1) , ...
            mvnrnd(zeros(n_covs,1),X_varcov{p},n{p})];
		% ground.X_beta{p}(:,:,tt) = [repmat(1,n{p},1) , ...
        %     mvnrnd(ones(n_covs,1),X_varcov{p},n{p})];
    end
    %%% Coordinates grid
    ground.coordinates{p} = [DSTEM_str_ground.(['Coords_' poll{p}]).Latitude(sites_idx{p}), ...
        DSTEM_str_ground.(['Coords_' poll{p}]).Longitude(sites_idx{p})];
    %%% Constant term
    ground.X_z{p} = ones(n{p}, 1);
    ground.X_z_name{p} = {['constant_' poll{p}]};
    %%% STEM_grid
    % contains all the information related to the sampling locations of
    % a single variable
    obj_stem_grid = cell(length(poll),1);
    obj_stem_grid{p} = stem_grid(ground.coordinates{p}, 'deg', 'sparse', 'point');
    %%% STEM_gridlist
    % collector of the stem_grid objects for all the variables;
    obj_stem_gridlist_p.add(obj_stem_grid{p});
    if spatialCV_step ~= 0
        %%% Cross-validation settings
        S_val{1,p} = 1:spatialCV_step:n{p};
        ground.S_val{1,p} = S_val{1,p};
    end
end



%% Transform data to functional shape
%%% Creating basic object/dataframe
out = num2cell(ground.X_beta{1}, [1 2]); %split into cell array keeping dimensions 1 and 2 together
out = vertcat(out{:}); %concatenate all the cells vertically
data_long = array2table(out);
data_long.Properties.VariableNames = ground.X_beta_name{p};
data_long.IDStat = repmat(sites_idx{p},T,1);
data_long.Date = repelem(time_stamps',n_sites{p},1);
c = DSTEM_str_ground.(['Coords_' poll{p}])(sites_idx{p},:);
c.IDStat = sites_idx{p};
data_long = outerjoin(data_long,c,'Keys','IDStat','MergeKeys',true,'type','left');
data_long.y = year(data_long.Date);
data_long.m = month(data_long.Date);
data_long.d = day(data_long.Date);
data_long.h = hour(data_long.Date);
% Necessary to order at the beginning as stem_misc.data_formatter has a
% bug in lines 1241-1245 and 1075 (See email from Philipp 18/06/2021)
data_long = sortrows(data_long,{'Latitude'},{'ascend'});
[~,~,gId] = unique( data_long(:,{'y','m','d','IDStat'}), 'rows');
data_long.gId = gId;

%%% Response variable
data_long.Y = repelem(0,size(data_long,1))';
Y = splitapply( @(x){x}, data_long.Y, gId );
Y = cellfun(@transpose,Y,'UniformOutput',false);
data_fhdgm = table(Y);
data_fhdgm.Y_name = repelem('Y',size(data_fhdgm,1))';
data_fhdgm.Y_name = cellstr(data_fhdgm.Y_name);

%%% Coordinates
% Latitude
Lati = splitapply( @(x){x}, data_long.Latitude, gId );
Lati = cellfun(@transpose,Lati,'UniformOutput',false);
data_fhdgm.Y_coordinate = cellfun(@nanmean,Lati,'UniformOutput',false);
data_fhdgm.Y_coordinate = cell2mat(data_fhdgm.Y_coordinate);
data_fhdgm.Properties.VariableUnits{'Y_coordinate'} = 'deg';
% Longitude
Longi = splitapply( @(x){x}, data_long.Longitude, gId );
Longi = cellfun(@transpose,Longi,'UniformOutput',false);
data_fhdgm.X_coordinate = cellfun(@nanmean,Longi,'UniformOutput',false);
data_fhdgm.X_coordinate = cell2mat(data_fhdgm.X_coordinate);
data_fhdgm.Properties.VariableUnits{'X_coordinate'} = 'deg';

%%% Date/Time
data_long.Date_day = data_long.Date;
data_long.Date_day.Format = 'dd-MM-yyyy';
d = splitapply( @(x){x}, data_long.Date_day, gId );
data_fhdgm.day = cellfun(@nanmean,d,'UniformOutput',false);
for i = 1:size(data_fhdgm,1)
    data_fhdgm.Time(i) = data_fhdgm.day{i};
end

%%% Profile
data_fhdgm.Profile = [1:1:size(data_fhdgm,1)]';

%%% gId (Grouping variable wrt hours of the day)
% gId_fun = splitapply( @(x){x}, data_long.h, gId );
% data_fhdgm.gId = cellfun(@transpose,gId_fun,'UniformOutput',false);

%%% Hour of the day
h = splitapply( @(x){x}, data_long.h, gId );
data_fhdgm.X_h_Hour = cellfun(@transpose,h,'UniformOutput',false);

%%% Constant
cons = splitapply( @(x){x}, data_long.Constant, gId );
data_fhdgm.X_beta_Cons = cellfun(@transpose,cons,'UniformOutput',false);

%%% Covariates
covs_fhdgm_names = strcat('X_beta_',X_names);
for i = 1:length(X_names)
    var = splitapply( @(x){x}, data_long{:,X_names{i}}, gId );
    var = cellfun(@transpose,var,'UniformOutput',false);
    data_fhdgm.X_beta_var = var;
    data_fhdgm.Properties.VariableNames{'X_beta_var'} = covs_fhdgm_names{i};
end

%%% Reorder columns
order = {'Profile','Y_name','Y','X_h_Hour','X_beta_Cons',...
    covs_fhdgm_names{:},...
    'Y_coordinate','X_coordinate','Time'};
[~, ind] = ismember(order, data_fhdgm.Properties.VariableNames);
data_fhdgm = data_fhdgm(:,order);



%% %%%%% Store original data information
ground.obj_stem_gridlist_p = obj_stem_gridlist_p;
ground.datestamp_begin = datestamp_begin;
ground.datestamp_end = datestamp_end;
DSTEM_ground = ground;



%% %%%%% STEM settings
%%% STEM_fda
obj_stem_fda = stem_fda(fda_setup);

%%% Model type
obj_stem_modeltype = stem_modeltype(model_type);

%%% create an object of class stem_data
input_data.stem_modeltype = obj_stem_modeltype;
input_data.data_table = data_fhdgm;
input_data.stem_fda = obj_stem_fda;
obj_stem_data = stem_data(input_data);

%%% DSTEM_par object creation
obj_stem_par_constraints = stem_par_constraints();
obj_stem_par_constraints.time_diagonal = 0;
obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);

%%% DSTEM_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);

%%% DSTEM_par export
DSTEM_obj_sim_stem_par = obj_stem_par;



%% %%%%% Simulation
% Create stem_sim object
DSTEM_obj_sim = stem_sim(obj_stem_model);
% Setup missing Y rate and pattern
DSTEM_obj_sim.nan_rate = DSTEM_str_sim_setup.nan_rate;
DSTEM_obj_sim.nan_pattern_par = DSTEM_str_sim_setup.nan_pattern_par;
% Setup the artificial parameters
DSTEM_obj_sim.stem_model.stem_par.beta = DSTEM_str_sim_setup.beta;
DSTEM_obj_sim.stem_model.stem_par_initial = DSTEM_obj_sim.stem_model.stem_par;
DSTEM_obj_sim.stem_model.stem_par_initial.theta_z = DSTEM_str_sim_setup.theta_z;
DSTEM_obj_sim.stem_model.stem_par_initial.v_z = DSTEM_str_sim_setup.v_z;
DSTEM_obj_sim.stem_model.stem_par_initial.sigma_eta = DSTEM_str_sim_setup.sigma_eta;
DSTEM_obj_sim.stem_model.stem_par_initial.G = DSTEM_str_sim_setup.G;
DSTEM_obj_sim.stem_model.stem_par_initial.sigma_eps = DSTEM_str_sim_setup.sigma_eps;
% Simulate Y
simulate(DSTEM_obj_sim,DSTEM_obj_sim.nan_rate,DSTEM_obj_sim.nan_pattern_par)


%% Save simulation settings/structure
% Save simulated Y values both in long and function form
[Y_long,Y_fun,Y_mat] = DSTEM_extract_Y_to_funY(DSTEM_obj_sim.stem_model,gId,0);
data_long.Y = Y_long;
input_data.data_table.Y = Y_fun;
DSTEM_ground.Y{1,p} = Y_mat;
input_data.data_long = data_long;

DSTEM_str_ground_sim.ARPA_stats_reg = DSTEM_str_ground.ARPA_stats_reg(unique(vertcat(sites{:,1})),:);
DSTEM_str_ground_sim.lon = DSTEM_str_ground.lon(unique(vertcat(sites{:,1})),:);
DSTEM_str_ground_sim.lat = DSTEM_str_ground.lat(unique(vertcat(sites{:,1})),:);
DSTEM_str_ground_sim.sites = size(DSTEM_str_ground_sim.ARPA_stats_reg,1);
DSTEM_str_ground_sim.Y_names = DSTEM_str_ground.poll;
DSTEM_str_ground_sim.X_names = X_names;
DSTEM_str_ground_sim.measure_type = 'simulated';
DSTEM_str_ground_sim.temporal_granularity = DSTEM_str_ground.temporal_granularity;
DSTEM_str_ground_sim.time_resampled = DSTEM_str_ground.time_resampled;
DSTEM_str_ground_sim.data_source = DSTEM_str_ground.data_source;
DSTEM_str_ground_sim.coordinate_unit = DSTEM_str_ground.coordinate_unit;
DSTEM_str_ground_sim.data_type = DSTEM_str_ground.data_type;
DSTEM_str_ground_sim.model_type = model_type;
d1 = datenum(d1);
d2 = datenum(d2);
DSTEM_str_ground_sim.date_time = d1:d2;
DSTEM_str_ground_sim.time_stamps = T;
for p = 1:length(poll)
    DSTEM_str_ground_sim.(['Coords_' poll{p}]) = DSTEM_str_ground.(['Coords_' poll{p}])(sites_idx{1,p},:);
    DSTEM_str_ground_sim.(poll{p}) = Y_mat;
    DSTEM_str_ground_sim.(['X_' poll{p}]) = ground.X_beta{p};
end

% Save further settings
DSTEM_str_ground_sim.fda_stem = obj_stem_fda;
DSTEM_str_ground_sim.input_data = input_data;
DSTEM_str_ground_sim.data_long = data_long;


return;

end