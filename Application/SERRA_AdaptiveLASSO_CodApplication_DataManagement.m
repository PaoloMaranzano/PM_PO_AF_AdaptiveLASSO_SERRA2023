%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Data management for application of 'Adaptive LASSO estimation for f-HDMG' %% %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% %%%%% Set working directory
% cd('C:/Users/paulm/Google Drive/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application')

%% %%%%% Auxiliary functions
addpath(genpath('../../../../VisitingLUH2021'));
% addpath(genpath('../../../SPASTA2021/Code/DSTEM_software'));
% addpath(genpath('../../../SPASTA2021/Code/Matlab_auxfuns'));
addpath(genpath('../../../DSTEM_Code_General/DSTEM_software'));
addpath(genpath('../../../DSTEM_Code_General/Matlab_auxfuns'));
addpath(genpath('../../../../SPASTA2021/Data'));

if ~exist('Ground')
    paper = "CSDA2021";
    run('Reshape_to_HDGM.m');
    run('Reshape_to_funHDGM.m');
end
clearvars -except crossval_step data data_backup data_fhdgm data_fhdgm_backup Ground Ground_backup log_transform poll standardization
clearvars -except crossval_step data data_fhdgm Ground log_transform poll standardization

%%% Sub-period selection: 01/03/2020 - 01/06/2020
d_start = datetime(2020,03,01,'TimeZone', 'Z');
d_end = datetime(2020,06,01,'TimeZone', 'Z');

%%% Filtering f-HDGM format according to the subperiod
idx_date = Ground.data_fhdgm.Time <= d_end & Ground.data_fhdgm.Time >= d_start;
Ground.data_fhdgm = Ground.data_fhdgm(idx_date,:);
idx_date = Ground.data_long.Date < d_end & Ground.data_long.Date >= d_start;
Ground.data_long = Ground.data_long(idx_date,:);

%%% Variable selection
vars_names = {'Temperature','Rainfall','WindU','WindV','RelHumid',...
    'Pressure','VegetationHigh','VegetationLow','GeopotHeight'};
for i = 1:length(vars_names)
    vars_names{i} = ['X_beta_' vars_names{i}];
end
Ground.data_fhdgm = Ground.data_fhdgm(:,{'Profile','Y_namelvl','Ylvl','X_h_Hour',...
    'X_beta_Cons',vars_names{:},'Y_coordinate','X_coordinate','Time'});
Ground.vars_names = vars_names;

%%% Adding extra information to f-HDGM format
% New profile index
Ground.data_fhdgm.Profile = [1:size(Ground.data_fhdgm,1)]';
% Change names for DSTEM
Ground.data_fhdgm.Properties.VariableNames{'Y_namelvl'} = 'Y_name';
Ground.data_fhdgm.Properties.VariableNames{'Ylvl'} = 'Y';

%%% Index variable for reshaping long to functional shape (after filtering)
[~,~,gId] = unique(Ground.data_long(:,{'y','m','d','IDStat'}), 'rows');
Ground.data_long.gId = gId;
Ground.data_long = Ground.data_long;

%%% Filtering original HDGM format according to the subperiod
t = datetime(Ground.date_time, 'ConvertFrom', 'datenum',...
    'Format', 'yyyy-MM-dd HH:mm','TimeZone', 'Z');
idx_date = t < d_end & t >= d_start;
Ground.NO2 = Ground.NO2(:,find(idx_date));
Ground.X_NO2 = Ground.X_NO2(:,:,find(idx_date));

%%% Reorder Y and X matrices according to latitude (DSTEM format)
% Passaggio aggiunto 25/03/2022 e da sistemare a monte in reshape_HDGM.m
Ground.ARPA_stats_reg = sortrows(Ground.ARPA_stats_reg,'Latitude');
ord = str2double(extractAfter(string(Ground.ARPA_stats_reg.New_cod_stz),'_'));
Ground.NO2 = Ground.NO2(ord,:);
Ground.X_NO2 = Ground.X_NO2(ord,:,:);

%%% Remove duplicate tables
clearvars data data_long data_fhdgm idx_date standardization t d_end d_start i gId vars_names;


if 1
    save('Application_data_ENV.mat','-v7.3')
end