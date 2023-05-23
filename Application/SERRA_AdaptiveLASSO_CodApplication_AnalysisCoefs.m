%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Application of PenLik to functional HDGM model %% %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% Change folder
cd('C:/Users/paulm/Google Drive/QA&COVID19/VisitingLUH2021/Code/funHDGM_sim/Application')

% Auxiliary functions and data
addpath(genpath('../../../../VisitingLUH2021'));
% addpath(genpath('../../../SPASTA2021/Code/DSTEM_software'));
% addpath(genpath('../../../SPASTA2021/Code/Matlab_auxfuns'));
addpath(genpath('../../../DSTEM_Code_General/DSTEM_software'));
addpath(genpath('../../../DSTEM_Code_General/Matlab_auxfuns'));
addpath(genpath('../../../../SPASTA2021/Data'));

% Load interesting coefficients
nb5 = readtable('VarCovApprox_nbasis5_SpPart2\Coefficients.csv');
nb7 = readtable('VarCovApprox_nbasis7_SpPart2\Coefficients.csv');
nb9 = readtable('VarCovApprox_nbasis9_SpPart2\Coefficients.csv');

% Treshold
tresh = 10^-04;

perc_MLE = [abs(sum(nb5.beta_MLE) < tresh) / size(nb5,1) *100 , ...
    sum(abs(nb7.beta_MLE) < tresh) / size(nb7,1) *100 , ...
    sum(abs(nb9.beta_MLE) < tresh) / size(nb9,1) * 100];
perc_minMAE = [sum(abs(nb5.beta_min_MAE) < tresh) / size(nb5,1) *100 , ...
    sum(abs(nb7.beta_min_MAE) < tresh) / size(nb7,1) *100 , ...
    sum(abs(nb9.beta_min_MAE) < tresh) / size(nb9,1) * 100];
perc_1seMAE = [sum(abs(nb5.beta_1se_min_MAE) < tresh) / size(nb5,1) *100 , ...
    sum(abs(nb7.beta_1se_min_MAE) < tresh) / size(nb7,1) *100 , ...
    sum(abs(nb9.beta_1se_min_MAE) < tresh) / size(nb9,1) * 100];
perc_minRMSE = [sum(abs(nb5.beta_min_RMSE) < tresh) / size(nb5,1) *100 , ...
    sum(abs(nb7.beta_min_RMSE) < tresh) / size(nb7,1) *100 , ...
    sum(abs(nb9.beta_min_RMSE) < tresh) / size(nb9,1) * 100];
perc_1seRMSE = [sum(abs(nb5.beta_1se_min_RMSE) < tresh) / size(nb5,1) *100 , ...
    sum(abs(nb7.beta_1se_min_RMSE) < tresh) / size(nb7,1) *100 , ...
    sum(abs(nb9.beta_1se_min_RMSE) < tresh) / size(nb9,1) * 100];

perc = [ perc_MLE ; perc_minMAE ; perc_1seMAE ; perc_minRMSE ; perc_1seRMSE ];
perc = array2table(perc);
perc.Properties.VariableNames = {'p=5','p=7','p=9'};
perc.Properties.RowNames = {'MLE','minMAE','1seminMAE','minRMSE','1seminRMSE'};

perc
writetable(perc,'PercZeros_VarCovApprox_SpPart2_10_4.xlsx')




plot(nb5.beta_min_MAE)

% Togliere costante dai grafici 
% Mettere tabella con i tempi, togliere lambda e mettere RMSE e MAE 1-se


nb5_red = nb5(:,{'Variable','Beta','beta_min_RMSE','beta_1se_min_RMSE'});
nb7_red = nb7(:,{'Variable','Beta','beta_min_RMSE','beta_1se_min_RMSE'});
nb9_red = nb9(:,{'Variable','Beta','beta_min_RMSE','beta_1se_min_RMSE'});

coefs = outerjoin(nb5_red,nb7_red,"Keys",{'Variable','Beta'},"MergeKeys",1);
coefs = outerjoin(coefs,nb9_red,"Keys",{'Variable','Beta'},"MergeKeys",1);

coefs.beta_min_RMSE_nb5_red(abs(coefs.beta_min_RMSE_nb5_red) < tresh) = 0;
coefs.beta_1se_min_RMSE_nb5_red(abs(coefs.beta_1se_min_RMSE_nb5_red) < tresh) = 0;
coefs.beta_min_RMSE_nb7_red(abs(coefs.beta_min_RMSE_nb7_red) < tresh) = 0;
coefs.beta_1se_min_RMSE_nb7_red(abs(coefs.beta_1se_min_RMSE_nb7_red) < tresh) = 0;
coefs.beta_min_RMSE(abs(coefs.beta_min_RMSE) < tresh) = 0;
coefs.beta_1se_min_RMSE(abs(coefs.beta_1se_min_RMSE) < tresh) = 0;

writetable(coefs,'nullcoefs.xlsx');







