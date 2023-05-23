%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%% DSTEM_fitted_given_beta %%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%%% Computes the fitted values of a HDGM given the vector of beta
%%% parameters

function [output] = DSTEM_fitted_given_beta(DSTEM_model,DSTEM_EM_opts,beta_vec)

if 0
    DSTEM_model = obj_temp_stem_model;
    DSTEM_EM_opts = obj_temp_stem_EM_options;
    beta_vec = beta_LQA(:,lam);
end

%%% Converting to stem_EM object
EM_st = stem_EM(DSTEM_model,DSTEM_EM_opts);

%%% Fitting E-step of EM algorithm to obtain the residuals
EM_st.stem_model.stem_par.beta = beta_vec;
[~,~,~,~,~,~,~,~,~,~,~,E_e_y1,~,~] = E_step(EM_st);

%%% Computing fitted values
y_hat = DSTEM_model.stem_data.Y;
y_hat(isnan(y_hat)) = 0;
y_hat = y_hat - E_e_y1;
%%% Computing residuals
res_hat = DSTEM_model.stem_data.Y - y_hat;

%%% Back transforming
blocks = [0 cumsum(DSTEM_model.dim)];
counter = 1;
for i=1:DSTEM_model.stem_data.stem_varset_p.nvar
    if DSTEM_model.stem_data.stem_varset_p.standardized
        s = DSTEM_model.stem_data.stem_varset_p.Y_stds{i};
        m = DSTEM_model.stem_data.stem_varset_p.Y_means{i};
    else
        s = 1;
        m = 0;
    end
    y_hat_back(blocks(counter)+1:blocks(counter+1),:) = y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m;
    y_back(blocks(counter)+1:blocks(counter+1),:) = DSTEM_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m;
    counter = counter+1;
end
res_hat_back = y_back - y_hat_back;

%%% Transform to matrix shape
n_st = size(DSTEM_model.stem_data.stem_varset_p.Y{1,1},1);
n_dd = size(DSTEM_model.stem_data.stem_varset_p.Y{1,1},2);
n_hh = length(DSTEM_model.stem_data.stem_varset_p.Y);

for h = 1:n_hh
    rows = (h-1)*n_st + 1:h*n_st;
    Y_hat{h,1} = y_hat(rows,:);
    R_hat{h,1} = res_hat(rows,:);
    Y_hat_back{h,1} = y_hat_back(rows,:);
    R_hat_back{h,1} = res_hat_back(rows,:);
    Y{h,1} = DSTEM_model.stem_data.Y(rows,:);
    Y_back{h,1} = y_back(rows,:);
end

for st = 1:n_st
    for dd = 1:n_dd
        for hh = 0:(n_hh-1)
            y1(hh+1,dd) = Y_hat{hh+1,1}(st,dd);
            r1(hh+1,dd) = R_hat{hh+1,1}(st,dd);
            y2(hh+1,dd) = Y_hat_back{hh+1,1}(st,dd);
            r2(hh+1,dd) = R_hat_back{hh+1,1}(st,dd);
            y3(hh+1,dd) = Y{hh+1,1}(st,dd);
            y4(hh+1,dd) = Y_back{hh+1,1}(st,dd);
        end
    end
    yhat(:,st) = y1(:);
    rhat(:,st) = r1(:);
    yhatback(:,st) = y2(:);
    rhatback(:,st) = r2(:);
    y(:,st) = y3(:);
    yback(:,st) = y4(:);
end

%%% Output
% Blocks shape
output.y_hat = y_hat;
output.res_hat = res_hat;
output.y_hat_back = y_hat_back;
output.res_hat_back = res_hat_back;
output.y = DSTEM_model.stem_data.Y;
output.y_back = y_back;
% Matrix shape
output.y_hat_mat = yhat;
output.res_hat_mat = yhat;
output.y_hat_back_mat = yhatback;
output.res_hat_back_mat = rhatback;
output.y_mat = y;
output.y_back_mat = yback;

end