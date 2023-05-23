function [] = DSTEM_HDGM_CV_PenLik_OptRes_Plot(...
    CV_perf_metrics,beta_true,fda_setup_true,fda_setup_estim,...
    vars_names,plot_list,...
    flag_1se,flag_1s,flag_lambda_large,...
    var_name_boxplot,boxplot_title,...
    sup_text_tit,...
    save_plot_flag,save_plot_path,save_plot_name)

% % fda_setup_true = fda_setup;
% % fda_setup_estim = input_fda;
% flag_1se = 1;
% flag_1s = 1;
% flag_lambda_large = 1;
% save_plot_flag = 0;
% % vars_names = {'Intercept','X1','X2','X3'};
% plot_list = {'A','B','C','E','F'};
% CV_perf_metrics = optimal_vals;

if isempty(beta_true) && isempty(fda_setup_true)
    fda_setup = fda_setup_estim;
    simulated = 0;
else
    fda_setup = fda_setup_true;
    beta_true = beta_true;
    simulated = 1;
end
input_fda = fda_setup_estim;
zoom = 30;

%%%%% Colors setting
col_sup_title = hex2rgb('#000000');   % Black
col_central = hex2rgb('#0033FF');
col_CI = hex2rgb('#FF0000');
col_min_MAE = hex2rgb('#999999');
col_1se_MAE = hex2rgb('#FF33FF');
col_1s_MAE = hex2rgb('#33CC00');
col_min_RMSE = hex2rgb('#339900');
col_1se_RMSE = hex2rgb('#FF6600');
col_1s_RMSE = hex2rgb('#9900FF');
col_MLE = hex2rgb('#000000');
col_sim = hex2rgb('#FF0000');
col_lambda_large = hex2rgb('#FF99FF');



%%% log(lambda)
lambda_seq = CV_perf_metrics.Lambda_seq;
log_lambda = log(lambda_seq);
dist = log_lambda(2) - log_lambda(3);
log_lambda(1) = log_lambda(2) - abs(dist);


%%%%%%%%%%%%%%%%%
%%%%% Plots %%%%%
%%%%%%%%%%%%%%%%%
tobs = fda_setup.spline_range(1):0.01:fda_setup.spline_range(2);
basis_obj = create_bspline_basis(fda_setup.spline_range, fda_setup.n_basis,...
    fda_setup.spline_order, fda_setup.spline_knots);
basis_obj_mat = eval_basis(tobs, basis_obj);
basis_mat = full(basis_obj_mat);



%% %%%%% Plot A: K-fold CV MAE and RMSE vs log(lambda)
% Useful quantities
MAE = CV_perf_metrics.Avg_before_opt.avgnsims_MAE;
RMSE = CV_perf_metrics.Avg_before_opt.avgnsims_RMSE;
p_MAE = CV_perf_metrics.Avg_before_opt.Optimal_vals.pos_min_MAE;
p_1se_MAE = CV_perf_metrics.Avg_before_opt.Optimal_vals.pos_oneSE_min_MAE;
p_1s_MAE = CV_perf_metrics.Avg_before_opt.Optimal_vals.pos_oneS_min_MAE;
p_RMSE = CV_perf_metrics.Avg_before_opt.Optimal_vals.pos_min_RMSE;
p_1se_RMSE = CV_perf_metrics.Avg_before_opt.Optimal_vals.pos_oneSE_min_RMSE;
p_1s_RMSE = CV_perf_metrics.Avg_before_opt.Optimal_vals.pos_oneS_min_RMSE;
m_MAE = nanmean(MAE,2);
s_MAE = sqrt(nanvar(MAE,[],2));
n_MAE = size(MAE,2);
se_MAE = sqrt(nanvar(MAE,[],2)) ./ sqrt(n_MAE);
m_RMSE = nanmean(RMSE,2);
s_RMSE = sqrt(nanvar(RMSE,[],2));
n_RMSE = size(RMSE,2);
se_RMSE = sqrt(nanvar(RMSE,[],2)) ./ sqrt(n_RMSE);
if any(contains(plot_list,{'A'}))
    f1 = figure('Name','plot A');
    f1.WindowState = 'maximized';
    t = sgtitle(sup_text_tit,'Color',col_sup_title,...
        'fontweight','bold','fontsize',20);
    %
    subplot(2,2,1)
    plot(log_lambda,m_MAE,'linewidth',2)
    hold on
    if flag_1se == 1
        xline(log_lambda(p_1se_MAE),'',lambda_seq(p_1se_MAE),'LineWidth',2,...
            'Color',col_1se_MAE);
        yline(m_MAE(p_1se_MAE),'',m_MAE(p_1se_MAE),'LineWidth',2,...
            'Color',col_1se_MAE);
    end
    if flag_1s == 1
        xline(log_lambda(p_1s_MAE),'',lambda_seq(p_1s_MAE),'LineWidth',2,...
            'Color',col_1s_MAE);
        yline(m_MAE(p_1s_MAE),'',m_MAE(p_1s_MAE),'LineWidth',2,...
            'Color',col_1s_MAE);
    end
    xline(log_lambda(p_MAE),'',lambda_seq(p_MAE),'LineWidth',2,...
        'Color',col_min_MAE)
    yline(m_MAE(p_MAE),'',m_MAE(p_MAE),'LineWidth',2,...
        'Color',col_min_MAE)
    xlabel('log(\lambda)')
    ylabel('10-fold average MAE')
    title('Adaptive LASSO: 10-fold average MAE vs log({\lambda})')
    %
    subplot(2,2,2)
    plot(log_lambda,m_RMSE,'linewidth',2)
    hold on
    if flag_1se == 1
        xline(log_lambda(p_1se_RMSE),'',lambda_seq(p_1se_RMSE),'LineWidth',2,...
            'Color',col_1se_RMSE);
        yline(m_RMSE(p_1se_RMSE),'',m_RMSE(p_1se_RMSE),'LineWidth',2,...
            'Color',col_1se_RMSE);
    end
    if flag_1s == 1
        xline(log_lambda(p_1s_RMSE),'',lambda_seq(p_1s_RMSE),'LineWidth',2,...
            'Color',col_1s_RMSE);
        yline(m_RMSE(p_1s_RMSE),'',m_RMSE(p_1s_RMSE),'LineWidth',2,...
            'Color',col_1s_RMSE);
    end
    xline(log_lambda(p_RMSE),'',lambda_seq(p_RMSE),'LineWidth',2,'Color',col_min_RMSE)
    yline(m_RMSE(p_RMSE),'',m_RMSE(p_RMSE),'LineWidth',2,'Color',col_min_RMSE)
    xlabel('log(\lambda)')
    ylabel('10-fold average RMSE')
    title('Adaptive LASSO: 10-fold average RMSE vs log({\lambda})')
    %
    subplot(2,2,3)
    errorbar(log_lambda(1:zoom),m_MAE(1:zoom),se_MAE(1:zoom),'-s','MarkerSize',10,...
        'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold on
    xline(log_lambda(p_MAE),'',lambda_seq(p_MAE),'LineWidth',2,...
        'Color',col_min_MAE)
    yline(m_MAE(p_MAE),'',m_MAE(p_MAE),'LineWidth',2,...
        'Color',col_min_MAE)
    if flag_1se == 1
        if p_1se_MAE <= zoom
            xline(log_lambda(p_1se_MAE),'',lambda_seq(p_1se_MAE),'LineWidth',2,...
                'Color',col_1se_MAE)
            yline(m_MAE(p_1se_MAE),'',m_MAE(p_1se_MAE),'LineWidth',2,...
                'Color',col_1se_MAE)
        end
    end
    if flag_1s == 1
        if p_1s_MAE <= zoom
            xline(log_lambda(p_1s_MAE),'',lambda_seq(p_1s_MAE),'LineWidth',2,...
                'Color',col_1s_MAE)
            yline(m_MAE(p_1s_MAE),'',m_MAE(p_1s_MAE),'LineWidth',2,...
                'Color',col_1s_MAE)
        end
    end
    xlabel('log(\lambda)')
    ylabel('10-fold average MAE and error bars (std. error)')
    title('(Zoom) Adaptive LASSO: 10-fold average MAE vs log({\lambda})')
    %
    subplot(2,2,4)
    errorbar(log_lambda(1:zoom),m_RMSE(1:zoom),se_RMSE(1:zoom),'-s','MarkerSize',10,...
        'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold on
    xline(log_lambda(p_RMSE),'',lambda_seq(p_RMSE),'LineWidth',2,'Color',col_min_RMSE)
    yline(m_RMSE(p_RMSE),'',m_RMSE(p_RMSE),'LineWidth',2,'Color',col_min_RMSE)
    if flag_1se == 1
        if p_1se_RMSE <= zoom
            xline(log_lambda(p_1se_RMSE),'',lambda_seq(p_1se_RMSE),'LineWidth',2,...
                'Color',col_1se_RMSE)
            yline(m_RMSE(p_1se_RMSE),'',m_RMSE(p_1se_RMSE),'LineWidth',2,...
                'Color',col_1se_RMSE)
        end
    end
    if flag_1s == 1
        if p_1s_RMSE <= zoom
            xline(log_lambda(p_1s_RMSE),'',lambda_seq(p_1s_RMSE),'LineWidth',2,...
                'Color',col_1s_RMSE)
            yline(m_RMSE(p_1s_RMSE),'',m_RMSE(p_1s_RMSE),'LineWidth',2,...
                'Color',col_1s_RMSE)
        end
    end
    xlabel('log(\lambda)')
    ylabel('10-fold average RMSE and error bars (std. error)')
    title('(Zoom) Adaptive LASSO: 10-fold average RMSE vs log({\lambda})')
    set(gcf, 'PaperUnits', 'centimeters');
    x_width=35 ;y_width=25;
    set(gcf, 'PaperPosition', [30 25 x_width y_width]);
    if save_plot_flag == 1
        if isempty(save_plot_path) == 0
            saveas(gcf,[save_plot_path save_plot_name '_A.png'])
        else
            saveas(gcf,[save_plot_name '_A.png'])
        end
    end
end



%% %%%%% Plot B: Full sample spline plot
% Useful quantities
nvars = length(vars_names);
tobs = fda_setup.spline_range(1):0.01:fda_setup.spline_range(2);
beta_LASSO_LQA = CV_perf_metrics.Avg_before_opt.Full_sample.beta_LASSO_LQA;
l_large = length(log_lambda)-10;

if any(contains(plot_list,{'B'}))
    f2 = figure('Name','plot B');
    f2.WindowState = 'maximized';
    sup_t1 = sgtitle({[sup_text_tit] ['Functional \beta coefficients over the 24 hours']},...
        'Color',col_sup_title,'fontweight','bold','fontsize',20);
end
for i = 1:nvars
    if simulated == 1
        % Compute the theoretical curve (true) for i-th variable
        betas = beta_true;
        nbasis_theo = (fda_setup.knots_number + 2*fda_setup.spline_order - (fda_setup.spline_order+1));
        bas_theo = (i-1)*nbasis_theo + 1:i*nbasis_theo;
        [basis_obj_theo,basis_obj_mat_theo,basis_mat_theo,yhat_theo] = simulate_bspline_basis(...
            fda_setup,betas(bas_theo),[],[]);
        leg_txt{2} = {['Theoretical (k=5) [' num2str(beta_true(bas_theo))' ']']};
    end
    % MLE
    betas = beta_LASSO_LQA(:,1) ;
    nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
    bas = (i-1)*nbasis + 1:i*nbasis;
    [basis_obj,basis_obj_mat,basis_mat,yhat_MLE] = simulate_bspline_basis(...
        input_fda,betas(bas),[],[]);
    leg_txt{1} = {'MLE'};
    % MAE lambda optimal
    betas = beta_LASSO_LQA(:,p_MAE) ;
    nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
    bas = (i-1)*nbasis + 1:i*nbasis;
    [basis_obj,basis_obj_mat,basis_mat,yhat_lambda_star_MAE] = simulate_bspline_basis(...
        input_fda,betas(bas),[],[]);
    leg_txt{3} = 'PenLik at \lambda^* (MAE)';
    % 1-SE MAE lambda optimal
    if flag_1se == 1
        betas = beta_LASSO_LQA(:,p_1se_MAE) ;
        nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
        bas = (i-1)*nbasis + 1:i*nbasis;
        [basis_obj,basis_obj_mat,basis_mat,yhat_lambda_star_1se_MAE] = simulate_bspline_basis(...
            input_fda,betas(bas),[],[]);
    end
    % 1-S MAE lambda optimal
    if flag_1s == 1
        betas = beta_LASSO_LQA(:,p_1s_MAE) ;
        nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
        bas = (i-1)*nbasis + 1:i*nbasis;
        [basis_obj,basis_obj_mat,basis_mat,yhat_lambda_star_1s_MAE] = simulate_bspline_basis(...
            input_fda,betas(bas),[],[]);
    end
    % RMSE lambda optimal
    betas = beta_LASSO_LQA(:,p_RMSE) ;
    nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
    bas = (i-1)*nbasis + 1:i*nbasis;
    [basis_obj,basis_obj_mat,basis_mat,yhat_lambda_star_RMSE] = simulate_bspline_basis(...
        input_fda,betas(bas),[],[]);
    leg_txt{4} = 'PenLik at \lambda^* (RMSE)';
    % 1-SE RMSE lambda optimal
    if flag_1se == 1
        betas = beta_LASSO_LQA(:,p_1se_RMSE) ;
        nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
        bas = (i-1)*nbasis + 1:i*nbasis;
        [basis_obj,basis_obj_mat,basis_mat,yhat_lambda_star_1se_RMSE] = simulate_bspline_basis(...
            input_fda,betas(bas),[],[]);
        leg_txt{5} = 'PenLik at 1-SE \lambda^* (MAE)';
        leg_txt{6} = 'PenLik at 1-SE \lambda^* (RMSE)';
    end
    % 1-S RMSE lambda optimal
    if flag_1s == 1
        betas = beta_LASSO_LQA(:,p_1s_RMSE) ;
        nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
        bas = (i-1)*nbasis + 1:i*nbasis;
        [basis_obj,basis_obj_mat,basis_mat,yhat_lambda_star_1s_RMSE] = simulate_bspline_basis(...
            input_fda,betas(bas),[],[]);
        leg_txt{7} = 'PenLik at 1-S \lambda^* (MAE)';
        leg_txt{8} = 'PenLik at 1-S \lambda^* (RMSE)';
    end
    % Lambda large
    if flag_lambda_large == 1
        l_large = length(log_lambda)-10;
        betas = beta_LASSO_LQA(:,l_large) ;
        nbasis = input_fda.knots_number + 2*input_fda.spline_order - (input_fda.spline_order+1);
        bas = (i-1)*nbasis + 1:i*nbasis;
        [basis_obj,basis_obj_mat,basis_mat,yhat_lambda_large] = simulate_bspline_basis(...
            input_fda,betas(bas),[],[]);
        leg_txt{9} = 'PenLik at \lambda large';
    end
    if any(contains(plot_list,{'B'}))
        % Plot
        if simulated == 1
            subplot(2,3,i)
        else
            subplot(4,3,i)
        end
        h(:,i) = plot(tobs,yhat_MLE,'Color',col_MLE,'linewidth',2)
        leg_txt{1} = {'MLE'};
        hold on
        if simulated == 1
            plot(tobs,yhat_theo,'Color',col_sim,'linewidth',2)
        end
        plot(tobs,yhat_lambda_star_MAE,'Color',col_min_MAE,'linewidth',2)
        plot(tobs,yhat_lambda_star_RMSE,'Color',col_min_RMSE,'linewidth',2)
        if flag_1se == 1
            plot(tobs,yhat_lambda_star_1se_MAE,'Color',col_1se_MAE,'linewidth',2)
            plot(tobs,yhat_lambda_star_1se_RMSE,'Color',col_1se_RMSE,'linewidth',2)
        end
        if flag_1s == 1
            plot(tobs,yhat_lambda_star_1s_MAE,'Color',col_1s_MAE,'linewidth',2)
            plot(tobs,yhat_lambda_star_1s_RMSE,'Color',col_1s_RMSE,'linewidth',2)
        end
        if flag_lambda_large == 1
            plot(tobs,yhat_lambda_large,'Color',col_lambda_large,'linewidth',2)
        end
        title(vars_names{1,i})
    end
end
if any(contains(plot_list,{'B'}))
    % Common legend
    Lgnd = legend(string(leg_txt(find(~(cellfun('isempty',leg_txt))))));
    Lgnd = legend(legendUnq(Lgnd));
    Lgnd.FontSize = 12;
    x_width = 0.10 ; y_width = 0.21;
    x_height = 0.6 ; y_height = 0.1;
    % Lgnd.Position = [x_height y_height x_width y_width];
    Lgnd.Position(1:2) = [.5-Lgnd.Position(3)/2 .025];
    % Common labels
    [ax1,h1]=suplabel('Time (hours)','x');
    [ax2,h2]=suplabel('Full sample: adaptive LASSO \beta Coefficients (spline)','y');
    set(h1,'FontSize',14)
    set(h2,'FontSize',14)
    hold off
    set(gcf, 'PaperUnits', 'centimeters');
    x_width = 40 ; y_width = 30;
    x_height = 35 ; y_height = 30;
    set(gcf, 'PaperPosition', [x_height y_height x_width y_width]);
    if save_plot_flag == 1
        if isempty(save_plot_path) == 0
            saveas(gcf,[save_plot_path save_plot_name '_B.png'])
        else
            saveas(gcf,[save_plot_name '_B.png'])
        end
    end
end

%% %%%%% Plot C: beta plot
if any(contains(plot_list,{'C'}))
    cols = string.empty(size(beta_LASSO_LQA,1),0);
    type = string.empty(size(beta_LASSO_LQA,1),0);
    if simulated == 1
        for v = 1:size(beta_LASSO_LQA,1)
            if beta_true(v) ~= 0
                cols(v) = "b";
                type(v) = "-";
            else
                cols(v) = "r";
                type(v) = "-.";
            end
        end
    else
        for v = 1:size(beta_LASSO_LQA,1)
            cols(v) = "b";
            type(v) = "-";
        end
    end
    f3 = figure('Name','plot C');
    f3.WindowState = 'maximized';
    sup_t = title(sup_text_tit,'Color',col_sup_title,...
        'fontweight','bold','fontsize',20);
    t = subtitle('Full sample: Adaptive LASSO \beta coefficients vs log(\lambda)');
    set(t,'fontweight','bold','fontsize',16);
    n_basis_int = size(beta_LASSO_LQA,1) / nvars;
    for v = (n_basis_int+1):size(beta_LASSO_LQA,1)
        hold on
        plot(log_lambda,beta_LASSO_LQA(v,:),type(v),'Color',cols(v),...
            'HandleVisibility','off');
    end
    xlabel('log(\lambda)')
    ylabel('Full sample: adaptive LASSO \beta Coefficients (spline)')
    xline(log_lambda(1),'',lambda_seq(1),'Color',col_MLE,'LineWidth',2,...
        'DisplayName','MLE')
    xline(log_lambda(p_MAE),'',lambda_seq(p_MAE),'LineWidth',2,...
        'Color',col_min_MAE,'DisplayName','PenLik at \lambda^* (MAE)')
    xline(log_lambda(p_RMSE),'',lambda_seq(p_RMSE),'LineWidth',2,...
        'Color',col_min_RMSE,'DisplayName','PenLik at \lambda^* (RMSE)')
    if flag_1se == 1
        xline(log_lambda(p_1se_MAE),'',lambda_seq(p_1se_MAE),'LineWidth',2,...
            'Color',col_1se_MAE,'DisplayName','PenLik at 1-SE \lambda^* (MAE)')
        xline(log_lambda(p_1se_RMSE),'',lambda_seq(p_1se_RMSE),'LineWidth',2,...
            'Color',col_1se_RMSE,'DisplayName','PenLik at 1-SE \lambda^* (RMSE)')
    end
    if flag_1s == 1
        xline(log_lambda(p_1s_MAE),'',lambda_seq(p_1s_MAE),'LineWidth',2,...
            'Color',col_1s_MAE,'DisplayName','PenLik at 1-S \lambda^* (MAE)')
        xline(log_lambda(p_1s_RMSE),'',lambda_seq(p_1s_RMSE),'LineWidth',2,...
            'Color',col_1s_RMSE,'DisplayName','PenLik at 1-S \lambda^* (RMSE)')
    end
    if flag_lambda_large == 1
        xline(log_lambda(l_large),'',lambda_seq(l_large),...
            'Color',col_lambda_large,'LineWidth',2,...
            'DisplayName','PenLik at \lambda large')
    end
    legend(string(leg_txt(find(~(cellfun('isempty',leg_txt))))),...
        'Location', 'Best')
    set(gcf, 'PaperUnits', 'centimeters');
    x_width = 30 ; y_width = 20;
    x_height = 30 ; y_height = 25;
    set(gcf, 'PaperPosition', [x_height y_height x_width y_width]);
    if save_plot_flag == 1
        if isempty(save_plot_path) == 0
            saveas(gcf,[save_plot_path save_plot_name '_C.png'])
        else
            saveas(gcf,[save_plot_name '_C.png'])
        end
    end
end



%% %%%%% Plot D: (avg. across sims) RMSE by fold
if any(contains(plot_list,{'D'}))
    f4 = figure('Name','plot D');
    f4.WindowState = 'maximized';
    sup_t = sgtitle(sup_text_tit,'Color',col_sup_title,...
        'fontweight','bold','fontsize',20)
    subplot(1,2,1)
    plot(log_lambda,CV_perf_metrics.Avg_before_opt.avgnsims_RMSE)
    xlabel('log(\lambda)')
    ylabel('RMSE')
    set(gcf, 'PaperUnits', 'centimeters');
    x_width = 30 ; y_width = 20;
    x_height = 30 ; y_height = 25;
    set(gcf, 'PaperPosition', [x_height y_height x_width y_width]);
    t2 = title('RMSE by CV fold vs log(\lambda)');
    set(t2,'fontweight','bold','fontsize',16);
    %
    subplot(1,2,2)
    plot(log_lambda, CV_perf_metrics.Avg_before_opt.avgnsims_MAE)
    xlabel('log(\lambda)')
    ylabel('MAE')
    set(gcf, 'PaperUnits', 'centimeters');
    x_width = 30 ; y_width = 20;
    x_height = 30 ; y_height = 25;
    set(gcf, 'PaperPosition', [x_height y_height x_width y_width]);
    t3 = title('MAE by CV fold vs log(\lambda)');
    set(t3,'fontweight','bold','fontsize',16);
    if save_plot_flag == 1
        if isempty(save_plot_path) == 0
            saveas(gcf,[save_plot_path save_plot_name '_D.png'])
        else
            saveas(gcf,[save_plot_name '_D.png'])
        end
    end
end



%% %%%%% Plot E: (avg. across sims) RMSE by fold
if any(contains(plot_list,{'E'}))
    var_boxplot = CV_perf_metrics.Avg_after_opt.Full_sample.beta_LASSO_LQA.(var_name_boxplot);
    f5 = figure('Name','plot E');
    f5.WindowState = 'maximized';
    hAx=gca;
    boxplot(var_boxplot','Labels',CV_perf_metrics.Coef_names)
    hAx.XAxis.TickLabelInterpreter='tex';
    hold on
    yline(1,':k','True value')
    yline(1,':k','\beta=1','LabelHorizontalAlignment','right',...
        'LabelVerticalAlignment','bottom')
    yline(0,':k','True value','LabelHorizontalAlignment','left')
    yline(0,':k','\beta=0','LabelHorizontalAlignment','left',...
        'LabelVerticalAlignment','bottom')
    annotation('textbox',[.2 0 0 .08],...
        'String','Intercept','EdgeColor','none','FontSize',16)
    annotation('rectangle',[.13 0.04 0.20 .08],'EdgeColor','black','LineWidth',2)
    annotation('textbox',[.4 0 0 .08],...
        'String','X_1','EdgeColor','none','FontSize',16)
    annotation('rectangle',[.33 0.04 0.19 .08],'EdgeColor','black','LineWidth',2)
    annotation('textbox',[.6 0 0 .08],...
        'String','X_2','EdgeColor','none','FontSize',16)
    annotation('rectangle',[.52 0.04 0.19 .08],'EdgeColor','black','LineWidth',2)
    annotation('textbox',[.8 0 0 .08],...
        'String','X_4','EdgeColor','none','FontSize',16)
    annotation('rectangle',[.71 0.04 0.19 .08],'EdgeColor','black','LineWidth',2)
    hold off
    set(gcf, 'PaperUnits', 'centimeters');
    sup_t = title(sup_text_tit);
    set(sup_t,'Color','black',...
        'fontweight','bold','fontsize',20);
%     set(sup_t,'Position',[0.53 -0.08 0],'Color','black',...
%         'fontweight','bold','fontsize',20);
    t = subtitle(boxplot_title);
    set(t,'fontweight','bold','fontsize',16);
    x_width = 50 ; y_width = 30;
    x_height = 35 ; y_height = 30;
    set(gcf, 'PaperPosition', [x_height y_height x_width y_width]);
    if save_plot_flag == 1
        if isempty(save_plot_path) == 0
            saveas(gcf,[save_plot_path save_plot_name '_E_' var_name_boxplot '.png'])
        else
            saveas(gcf,[save_plot_name '_E_' var_name_boxplot '.png'])
        end
    end
end



%% %%%%% Figure F: RMSE and MAE by simulation
if any(contains(plot_list,{'F'}))
    m_RMSE = CV_perf_metrics.Avg_after_opt.avgKfold_RMSE;
    m_MAE = CV_perf_metrics.Avg_after_opt.avgKfold_MAE;
    alpha = 0.05;
    %%% RMSE
    central_m_RMSE = mean(m_RMSE,2);
    s_m_RMSE = sqrt(nanvar(m_RMSE,[],2));
    se_m_RMSE = s_m_RMSE / size(m_RMSE,2);
    upper_m_RMSE = central_m_RMSE + tinv(1-alpha/2,length(central_m_RMSE)-1)*se_m_RMSE;
    lower_m_RMSE = central_m_RMSE - tinv(1-alpha/2,length(central_m_RMSE)-1)*se_m_RMSE;
    %%% MAE
    central_m_MAE = mean(m_MAE,2);
    s_m_MAE = sqrt(nanvar(m_MAE,[],2));
    se_m_MAE = s_m_MAE / size(m_MAE,2);
    upper_m_MAE = central_m_MAE + tinv(1-alpha/2,length(central_m_MAE)-1)*se_m_MAE;
    lower_m_MAE = central_m_MAE - tinv(1-alpha/2,length(central_m_MAE)-1)*se_m_MAE;
    %%% Plot
    f6 = figure('Name','plot F');
    f6.WindowState = 'maximized';
    sup_t = sgtitle(sup_text_tit,'Color',col_sup_title,...
        'fontweight','bold','fontsize',20)
    % RMSE
    subplot(1,2,1)
    for run = 1:size(m_RMSE,2)
        hold on
        plot(log_lambda,m_RMSE(:,run),'black')
    end
    p2 = plot(log_lambda,central_m_RMSE,'','LineWidth',3,...
        'DisplayName','Average curve','Color',col_min_RMSE);
    p3 = plot(log_lambda,upper_m_RMSE,'',...
        'LineWidth',2,'LineStyle','--',...
        'DisplayName','95% CI','Color',col_min_RMSE);
    p4 = plot(log_lambda,lower_m_RMSE,'',...
        'LineWidth',2,'LineStyle','--',...
        'DisplayName','95% CI','Color',col_min_RMSE);
    xlabel('log(\lambda)')
    ylabel('Average 10-fold RMSE')
    t = title(['Average 10-fold RMSE across ' num2str(size(m_RMSE,2)) ' simulations']);
    set(t,'fontweight','bold','fontsize',16);
    legend([p2 p3],'Location','southeast')
    x_width = 50 ; y_width = 30;
    x_height = 35 ; y_height = 30;
    set(gcf, 'PaperPosition', [x_height y_height x_width y_width]);
    % MAE
    subplot(1,2,2)
    for run = 1:size(m_MAE,2)
        hold on
        plot(log_lambda,m_MAE(:,run),'black')
    end
    p2 = plot(log_lambda,central_m_MAE,'','LineWidth',3,...
        'DisplayName','Average curve','Color',col_min_MAE);
    p3 = plot(log_lambda,upper_m_MAE,'',...
        'LineWidth',2,'LineStyle','--',...
        'DisplayName','95% CI','Color',col_min_MAE);
    p4 = plot(log_lambda,lower_m_MAE,'',...
        'LineWidth',2,'LineStyle','--',...
        'DisplayName','95% CI','Color',col_min_MAE);
    xlabel('log(\lambda)')
    ylabel('Average 10-fold MAE')
    t = title(['Average 10-fold MAE across ' num2str(size(m_MAE,2)) ' simulations']);
    set(t,'fontweight','bold','fontsize',16);
    legend([p2 p3],'Location','southeast')
    x_width = 50 ; y_width = 30;
    x_height = 35 ; y_height = 30;
    set(gcf, 'PaperPosition', [x_height y_height x_width y_width]);
    if save_plot_flag == 1
        if isempty(save_plot_path) == 0
            saveas(gcf,[save_plot_path save_plot_name '_F.png'])
        else
            saveas(gcf,[save_plot_name '_F.png'])
        end
    end
end


end