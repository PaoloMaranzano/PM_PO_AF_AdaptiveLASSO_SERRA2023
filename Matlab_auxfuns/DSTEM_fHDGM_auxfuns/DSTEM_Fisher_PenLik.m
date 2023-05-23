%% %%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%% fisher_DSTEM %%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%% %%

%%% Penalized maximum likelihood estimates of a generic HDGM using
%%% Newton-Raphson algorithm with Levenburg-Marquardt adjustment for the
%%% parameters. It allows to use several penalty functions.
%%% See PENALIZED toolbox (William McIlhagga, 2014)

function [beta,InfoCrit,fval,pen,Lgrad,flag] = DSTEM_Fisher_PenLik(...
    DSTEM_model,lambda,opts)

% FISHER(model,beta,lambda,alpha,opts) maximizes the penalized
% likelihood
%
%    logL(model,beta)/nobs-sum(wt.*penalty(lambda,beta))
%
% using a Fisher scoring + Levenburg-Marquardt trust-region algorithm, with a warm
% start at the given beta. fisher is an internal routine and should 
% not be called outside a wrapper.
%
% The algorithm deals with singularities in the penalty function at beta=0
% but ignores them every where else. 
%
% Inputs:
%   model  : an object specifying the model. 
%   beta   : the initial value of beta
%   lambda : current penalty weight
%   penalty: a penalty function 
%   opts   : an options structure 
%
% Outputs:
%   beta  : the best fitted beta
%   fval  : the function value at beta
%   pen   : the penalty value at beta
%   Lgrad : the gradient of the likelihood at beta
%   flag  : how the maximization ended

if 0
    lam = 10;
    lambda = lambda_seq(lam);
    DSTEM_model = obj_temp_stem_model;
end

penalty = opts.penalty;
n = DSTEM_model.stem_data.stem_varset_p.N * DSTEM_model.stem_data.stem_varset_p.T;
nobs = DSTEM_model.stem_data.stem_varset_p.N * DSTEM_model.stem_data.stem_varset_p.T;
trust = opts.trustrgn;          % Omega weight of Levenburg-Marquardt
wt = lambda*opts.penaltywt;
flag  = '';

beta = DSTEM_model.stem_EM_result.stem_par.beta;
npars = length(beta);
logL_MLE = DSTEM_model.stem_EM_result.logL;
varcov_MLE = DSTEM_model.stem_EM_result.stem_par.varcov(1:npars,1:npars);
F_MLE = inv(varcov_MLE)/nobs;

for iter=1:opts.coreiter
    
    % save beta for convergence calcs.
    oldbeta = beta;
    
    % calculate fval, and unpenalized gradient of model at beta
    pen  = sum( wt.*call_penalty(penalty,'',beta,lambda,opts) );
    if iter == 1
        pen  = sum( wt.*call_penalty(penalty,'',beta,lambda,opts) );
        logl = logL_MLE;
    else
        pen = pen;
        logl = logl;
    end
    fval = logl/nobs - pen;
    grad = repelem(0,length(oldbeta),1);
    
    % save the gradient of model for export - used in convergence
    % calculations
    Lgrad = grad;
    
    % find additions to the active set.
    inactive = find(beta==0);
    if ~isempty(inactive)
        potential = abs( checksubgrad(grad,penalty,beta,lambda,opts) );
        [p,idx] = sort(-potential);
        idx = idx(1:min(length(idx),opts.maxnewvars));
        idx(potential(idx)<=0)=[];
        % add infinitesimal change to beta over the new active entries
        beta(inactive(idx)) = eps*sign(grad(inactive(idx)));
    end
    
    % work out new active set
    active = beta~=0;
    if sum(active)==0
        % probably because lambda is too high
        flag = 'empty active';
        break;
    end
    
    % calculate the gradient over the active set
    grad(active) = grad(active) - wt(active).*restrictpenalty(penalty,'deriv',beta,lambda,opts);
    grad(~active)= 0;

    % compute the info matrix
    Info = F_MLE(find(active),find(active));
    
    % extract diagonal of Info; this will give the size of the trust region
    % ellipse
    DI = diag(diag(Info));
    
    % add in second derivatives of the penalty to the Info matrix.
    % Eq. 6: F/n + lambda*pi_cap
    p2 = restrictpenalty(penalty,'2ndderiv',beta,lambda,opts);
    if any(p2~=0)
        Info = Info + diag(wt(active).*p2);
    end
    
    % update beta using LM/trust region algorithm
    oldfval = fval;
    LM_step = zeros(size(beta));
    solveopts.SYM=true;
    
    for t=1:opts.trustiter
        % compute Levenberg-Marquardt step.
        % Solve the system eq. 7: F/n + lambda*pi_cap + trust*F/n
        % LM_step is the second term of Eq. 8
        newInfo = Info + trust*DI;
        [LM_step(active),rflag] = linsolve(newInfo,grad(active),solveopts);

        % project step onto allowable model region
        % Eq. 9: beta_t + LM_step for the active set
        b = beta + LM_step;
        
        % project onto the allowable orthant
        if penalty('project',beta,lambda,opts)
            bad = wt~=0 & sign(b)~=sign(beta);
            b(bad) = 0;
        end

        % calculate value of the penalized objective function at new point
        newpen = sum(wt.*call_penalty(penalty,'',b,lambda,opts));
        DSTEM_model.stem_par.beta = b;
        DSTEM_model.set_logL;
        newlogL = DSTEM_model.stem_EM_result.logL;
        newfval = newlogL/nobs - newpen;
        % Calculate the change in the penalized likelihood
        change  = newfval - fval;
        
        % calculate change from quadratic model. Note minus sign
        % because Info is negative of hessian.
        step = b(active) - beta(active);
        approx_change = step'*grad(active) - 0.5*(step'*Info*step);

        % adjust trust region
        if change>=0                    % some improvement
            rho = change/approx_change;
            if rho>0.9                  % good step, increase trust region
                trust = trust/2;
            elseif rho<0.1              % poor step, reduce trust region
                trust = trust*2;
            end
            % accept change & exit loop
            beta = b;
            fval = newfval;
            pen  = newpen;
            logl = newlogL;
            break;
        else                    % change<0
            trust = trust*64;   % very poor step, massively reduce trust region and try again
        end
    end
    coefs_iter_LQA(:,iter) = b;
    % Calculate information criteria
    [AIC,AICc,BIC,HQIC] = DSTEM_InfoCrit(DSTEM_model);
    InfoCrit = [AIC,AICc,BIC,HQIC];

    % converged when change in fval is small and there has been improvement
    convg = abs(oldfval-fval) < opts.devchange;
    convb = norm(beta-oldbeta)< opts.betathresh*sum(active);
    if convg && convb && oldfval<=fval
        flag='converged';
        break
    end

    if t==opts.trustiter
        flag = 'trust iterations exceeded';
    end
end

% and return
if iter==opts.coreiter && isempty(flag)
    flag = 'fisher iterations exceeded';
end
