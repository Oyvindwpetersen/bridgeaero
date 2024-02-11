function [y_pred,yr_pred,yi_pred,std_yr_pred,std_yi_pred,std_yr_obs,std_yi_obs,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,pred_matrix,y,model)

%%


%%

[Ka,Ky,beta,xt_uni,St,D_glob]=ad_gpr(test_matrix,y,model);

[D_glob_pred,xp_uni]=rf_matrix_multi(pred_matrix,model.hyp.d);
[Sp]=restack_a(model.na,length(xp_uni));

[Ka_star_t_blk]=kernel_cov(xp_uni,xt_uni,model.kernel,model.hyp);
Ka_star_t=blkdiag2(Ka_star_t_blk{:});
Ka_star_t=Sp*Ka_star_t*St.';

if strcmpi(model.basis,'constant')
    a_mean_pred=repmat(model.abar,length(xp_uni),1);
    a_pred=a_mean_pred+Ka_star_t*beta;
else
    a_pred=Ka_star_t*beta;
end

y_pred=D_glob_pred*a_pred;

yr_pred=y_pred(1:length(y_pred)/2);

yi_pred=y_pred((length(y_pred)/2+1):end);

%% Uncertainty

[Ka_star_star_blk]=kernel_cov(xp_uni,xp_uni,model.kernel,model.hyp);
Ka_star_star=blkdiag2(Ka_star_star_blk{:});
Ka_star_star=Sp*Ka_star_star*Sp.';

% Uncertainty of a
cov_a_pred=Ka_star_star-Ka_star_t*(D_glob).'/Ky*(D_glob)*Ka_star_t.';

% Uncertainty of prediction
cov_y_pred=(D_glob_pred)*cov_a_pred*(D_glob_pred).';

n=size(cov_y_pred,1);

std_yr_pred=diag(cov_y_pred(1:n/2,1:n/2)).^0.5;
std_yi_pred=diag(cov_y_pred(n/2+1:end,n/2+1:end)).^0.5;

% Uncertainty of a noisy observation

N=noisecov(pred_matrix,model);

cov_y_obs=cov_y_pred+N;

std_yr_obs=diag(cov_y_obs(1:n/2,1:n/2)).^0.5;
std_yi_obs=diag(cov_y_obs(n/2+1:end,n/2+1:end)).^0.5;



