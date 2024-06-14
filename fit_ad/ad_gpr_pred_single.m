function [yr_p,yi_p,std_yr_p,std_yi_p,std_yr_p_obs,std_yi_p_obs,ap,cov_ap]=ad_gpr_pred_single(test_matrix,pred_matrix,yr,yi,model)
%% Predict using GPR model
%
% Inputs:
% test_matrix: [M,2] matrix with K and x as columns
% pred_matrix: [N,2] matrix with K and x as columns
% yr_t: [M,1] vector with K^2*AD_stiffness
% yr_i: [M,1] vector with K^2*AD_damping
% model: struct with GPR model
%
% Outputs:
% yr_p: [N,1] vector with K^2*AD_stiffness
% yi_p: [N,1] vector with K^2*AD_damping
% std_yr_p: [N,1] vector with SD of yr_p
% std_yi_p: [N,1] vector with SD of yi_p
% std_yr_p_obs: [N,1] vector with SD of yr_p (+noise)
% std_yi_p_obs: [N,1] vector with SD of yr_p (+noise)
% ap: predicted a coefficients
%

%%

[Ka,Ky,xt_uni,St,D_glob]=ad_gpr(test_matrix,yr,yi,model);

y=[yr;yi];

% Difference between measurement and mean
if isempty(model.idx.abar)
    e=y;
else
    m=repmat(model.hyp.abar,length(xt_uni),1);
    e=y-D_glob*m;
end

beta=(D_glob).'/Ky*e;

[D_glob_pred,xp_uni]=rf_matrix_multi(pred_matrix,model.hyp.d);
[Sp]=restack_a(model.na,length(xp_uni));

[Ka_star_t_blk]=kernel_cov(xp_uni,xt_uni,model.kernel,model.hyp);
Ka_star_t=blkdiag2(Ka_star_t_blk{:});
Ka_star_t=Sp*Ka_star_t*St.';

if strcmpi(model.basis,'constant')
    a_mean_pred=repmat(model.hyp.abar,length(xp_uni),1);
    ap=a_mean_pred+Ka_star_t*beta;
else
    ap=Ka_star_t*beta;
end

yp=D_glob_pred*ap;

yr_p=yp(1:length(yp)/2);

yi_p=yp((length(yp)/2+1):end);

%% Uncertainty

[Ka_star_star_blk]=kernel_cov(xp_uni,xp_uni,model.kernel,model.hyp);
Ka_star_star=blkdiag2(Ka_star_star_blk{:});
Ka_star_star=Sp*Ka_star_star*Sp.';

% Uncertainty of a
cov_ap=Ka_star_star-Ka_star_t*(D_glob).'/Ky*(D_glob)*Ka_star_t.';

% Uncertainty of prediction
cov_yp=(D_glob_pred)*cov_ap*(D_glob_pred).';

n=size(cov_yp,1);

std_yr_p=diag(cov_yp(1:n/2,1:n/2)).^0.5;
std_yi_p=diag(cov_yp(n/2+1:end,n/2+1:end)).^0.5;

% Uncertainty of a noisy observation
N=noisecov(pred_matrix,model);

cov_y_obs=cov_yp+N;

std_yr_p_obs=diag(cov_y_obs(1:n/2,1:n/2)).^0.5;
std_yi_p_obs=diag(cov_y_obs(n/2+1:end,n/2+1:end)).^0.5;



