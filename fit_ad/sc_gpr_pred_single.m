function [yp,std_y_p,std_y_p_obs,ap,cov_ap]=sc_gpr_pred_single(test_matrix,pred_matrix,y,model)
%% Predict using GPR model
%
% Inputs:
% test_matrix: [M,2] matrix with alpha and x as columns
% pred_matrix: [M,2] matrix with alpha and x as columns
% y: [M,1] vector with outputs
% model: struct with GPR model
%
% Outputs:
% yp: [N,1] vector with K^2*AD_stiffness
% std_y_p: [N,1] vector with SD of yr_p
% std_y_p_obs: [N,1] vector with SD of yr_p (+noise)
% ap: predicted a coefficients
% cov_ap: covariance of a coefficients
%

%%

[Ka,Ky,xt_uni,St,T_glob]=sc_gpr(test_matrix,y,model);

% Difference between measurement and mean
e=y;

beta=(T_glob).'/Ky*e;

[T_glob_pred,xp_uni]=poly_matrix_multi(pred_matrix,model.p);

[Sp]=restack_a(model.p+1,length(xp_uni));

Ka_star_t_blk=kernel_cov(xp_uni,xt_uni,model.kernel,model.hyp);
Ka_star_t=blkdiag2(Ka_star_t_blk{:});
Ka_star_t=Sp*Ka_star_t*St.';

ap=Ka_star_t*beta;

yp=T_glob_pred*ap;

%% Uncertainty

Ka_star_star_blk=kernel_cov(xp_uni,xp_uni,model.kernel,model.hyp);
Ka_star_star=blkdiag2(Ka_star_star_blk{:});
Ka_star_star=Sp*Ka_star_star*Sp.';

% Uncertainty of a
cov_ap=Ka_star_star-Ka_star_t*(T_glob).'/Ky*(T_glob)*Ka_star_t.';

% Uncertainty of prediction
cov_yp=(T_glob_pred)*cov_ap*(T_glob_pred).';
std_y_p=diag(cov_yp).^0.5;

% Uncertainty of a noisy observation
N=model.hyp.sigma_v.^2*eye(size(pred_matrix,1));
cov_y_obs=cov_yp+N;
std_y_p_obs=diag(cov_y_obs).^0.5;
