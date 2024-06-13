function [theta0,theta_lb,theta_ub]=hypbounds(model)
%% Obtain hyperparameters bounds from model
%
% Inputs:
% model: struct with GPR model
%
% Outputs:
% theta0: initial guess
% theta_lb: lower bound
% theta_ub: upper bound
%

%%

theta0(model.idx.d,1)=model.ini.d;
theta0(model.idx.sigma_v,1)=model.ini.sigma_v;
theta0(model.idx.sigma,1)=model.ini.sigma;
theta0(model.idx.L,1)=model.ini.L;
theta0(model.idx.alpha,1)=model.ini.alpha0;
theta0(model.idx.abar,1)=model.ini.abar0;

theta_lb(model.idx.d,1)=model.lb.d;
theta_lb(model.idx.sigma_v,1)=model.lb.sigma_v;
theta_lb(model.idx.sigma,1)=model.lb.sigma;
theta_lb(model.idx.L,1)=model.lb.L;
theta_lb(model.idx.alpha,1)=model.lb.alpha;
theta_lb(model.idx.abar,1)=model.lb.abar;

theta_ub(model.idx.d,1)=model.ub.d;
theta_ub(model.idx.sigma_v,1)=model.ub.sigma_v;
theta_ub(model.idx.sigma,1)=model.ub.sigma;
theta_ub(model.idx.L,1)=model.ub.L;
theta_ub(model.idx.alpha,1)=model.ub.alpha;
theta_ub(model.idx.abar,1)=model.ub.abar;
