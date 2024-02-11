function [model,logL_opt,neg_logL_grad,neg_logL_grad_norm]=ad_gpr_opt(test_matrix,yt_r,yt_i,model,varargin)

%% Optimization
%
% Inputs:

%%

p=inputParser;
addParameter(p,'globalsearch',false,@islogical)

parse(p,varargin{:})

globalsearch=p.Results.globalsearch;

%%

y_std=std([yt_r;yt_i]);
L_char=max(test_matrix(:,2))-min(test_matrix(:,2));
y_mean=mean([yt_r;yt_i]);

if isempty(model.ini.sigma_v); model.ini.sigma_v=expandcol(y_std*1e-2,model.nv); end
if isempty(model.lb.sigma_v); model.lb.sigma_v=expandcol(y_std*1e-3,model.nv); end
if isempty(model.ub.sigma_v); model.ub.sigma_v=expandcol(y_std*1e0,model.nv); end

if isempty(model.ini.sigma); model.ini.sigma=expandcol(y_std*1e0,model.na); end
if isempty(model.lb.sigma); model.lb.sigma=expandcol(y_std*1e-1,model.na); end
if isempty(model.ub.sigma); model.ub.sigma=expandcol(y_std*1e1,model.na); end

if isempty(model.ini.L); model.ini.L=expandcol(L_char*1e0,model.na); end
if isempty(model.lb.L); model.lb.L=expandcol(L_char*1e-1,model.na); end
if isempty(model.ub.L); model.ub.L=expandcol(L_char*1e1,model.na); end

if ~isempty(model.idx.alpha)
if isempty(model.ini.alpha); model.ini.alpha=expandcol(1e0,model.na); end
if isempty(model.lb.alpha); model.lb.alpha=expandcol(1e-1,model.na); end
if isempty(model.ub.alpha); model.ub.alpha=expandcol(1e1,model.na); end
end

if ~isempty(model.idx.abar)
if isempty(model.ini.abar); model.ini.abar=expandcol(y_mean*1e0,model.na); end
if isempty(model.lb.abar); model.lb.abar=expandcol(y_mean*1e-1,model.na); end
if isempty(model.ub.abar); model.ub.abar=expandcol(y_mean*1e1,model.na); end
end

%% Assign initial and bounds

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

if any(theta_lb>=theta_ub)

    idx=find(theta_lb>theta_ub)
    error('theta lower bounds exceed upper bounds for elements index displayed above');

end

%%

% D_glob=diag(D_r1,D_r2,...,D_i1,D_i2,...)

% Kernel is created for restacked a

% | y_r(x1,K1) | = | D_r(K1)   0     | | a(x1) |
% | y_r(x2,K2) |   | 0       D_r(K2) | | a(x2) |
% | y_i(x1,K1) |   | D_i(K1)    0    |
% | y_i(x2,K2) |   |  0      D_i(K2) |
%

%%

y=[yt_r ; yt_i];

fun_obj= @(theta) ad_gpr_loglik_wrapper(test_matrix,y,model,theta);

options = optimoptions('fmincon','Display','iter','TypicalX',theta0,'ConstraintTolerance',1e-3,'OptimalityTolerance',1e-3,'StepTolerance',1e-3... ); %,...
    ,'CheckGradients',false,'Algorithm','interior-point','SpecifyObjectiveGradient',true);
% interior-point
% trust-region-reflective


% Ax<b
% d[i+1]-d[i]>delta_d -> -d[i+1]+d[i]<delta_d
% sigma_a> f*sigma_v -> -sigma_a+f*sigma_v<0 

% Constraint to separate poles d
A1=zeros(model.nd-1,length(theta0));
for k=1:(model.nd-1)
    A1(k,k-1+[1:2])=[1 -1];
end
b1=-model.delta_d*ones(model.nd-1,1);

f=1;

count=0;
A2=zeros(length(model.idx.sigma)*length(model.idx.sigma_v),length(theta0));
for k1=1:length(model.idx.sigma_v)
    for k2=1:length(model.idx.sigma)
        count=count+1;
        A2(count,[model.idx.sigma_v(k1) model.idx.sigma(k2) ])=[f -1];
    end
end
b2=zeros(size(A2,1),1);

A=[A1;A2];
b=[b1;b2];

if ~globalsearch

    [theta_opt,logL_opt]=fmincon(fun_obj,theta0,A,b,[],[],theta_lb,theta_ub,[],options);

end

%%

if globalsearch
    
    problem = createOptimProblem('fmincon','objective',fun_obj,'x0',theta0,...
        'lb',theta_lb,'ub',theta_ub,'Aineq',A,'bineq',b,'options',options);
    
    gs=GlobalSearch;
    gs.BasinRadiusFactor=0.5;
    gs.FunctionTolerance=1e-3;
    gs.XTolerance=1e-3;
    gs.StartPointsToRun='bounds-ineqs';
    gs.NumStageOnePoints=200; %200
    gs.NumTrialPoints=1000; %1000
    
    % t0=tic;
    [theta_opt_glob,logL_opt,~,~,manymins]=run(gs,problem);
    % t1=toc(t0);
    
    theta_opt=theta_opt_glob;
    
end

%%

d_opt=theta_opt(model.idx.d);
sigma_v_opt=theta_opt(model.idx.sigma_v);
sigma_opt=theta_opt(model.idx.sigma);
L_opt=theta_opt(model.idx.L);
alpha_opt=theta_opt(model.idx.alpha);
abar_opt=theta_opt(model.idx.abar);

model.hyp.d=d_opt;
model.hyp.sigma_v=sigma_v_opt;
model.hyp.sigma=sigma_opt;
model.hyp.L=L_opt;
model.hyp.alpha=alpha_opt;
model.hyp.abar=abar_opt;

model=orderstruct(model);

[neg_logL,neg_logL_grad]= ad_gpr_loglik_wrapper(test_matrix,y,model,theta_opt);
neg_logL_grad_norm=neg_logL_grad./theta_opt;