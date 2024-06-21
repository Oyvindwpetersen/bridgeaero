function [model,neg_logL,neg_logL_grad,theta_opt]=poly_gpr_opt(test_matrix,y,model,varargin)

%% Optimization of GPR model for poly model
%
% Inputs:
% test_matrix: [M,2] matrix with alpha and x as columns
% y: [M,1] vector with outputs
% model: struct with GPR model
%
% Outputs:
% model: struct with GPR model (optimized)
% neg_logL: -log(L)
% neg_logL_grad: gradient of -log(L) wrt. hyperparameters
% theta_opt: optimized hyperparameters
%

%%

p=inputParser;
addParameter(p,'globalsearch',false,@islogical)

parse(p,varargin{:})

globalsearch=p.Results.globalsearch;


%% Assign default initial values and bounds for hyperparameter

model=hypdefault2model_poly(test_matrix,y,model);

%% Build initial and bounds

theta0(model.idx.sigma_v,1)=model.ini.sigma_v;
theta0(model.idx.sigma,1)=model.ini.sigma;
theta0(model.idx.L,1)=model.ini.L;

theta_lb(model.idx.sigma_v,1)=model.ini.sigma_v;
theta_lb(model.idx.sigma,1)=model.lb.sigma;
theta_lb(model.idx.L,1)=model.lb.L;

theta_ub(model.idx.sigma_v,1)=model.ub.sigma_v;
theta_ub(model.idx.sigma,1)=model.ub.sigma;
theta_ub(model.idx.L,1)=model.ub.L;

%%

% D_glob=diag(D_r1,D_r2,...,D_i1,D_i2,...)

% Kernel is created for restacked a

% | y_r(x1,K1) | = | D_r(K1)   0     | | a(x1) |
% | y_r(x2,K2) |   | 0       D_r(K2) | | a(x2) |
% | y_i(x1,K1) |   | D_i(K1)    0    |
% | y_i(x2,K2) |   |  0      D_i(K2) |
%

%%

% for idx1=1:n1
%     for idx2=1:n2
%         y{idx1,idx2}=[yr_t{idx1,idx2} ; yr_i{idx1,idx2}];
%     end
% end

% [neg_logL_test,neg_logL_grad_test]=ad_gpr_loglik_wrapper(test_matrix,yr_t,yr_i,model,theta0);

fun_obj= @(theta) poly_gpr_loglik_wrapper(test_matrix,y,model,theta);

options = optimoptions('fmincon','Display','iter','TypicalX',theta0,'ConstraintTolerance',1e-3,'OptimalityTolerance',1e-3,'StepTolerance',1e-3,...
    'CheckGradients',false,'Algorithm','interior-point','SpecifyObjectiveGradient',true);

%% Perform optimization

if ~globalsearch
    [theta_opt,logL_opt]=fmincon(fun_obj,theta0,[],[],[],[],theta_lb,theta_ub,[],options);
end

if globalsearch

    problem = createOptimProblem('fmincon','objective',fun_obj,'x0',theta0,...
        'lb',theta_lb,'ub',theta_ub,'options',options);

    gs=GlobalSearch;
    gs.BasinRadiusFactor=0.5;
    gs.FunctionTolerance=1e-3;
    gs.XTolerance=1e-3;
    gs.StartPointsToRun='bounds-ineqs';
    gs.NumStageOnePoints=200; % 200;
    gs.NumTrialPoints=1000; % 1000;

    [theta_opt,logL_opt,~,~,manymins]=run(gs,problem);

end

%% Assign

model=hyp2model_poly(model,theta_opt);

[neg_logL,neg_logL_grad]=poly_gpr_loglik_wrapper(test_matrix,y,model,theta_opt);

