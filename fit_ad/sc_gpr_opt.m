function [model,neg_logL,neg_logL_grad,theta_opt]=sc_gpr_opt(test_matrix,y,model,varargin)

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

% Kernel is created for restacked a
%
% | y(alpha1,x1) | = | T(alpha1)   0       | | a(x1) |
% | y(alpha2,x2) |   | 0         T(alpha2) | | a(x2) |
%

%%

fun_obj= @(theta) sc_gpr_loglik_wrapper(test_matrix,y,model,theta);

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
    gs.NumStageOnePoints=20; % 200;
    gs.NumTrialPoints=100; % 1000;

    [theta_opt,logL_opt,~,~,manymins]=run(gs,problem);

end

%% Assign

model=hyp2model_poly(model,theta_opt);

[neg_logL,neg_logL_grad]=sc_gpr_loglik_wrapper(test_matrix,y,model,theta_opt);

