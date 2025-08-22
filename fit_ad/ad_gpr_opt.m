function [model,logL_opt,theta_opt]=ad_gpr_opt(test_matrix,yr_t,yr_i,model,varargin)

%% Optimization of GPR model for ADs
%
% Inputs:
% test_matrix: [M,2] matrix with K and x as columns
% yr_t: [M,1] vector with K^2*AD_stiffness
% yr_i: [M,1] vector with K^2*AD_damping
% model: struct with GPR model
%
% Outputs:
% model: struct with GPR model (optimized)
% logL_opt: struct with -log(L), gradient and hessian of -log(L) wrt. hyperparameters
% theta_opt: optimized hyperparameters
%

%%

p=inputParser;
addParameter(p,'globalsearch',false,@islogical)

parse(p,varargin{:})

globalsearch=p.Results.globalsearch;

%%  Ensure cell

cell_check=[iscell(test_matrix) iscell(yr_t) iscell(yr_i)];

if all(cell_check)
    % OK
    cell_out=true;
elseif ~any(cell_check)
    cell_out=false;
    test_matrix={test_matrix};
    yr_t={yr_t};
    yr_i={yr_i};
else
    error('All or none of inputs must be cell');
end

n1=size(test_matrix,1);
n2=size(test_matrix,2);

% If model is struct, then copy to all cells
if ~iscell(model)
    model_tmp=model; clear model
    for idx1=1:n1
        for idx2=1:n2
            model{idx1,idx2}=model_tmp;
        end
    end
end

%% Assign default initial values and bounds for hyperparameter

for idx1=1:n1
    for idx2=1:n2
        model{idx1,idx2}=hypdefault2model(test_matrix{idx1,idx2},yr_t{idx1,idx2},yr_i{idx1,idx2},model{idx1,idx2});
    end
end

%% Build initial and bounds

offset=model{1,1}.nd;

theta0=[]; theta_lb=[]; theta_ub=[]; 
for idx1=1:n1
    for idx2=1:n2

        [theta0_loc,theta_lb_loc,theta_ub_loc]=hypbounds(model{idx1,idx2});
        
        model{idx1,idx2}.idx.glob=[ [1:model{1,1}.nd] offset+[1:length(theta0_loc(model{1,1}.nd+1:end))] ];

        offset=offset+length(theta0_loc(model{1,1}.nd+1:end));

        theta0(model{idx1,idx2}.idx.glob,1)=theta0_loc;
        theta_lb(model{idx1,idx2}.idx.glob,1)=theta_lb_loc;
        theta_ub(model{idx1,idx2}.idx.glob,1)=theta_ub_loc;

        if any(theta_lb_loc>=theta_ub_loc)

            idx1
            idx2
            idx=find(theta_lb_loc>theta_ub_loc)
            error('theta lower bounds exceed upper bounds for elements index displayed above');

        end

    end
end

%%

% D_glob=diag(D_r1,D_r2,...,D_i1,D_i2,...)

% Kernel is created for restacked a

% | y_r(K1,x1) | = | D_r(K1)   0     | | a(x1) |
% | y_r(K2,x2) |   | 0       D_r(K2) | | a(x2) |
% | y_i(K1,x1) |   | D_i(K1)    0    |
% | y_i(K2,x2) |   |  0      D_i(K2) |
%

%%

fun_obj= @(theta) ad_gpr_loglik_wrapper(test_matrix,yr_t,yr_i,model,theta);

options = optimoptions('fmincon','Display','iter','TypicalX',theta0,'ConstraintTolerance',1e-3,'OptimalityTolerance',1e-3,'StepTolerance',1e-3,...
    'CheckGradients',false,'Algorithm','interior-point','SpecifyObjectiveGradient',true);

% Ax<b
% d[i+1]-d[i]>delta_d -> -d[i+1]+d[i]<delta_d
% sigma_a> f*sigma_v -> -sigma_a+f*sigma_v<0

% Constraint to separate poles d
A1=zeros(model{1,1}.nd-1,length(theta0));
for k=1:(model{1,1}.nd-1)
    A1(k,k-1+[1:2])=[1 -1];
end
b1=-model{1,1}.delta_d*ones(model{1,1}.nd-1,1);

for idx1=1:n1
    for idx2=1:n2

        f=1;
        count=0;
        A2=zeros(length(model{idx1,idx2}.idx.sigma)*length(model{idx1,idx2}.idx.sigma_v),length(theta0));
        for k1=1:length(model{idx1,idx2}.idx.sigma_v)
            for k2=1:length(model{idx1,idx2}.idx.sigma)
                count=count+1;
                A2(count,model{idx1,idx2}.idx.glob([model{idx1,idx2}.idx.sigma_v(k1) model{idx1,idx2}.idx.sigma(k2) ]))=[f -1];
            end
        end
        b2=zeros(size(A2,1),1);

        A2_all{idx1,idx2}=A2;
        b2_all{idx1,idx2}=b2;

    end
end

A=[A1;cell2mat(A2_all(:))];
b=[b1;cell2mat(b2_all(:))];

% options2=options; options2.MaxIterations=10;
% rng(0);
% for k=1:20
%     d0=sort(rand(3,1).*(model{1,1}.ub.d-model{1,1}.lb.d)+model{1,1}.lb.d)
%     d0_all(k,:)=d0;
% 
%     theta0_2=theta0; theta0_2(1:3)=d0;
%     [theta_opt,logL_opt(k)]=fmincon(fun_obj,theta0_2,A,b,[],[],theta_lb,theta_ub,[],options2);
% end

%% Perform optimization

if ~globalsearch
    [theta_opt,logL_opt]=fmincon(fun_obj,theta0,A,b,[],[],theta_lb,theta_ub,[],options);
end

if globalsearch

    problem = createOptimProblem('fmincon','objective',fun_obj,'x0',theta0,...
        'lb',theta_lb,'ub',theta_ub,'Aineq',A,'bineq',b,'options',options);

    gs=GlobalSearch;
    gs.BasinRadiusFactor=0.5;
    gs.FunctionTolerance=1e-3;
    gs.XTolerance=1e-3;
    gs.StartPointsToRun='bounds-ineqs';
    gs.NumStageOnePoints=100; % 200;
    gs.NumTrialPoints=200; % 1000;

    [theta_opt,logL_opt,~,~,manymins]=run(gs,problem);

end

%% Assign

for idx1=1:n1
    for idx2=1:n2

        theta_loc=theta_opt(model{idx1,idx2}.idx.glob);
        model{idx1,idx2}=hyp2model(model{idx1,idx2},theta_loc);

        model{idx1,idx2}=orderstruct(model{idx1,idx2});
    end
end

[neg_logL,neg_logL_grad,neg_logL_hess,neg_logL_grad_modelfit,neg_logL_grad_complexity]=ad_gpr_loglik_wrapper(test_matrix,yr_t,yr_i,model,theta_opt);

% Assign log(L) to struct
logL_opt=struct();
logL_opt.info='Optimized negative log(L)';
logL_opt.val=neg_logL;
logL_opt.grad=neg_logL_grad;
logL_opt.hess=neg_logL_hess;
logL_opt.gradm=neg_logL_grad_modelfit;
logL_opt.gradc=neg_logL_grad_complexity;

n_par=length(model{1,1}.idx.glob);

N_points=0;
for idx1=1:n1
    for idx2=1:n2
        N_points=N_points+size(test_matrix{idx1,idx2},1);
    end
end

logL_opt.AIC=-2*log(-neg_logL)+2*n_par;
logL_opt.BIC=-2*log(-neg_logL)+n_par*log(N_points);

%%

if cell_out==false
    model=model{1,1};
end
