function [neg_logL,neg_logL_grad]=poly_gpr_loglik(test_matrix,y,model,varargin)

%% Calculation of log likelihood and gradient
%
% Inputs:
% test_matrix: [M,2] matrix with alpha and x as columns
% y: [M,1] vector with outputs
% model: struct with GPR model
%
% Outputs:
% neg_logL: -log(L)
% neg_logL_grad: gradient of -log(L) wrt. hyperparameters
%

%% Parse inputs

p=inputParser;
addParameter(p,'gradient',true,@islogical)

parse(p,varargin{:})

gradient=p.Results.gradient;

%%

[Ka,Ky,xt_uni,St,T_glob,N,Ka_grad,N_grad]=poly_gpr(test_matrix,y,model);

%% Likelihood to be maximized

% Original expressions
% logL_bias=-0.5*y.'/Ky*y;
% logL_var=-0.5*log(det(Ky))

% Difference between measurement and mean
e=y;

% Model bias (fit penalty)
Ky_inv_e=Ky\e;
logL_bias=-0.5*e.'*Ky_inv_e;

% Model variance (complexity penalty)
[log_det_Ky_chol,L]=logdet(Ky,'chol');
logL_var=-0.5*log_det_Ky_chol;

logL=logL_bias+logL_var;
neg_logL=-logL;

if gradient==false
    neg_logL_grad=[];
    return
end

%% Gradient

% Number of hyperparameters
n_par=length(model.idx.sigma_v)+length(model.idx.sigma)+length(model.idx.L);

% Gradient of covariance matrix
Ky_grad=cell(n_par,1);

for j=1:n_par

    Ky_grad{j}=sparse(length(y),length(y));

    if any(j==model.idx.sigma_v)
        idx=j;
        Ky_grad{j}=N_grad{idx};
    elseif any(j==model.idx.sigma) | any(j==model.idx.L)
        idx=j-length(model.idx.sigma_v);
        Ky_grad{j}=T_glob*Ka_grad{idx}*T_glob.';
    end
    
end

logL_grad=nan*ones(n_par,1);

Ky_inv=Ky\eye(size(Ky));

for j=1:n_par

    % Original expressions
    % logL_grad(j,1)=0.5*y.'/Ky*Ky_grad{j}/Ky*y-0.5*trace(Ky\Ky_grad{j}); continue;
    % logL_grad(j,1)=0.5*(Ky_inv_y).'*Ky_grad{j}*Ky_inv_y-0.5*trace(Ky\Ky_grad{j});

    if any(j==model.idx.sigma_v) || any(j==model.idx.sigma) || any(j==model.idx.L)

        % logL_grad(j,1)=0.5*(Ky_inv_e).'*Ky_grad{j}*Ky_inv_e-0.5*trace(Ky_inv*Ky_grad{j});

        term1=0.5*(Ky_inv_e).'*Ky_grad{j}*Ky_inv_e;
        % term2=-0.5*trace(Ky_inv*Ky_grad{j});

        % trace(A*B)=sum(sum(A.'*B^T));
        term2=-0.5*sum(sum(Ky_inv.*Ky_grad{j}));
        logL_grad(j,1)=term1+term2;

    end

end

if any(isnan(logL_grad))
    logL_grad
    error('Gradient is NaN, aborting');
end

if any(isinf(logL_grad))
    logL_grad
    error('Gradient is Inf, aborting');
end

neg_logL_grad=-logL_grad;
