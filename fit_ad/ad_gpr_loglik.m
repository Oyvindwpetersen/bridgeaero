function [neg_logL,neg_logL_grad,neg_logL_grad_modelfit,neg_logL_grad_complexity]=ad_gpr_loglik(test_matrix,yr,yi,model,varargin)

%% Calculation of log likelihood and gradient
%
% Inputs:
% test_matrix: [M,2] matrix with K and x as columns
% yr: [M,1] vector with K^2*AD_stiffness
% yi: [M,1] vector with K^2*AD_damping
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

[Ka,Ky,xt_uni,Sa,D_glob,N,Ka_grad,D_glob_grad,N_grad,abar_grad]=ad_gpr(test_matrix,yr,yi,model);

y=[yr;yi];

%% Likelihood to be maximized

% Original expressions
% logL_bias=-0.5*y.'/Ky*y;
% logL_var=-0.5*log(det(Ky))

% Difference between measurement and mean
if isempty(model.idx.abar)
    e=y;
else
    m=repmat(model.hyp.abar,length(xt_uni),1);
    e=y-D_glob*m;
end

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
n_par=length(model.idx.d)+length(model.idx.sigma_v)+length(model.idx.sigma)+length(model.idx.L)+length(model.idx.alpha)+length(model.idx.abar);

% Gradient of covariance matrix and mean
Ky_grad=cell(n_par,1);
m_grad=cell(n_par,1);

for j=1:n_par

    Ky_grad{j}=sparse(length(y),length(y));
    m_grad{j}=sparse(model.na*length(xt_uni),1);

    if any(j==model.idx.d)
        % Ky_grad{j}=D_glob_grad{j}*Ka*D_glob.'+D_glob*Ka*D_glob_grad{j}.';

        tmp=full(D_glob_grad{j}*Ka*D_glob.');
        Ky_grad{j}=tmp+tmp.';
        
    elseif any(j==model.idx.sigma_v)
        idx=j-length(model.idx.d);
        Ky_grad{j}=N_grad{idx};
    elseif any(j==model.idx.sigma) | any(j==model.idx.L) | any(j==model.idx.alpha)
        idx=j-length(model.idx.d)-length(model.idx.sigma_v);
        Ky_grad{j}=D_glob*Ka_grad{idx}*D_glob.';
    elseif any(j==model.idx.abar)
        idx=j-length(model.idx.d)-length(model.idx.sigma_v)-length(model.idx.sigma)-length(model.idx.L)-length(model.idx.alpha);
        m_grad{j}=repmat(abar_grad{idx},length(xt_uni),1);
    end
    
end

logL_grad=nan*ones(n_par,1);
logL_grad_modelfit=nan*ones(n_par,1);
logL_grad_complexity=nan*ones(n_par,1);


Ky_inv=Ky\eye(size(Ky));

for j=1:n_par

    % Original expressions
    % logL_grad(j,1)=0.5*y.'/Ky*Ky_grad{j}/Ky*y-0.5*trace(Ky\Ky_grad{j});
    % logL_grad(j,1)=0.5*(Ky_inv_y).'*Ky_grad{j}*Ky_inv_y-0.5*trace(Ky\Ky_grad{j});

    if any(j==model.idx.d) || any(j==model.idx.sigma_v) || any(j==model.idx.sigma) || any(j==model.idx.L) || any(j==model.idx.alpha)

        % logL_grad(j,1)=0.5*(Ky_inv_e).'*Ky_grad{j}*Ky_inv_e-0.5*trace(Ky_inv*Ky_grad{j});

        logL_grad_modelfit(j,1)=0.5*(Ky_inv_e).'*Ky_grad{j}*Ky_inv_e;
        % term2=-0.5*trace(Ky_inv*Ky_grad{j});

        % trace(A*B)=sum(sum(A.'*B^T));
        logL_grad_complexity(j,1)=-0.5*sum(sum(Ky_inv.*Ky_grad{j}));
        logL_grad(j,1)=logL_grad_modelfit(j,1)+logL_grad_complexity(j,1);

    elseif any(j==model.idx.abar)
        logL_grad(j,1)=-0.5*(-D_glob*m_grad{j}).'*Ky_inv_e*2;
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
neg_logL_grad_modelfit=-logL_grad_modelfit;
neg_logL_grad_complexity=-logL_grad_complexity;



