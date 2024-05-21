function [neg_logL,neg_logL_grad]=ad_gpr_loglik(test_matrix,y,model,varargin)

%% Optimization
%
% Inputs:

%%

p=inputParser;
addParameter(p,'gradient',true,@islogical)

parse(p,varargin{:})

gradient=p.Results.gradient;

%%

[Ka,Ky,beta,xt_uni,Sa,D_glob,N,Ka_grad,D_glob_grad,N_grad,abar_grad]=ad_gpr(test_matrix,y,model);

%% Likelihood to be maximized

% logL_bias=-0.5*y.'/Ky*y;
% logL_var=-0.5*log(det(Ky))

if isempty(model.idx.abar)
    e=y;
else
    m=repmat(model.hyp.abar,length(xt_uni),1);
    e=y-D_glob*m;
end

Ky_inv_e=Ky\e;
logL_bias=-0.5*e.'*Ky_inv_e;

[log_det_Ky_chol,L]=logdet(Ky,'chol');
logL_var=-0.5*log_det_Ky_chol;

logL=logL_bias+logL_var;

if gradient==false
    neg_logL=-logL;
    neg_logL_grad=[];
    return
end

%%

n_par=length(model.idx.d)+length(model.idx.sigma_v)+length(model.idx.sigma)+length(model.idx.L)+length(model.idx.alpha)+length(model.idx.abar);

Ky_grad=cell(n_par,1);
m_grad=cell(n_par,1);

for j=1:n_par

    Ky_grad{j}=sparse(length(y),length(y));
    m_grad{j}=sparse(model.na*length(xt_uni),1);

    if any(j==model.idx.d)
        % Ky_grad{j}=D_glob_grad{j}*Ka*D_glob.'+D_glob*Ka*D_glob_grad{j}.';

        tmp=full(D_glob_grad{j}*Ka*D_glob.');
        Ky_grad{j}=tmp+tmp.';

        if strcmpi(model.noise,'model_a')
            error('Not implemented yet');
            % Ky_grad{j}=Ky_grad{j}+D_glob_grad{j}*N*D_glob.'+D_glob*model.hyp.sigma_v(3)^2*D_glob_grad{j}.';
        end

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

for j=1:n_par

    % logL_grad(j,1)=0.5*y.'/Ky*Ky_grad{j}/Ky*y-0.5*trace(Ky\Ky_grad{j});
    % logL_grad(j,1)=0.5*(Ky_inv_y).'*Ky_grad{j}*Ky_inv_y-0.5*trace(Ky\Ky_grad{j});

    if any(j==model.idx.d) | any(j==model.idx.sigma_v) | any(j==model.idx.sigma) | any(j==model.idx.L) | any(j==model.idx.alpha)
        logL_grad(j,1)=0.5*(Ky_inv_e).'*Ky_grad{j}*Ky_inv_e-0.5*trace(L.' \ (L \Ky_grad{j}));
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

neg_logL=-logL;
neg_logL_grad=-logL_grad;
