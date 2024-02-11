function [Ka_blk,Ka_grad_tmp]=kernel_cov(X1,X2,kernel_type,hyp_struct)

%% Covariance matrix from kernel function
%
% Inputs:
% X1: [n1,m] matrix with inputs as columns
% X2: [n2,m] matrix with inputs as columns
% model: struct with hyp
%
% Outputs:
% Ka_blk: output vector
% Ka_grad_tmp: cell with gradients of Ka wrt. hyperparameters
%

%%

if nargout>1
    return_grad=true;
else
    return_grad=false;
end

na=length(hyp_struct.sigma);

Ka_blk=cell(1,na);
Ka_grad_tmp=cell(2,na);

for k=1:na

    if strcmpi(kernel_type,'se')

        if return_grad
            [Ka_blk{k},Ka_grad_tmp(1:2,k)]=kernel_se(X1,X2,hyp_struct.sigma(k),hyp_struct.L(k));
        else
            [Ka_blk{k}]=kernel_se(X1,X2,hyp_struct.sigma(k),hyp_struct.L(k));
        end

    elseif strcmpi(kernel_type,'matern52')

        if return_grad
            [Ka_blk{k},Ka_grad_tmp(1:2,k)]=kernel_matern52(X1,X2,hyp_struct.sigma(k),hyp_struct.L(k));
        else
            [Ka_blk{k}]=kernel_matern52(X1,X2,hyp_struct.sigma(k),hyp_struct.L(k));
        end


    elseif strcmpi(kernel_type,'rq')

        Ka_grad_tmp{3,k}={};

        if return_grad
            [Ka_blk{k},Ka_grad_tmp(1:3,k)]=kernel_rq(X1,X2,hyp_struct.sigma(k),hyp_struct.L(k),hyp_struct.alpha(k));
        else
            [Ka_blk{k}]=kernel_rq(X1,X2,hyp_struct.sigma(k),hyp_struct.L(k),hyp_struct.alpha(k));
        end

    else
        error(['Kernel option not permissible: ' kernel_type]);
    end

end

