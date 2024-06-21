function [Ka,Ky,xt_uni,St,T_glob,N,Ka_grad,N_grad]=poly_gpr(test_matrix,y,model)
%% Basis for Gaussian process regression for poly
%
% Inputs:
% test_matrix: [M,2] matrix with alpha and x as columns
% y: [M,1] vector with outputs
% model: struct with GPR model
%
% Outputs:
% Ka: covariance matrix for a
% Ky: covariance matrix for y
% xt_uni: unique values for x
% St: repopulation matrix for test data
% T_glob: polynomial matrix for test data
% N: covariance matrix for test data
% Ka_grad: gradient of Ka
% N_grad: gradient of N
%

%%

if size(test_matrix,1)~=length(y)
    error('Test matrix or data vector y has wrong size');
end

[T_glob,xt_uni]=poly_matrix_multi(test_matrix,model.p);

[St,St_inv]=restack_a(model.p+1,length(xt_uni));

%%

% Kernel for a_restack
% a_restack=[a1(x1) a1(x2) .... a1(xN) , a2(x1) a2(x2) ... a2(xN) , a3(x1) a3(x2) ... a3(xN)];

[Ka_blk,Ka_grad_tmp]=kernel_cov(xt_uni,xt_uni,model.kernel,model.hyp);
Ka=blkdiag2(Ka_blk{:},'sparse');
Ka=St*Ka*St.';

Ka_grad=cell(size(Ka_grad_tmp,1)*(model.p+1),1);
idx=0;
for n=1:size(Ka_grad_tmp,1)
    for k=1:(model.p+1)

        idx=idx+1;
        range=(k-1)*length(xt_uni)+[1:length(xt_uni)];
        
        S_tmp=St(:,range);

        Ka_grad{idx}=S_tmp*sparse(Ka_grad_tmp{n,k})*S_tmp.';
    end
end

% Noise covariance
N=model.hyp.sigma_v.^2*eye(size(test_matrix,1));
N_grad{1}=2*model.hyp.sigma_v*eye(size(test_matrix,1));

% Covariance
Ky=(T_glob*Ka*T_glob.'+N);


