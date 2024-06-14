function [Ka,Ky,xt_uni,St,D_glob,N,Ka_grad,D_glob_grad,N_grad,abar_grad]=ad_gpr(test_matrix,yr,yi,model)
%% Basis for Gaussian process regression for aerodynamic derivatives
%
% Inputs:
% test_matrix: [M,2] matrix with K and x as columns
% yr: [M,1] vector with K^2*AD_stiffness
% yi: [M,1] vector with K^2*AD_damping
% model: struct with GPR model
%
% Outputs:
% Ka: covariance matrix for a
% Ky: covariance matrix for y
% xt_uni: unique values for x
% St: repopulation matrix for test data
% D_glob: regression matrix for test data
% N: covariance matrix for test data
% Ka_grad: gradient of Ka
% D_glob_grad: gradient of D_glob
% N_grad: gradient of N
%

%%

y=[yr;yi];

if size(test_matrix,1)~=length(y)/2
    error('Test matrix or data vector y has wrong size');
end

[D_glob,xt_uni,D_glob_grad]=rf_matrix_multi(test_matrix,model.hyp.d);

[St,St_inv]=restack_a(model.na,length(xt_uni));

%%

% Kernel for a_restack
% a_restack=[a1(x1) a1(x2) .... a1(xN) , a2(x1) a2(x2) ... a2(xN) , a3(x1) a3(x2) ... a3(xN)];

[Ka_blk,Ka_grad_tmp]=kernel_cov(xt_uni,xt_uni,model.kernel,model.hyp);
Ka=blkdiag2(Ka_blk{:},'sparse');
Ka=St*Ka*St.';

Ka_grad=cell(size(Ka_grad_tmp,1)*model.na,1);
idx=0;
for n=1:size(Ka_grad_tmp,1)
    for k=1:model.na

        idx=idx+1;
        range=(k-1)*length(xt_uni)+[1:length(xt_uni)];
        
        S_tmp=St(:,range);

        Ka_grad{idx}=S_tmp*sparse(Ka_grad_tmp{n,k})*S_tmp.';
    end
end

% Noise covariance
[N,N_grad]=noisecov(test_matrix,model);

% Covariance
Ky=(D_glob*Ka*D_glob.'+N);

%% Mean value

if strcmpi(model.basis,'constant')
    for k=1:model.na
        tmp=sparse(model.na,1); tmp(k)=1;
        abar_grad{k,1}=tmp;
    end
else
    abar_grad={};
end

