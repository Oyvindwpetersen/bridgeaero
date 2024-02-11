function [Ka,Ky,beta,xt_uni,St,D_glob,N,Ka_grad,D_glob_grad,N_grad,abar_grad]=ad_gpr(test_matrix,y,model)

%%
%
% Inputs:
% test_matrix: column vector
% y: column vector
% model: struct with settings
%
% Outputs:
% Ka: covariance matrix for a
% Ky: covariance matrix for y
% beta:
% xt_uni: unique values for
% Sa:
% D_glob:
% N:
% Ka_grad:
% D_glob_grad:
% N_grad:


%%

if size(test_matrix,1)~=length(y)/2
    error('Test matrix or data vector y has wrong size');
end

[D_glob,xt_uni,D_glob_grad]=rf_matrix_multi(test_matrix,model.hyp.d);

[St,St_inv]=restack_a(model.na,length(xt_uni));

%%

[N,N_grad]=noisecov(test_matrix,model);

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

if strcmpi(model.basis,'constant')
    for k=1:model.na
        tmp=sparse(model.na,1); tmp(k)=1;
        abar_grad{k,1}=tmp;
    end
else
    abar_grad={};
end

%%

Ky=(D_glob*Ka*D_glob.'+N);

if strcmpi(model.basis,'constant')
    tmp=repmat(model.abar,length(xt_uni),1);
    e=y-D_glob*tmp;
else
    e=y;
end

beta=(D_glob).'/Ky*e;

%%

