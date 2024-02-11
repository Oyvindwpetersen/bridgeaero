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

% na=model.na;

% d=model.hyp.d;
% sigma_v=model.hyp.sigma_v;
% sigma=model.hyp.sigma;
% L=model.hyp.L;
% alpha_a=model.hyp.alpha;
% abar=model.abar;

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


% OLD
% [Ka_blk,Ka_grad_tmp]=kernel_cov_old(xt_uni,xt_uni,model);
% Ka=blkdiag2(Ka_blk{:},'sparse');
% Ka=St*Ka*St.';
% 
% idx=0;
% for n=1:length(Ka_grad_tmp{1})
% 
%     offset_col=0;
%     offset_row=0;
% 
%     for k=1:model.na
% 
%         idx=idx+1;
%         range_col=offset_col+[1:size(Ka_grad_tmp{k}{n},1)];
%         range_row=offset_row+[1:size(Ka_grad_tmp{k}{n},1)];
% 
%         tmp=sparse(size(Ka,1),size(Ka,2));
%         tmp(range_col,range_row)=Ka_grad_tmp{k}{n};
% 
%         Ka_grad_old{idx}=St*tmp*St.';
% 
%         offset_col=range_col(end);
%         offset_row=range_row(end);
% 
%     end
% end

% NEW
[Ka_blk,Ka_grad_tmp]=kernel_cov(xt_uni,xt_uni,model.kernel,model.hyp);
Ka=blkdiag2(Ka_blk{:},'sparse');
Ka=St*Ka*St.';

Ka_grad=cell(size(Ka_grad_tmp,1)*model.na,1);
idx=0;
for n=1:size(Ka_grad_tmp,1)
    for k=1:model.na

        idx=idx+1;
        range=(k-1)*length(xt_uni)+[1:length(xt_uni)];
        % 
        % tmp=sparse(size(Ka,1),size(Ka,2));
        % tmp(range,range)=Ka_grad_tmp{n,k};
        % Ka_grad{idx}=St*tmp*St.';
        % 
        % 
        % S1=[sparse((k-1)*length(xt_uni),length(xt_uni)) ; speye(length(xt_uni)) ; sparse((model.na-k)*length(xt_uni),length(xt_uni)) ];
        % S=St*S1;
        % Ka_grad{idx}=St*S1*sparse(Ka_grad_tmp{n,k})*S2*St.';
        % Ka_grad{idx}=S*sparse(Ka_grad_tmp{n,k})*S.';

        S_test=St(:,range);

        Ka_grad{idx}=S_test*sparse(Ka_grad_tmp{n,k})*S_test.';
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

