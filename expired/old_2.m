


%%



K_vec=linspace(0,2,30);
X_vec=linspace(3,6.5,40);


d=theta_opt(idx_d)
sigma_v=theta_opt(idx_sigma_v)
sigma_a=theta_opt(idx_sigma)
La=theta_opt(idx_L)

plot_matrix=gridvec(K_vec,X_vec)


%%


function pred_ad(alpha)

na=length(sigma_a);

% Find unique x
[xp_uni,~,ic]=unique(pred_matrix(:,2));
nx=length(xp_uni);
[Sa,Sa_inv]=restack_a(na,nx);

% Create bins of index for each unique x
bins=cell(1,nx);
for k=1:length(xp_uni)
    bins{k}=find(ic==k);
end

[D_glob]=rf_matrix_multi(pred_matrix,d,bins);

I_mat=[speye(na*nx);speye(na*nx)];

D_glob_t=D_glob*I_mat*Sa;


Ka_res_blk2=cell(1,na);
for k=1:na
    Ka_res_blk2{k}=kernel_se(x_pred(k1),xp_uni,sigma_a(k),La(k));
end
Ka_res2=blkdiag2(Ka_res_blk2{:},'sparse');

a_res_pred=Ka_res2*alpha;

x_plot=[x_plot ; x_pred(k1)*ones(length(K_pred),1)];
K_plot=[K_plot ; K_pred];

Dp=matrix(d,K_pred);

y_pred=Dp*a_res_pred;

%%

% Kernel for a_restack
% a_restack=[a1(x1) a1(x2) .... a1(xN) , a2(x1) a2(x2) ... a2(xN) , a3(x1) a3(x2) ... a3(xN)];


% x_pred=[4.5:0.01:6.5].';
% K_pred=[0:0.02:2].';

x_plot=[];
K_plot=[];
b_real_plot=[];
b_imag_plot=[];


for k1=1:length(x_pred)

    Ka_res_blk2=cell(1,na);
    for k=1:na
        Ka_res_blk2{k}=kernel_se(x_pred(k1),xp_uni,sigma_a(k),La(k));
    end
    Ka_res2=blkdiag2(Ka_res_blk2{:},'sparse');

    a_res_pred=Ka_res2*Dr.'/(Dr*Ka_res*Dr.'+N)*y;

    x_plot=[x_plot ; x_pred(k1)*ones(length(K_pred),1)];
    K_plot=[K_plot ; K_pred];

    Dp=matrix(d,K_pred);

    y_pred=Dp*a_res_pred;


    b_real_plot=[b_real_plot ; y_pred(1:101)];
    b_imag_plot=[b_imag_plot ; y_pred(102:end)];


end
% a_pred=
close all


figure(); hold on; grid on;
scatter3(test_matrix(:,1),test_matrix(:,2),b_real,'ob');
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([K_plot x_plot],b_real_plot);


figure(); hold on; grid on;
scatter3(test_matrix(:,1),test_matrix(:,2),b_imag,'ob');
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([K_plot x_plot],b_imag_plot);

% figure();
% scatter3(test_matrix(:,1),test_matrix(:,2),b_imag,'ob');

tilefigs