%%

% idx_d=[1:nd];
% idx_sigma_v=idx_d(end)+1;
% idx_sigma=idx_sigma_v(end)+[1:na];
% idx_L=idx_sigma(end)+[1:na];
% 
% theta0(idx_d,1)=model.ini.d;
% theta0(idx_sigma_v,1)=model.ini.sigma_v;
% theta0(idx_sigma,1)=model.ini.sigma;
% theta0(idx_L,1)=model.ini.L;

theta=[2.6 ...
    1.5 ...
    3 5 10 ...
    2 0.5 6]';

theta_lb=-inf*theta;
theta_ub=inf*theta;

[d,sigma_v,sigma_a,L_a]=theta_split(theta,hypermeta.idx_d,hypermeta.idx_sigma_v,hypermeta.idx_sigma,hypermeta.idx_L);


[Ka,Ky,alpha,xt_uni,Sa,I2,D_glob,D_glob_grad,Ka_grad]=ad_gpr(test_matrix,y,d,sigma_v,sigma_a,L_a);

%%

dp=1e-8;
[~,Ky_perturb,~,~,~,~,D_glob_perturb,~,~]=ad_gpr(test_matrix,y,d+dp,sigma_v,sigma_a,L_a);

D_glob_grad_num=(D_glob_perturb-D_glob)/dp

ratio=D_glob_grad_num./D_glob_grad{1}-1


%%

D_glob_t=(D_glob*I2*Sa)

D_glob_grad_t{1}=(D_glob_grad{1}*I2*Sa)
j=1;

Ky_grad=D_glob_grad_t{j}*Ka*D_glob_t.'+D_glob_t*Ka*D_glob_grad_t{j}.';


Ky_grad_num=(Ky_perturb-Ky)/dp

ratio=Ky_grad_num./Ky_grad-1


% [Ka,Ky,alpha,xt_uni,Sa,I2,D_glob,D_glob_grad,Ka_grad]=ad_gpr(test_matrix,y,d,sigma_v,sigma_a,L_a);
