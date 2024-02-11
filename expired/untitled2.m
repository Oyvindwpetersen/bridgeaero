%%

clc

nd=1
na=nd+2;

idx_d=[1:nd];
idx_sigma_v=idx_d(end)+1;
idx_sigma=idx_sigma_v(end)+[1:na];
idx_L=idx_sigma(end)+[1:na];

theta0=[]
theta0(idx_d,1)=model.ini.d;
theta0(idx_sigma_v,1)=model.ini.sigma_v;
theta0(idx_sigma,1)=model.ini.sigma;
theta0(idx_L,1)=model.ini.L;

theta_lb(idx_d,1)=model.d_lb;
theta_lb(idx_sigma_v,1)=model.sigma_v_lb;
theta_lb(idx_sigma,1)=model.sigma_lb;
theta_lb(idx_L,1)=model.L_lb;

theta_ub(idx_d,1)=model.d_ub;
theta_ub(idx_sigma_v,1)=model.sigma_v_ub;
theta_ub(idx_sigma,1)=model.sigma_ub;
theta_ub(idx_L,1)=model.L_ub;

hypermeta.idx_d=idx_d;
hypermeta.idx_sigma_v=idx_sigma_v;
hypermeta.idx_sigma=idx_sigma;
hypermeta.idx_L=idx_L;

y=[Fr ; Fi];
theta=theta0;

[logL,logL_grad]=ad_gpr_loglik(test_matrix,y,hypermeta,theta)

%%


logLp=[];
for k=1:length(theta)
    thetap=theta;
    thetap(k)=thetap(k)+1e-3;
    
    [logLp(k,1)]=ad_gpr_loglik(test_matrix,y,hypermeta,thetap);

end



logL_grad_num=(logLp-logL)./1e-3

delta=logL_grad_num-logL_grad


ratio=delta./logL_grad