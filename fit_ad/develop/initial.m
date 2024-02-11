%%

n_check=1000

d0_all=rand(model.nd,n_check).*(model.ub.d-model.lb.d)+model.lb.d;
d0_all=sort(d0_all,1);

sigma_v0_all=10.^( rand(model.nv,n_check).*(-3+1)-1)*std(y)

sigma0_all=rand(model.na,n_check)*((10-0.1)+0.1)*std(y);

L_span=max(test_matrix(:,2))-min(test_matrix(:,2));
L0_all=rand(model.na,n_check)*((2-0.1)+0.1)*L_span;




d0_all(:,1)=model.ini.d;
sigma_v0_all(:,1)=model.ini.sigma_v;
sigma0_all(:,1)=model.ini.sigma;
L0_all(:,1)=model.ini.L;


for k=1:n_check

    hyper_check=hyper;
    hyper_check.d=d0_all(:,k);
    hyper_check.sigma_v=sigma_v0_all(:,k);
    hyper_check.sigma=sigma0_all(:,k);
    hyper_check.L=L0_all(:,k);

    hyper_check.alpha=[];
    hyper_check.abar=[];
    

    neglogL_all(k)=ad_gpr_loglik(test_matrix,y,hyper_check,'gradient',false);

end



fun_obj= @(theta) ad_gpr_loglik_wrapper(test_matrix,y,hyper


