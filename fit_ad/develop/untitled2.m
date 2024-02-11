%%


K_check_all=linspace(0,2,100)
for k=1:length(K_check_all)
K_check=0.01
K_check=K_check_all(k)
plot_matrix=[K_check 4]

[yp,yp_r,yp_i,std_yp_r,std_yp_i,std_yr_obs,std_yi_obs,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,plot_matrix,[yt_r;yt_i],model); 


[Dr,Di]=rf_matrix(model.d,K_check);


C=[Dr;Di]*cov_a_pred*[Dr;Di].'

c=cov2corr(C)


c21(k)=c(2,1)

end