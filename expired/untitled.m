%%



[Ka,Ky,alpha,xt_uni,Sa,I2,D_glob,D_glob_grad,Ka_grad]=ad_gpr(test_matrix,y,d,sigma_v,sigma_a,L_a);


d2=d;
d2(1)=d2(1)+1e-6;


[Ka2,Ky2,alpha2,xt_uni2,Sa2,I22,D_glob2,D_glob_grad2,Ka_grad2]=ad_gpr(test_matrix,y,d2,sigma_v,sigma_a,L_a);









dD_num=full((D_glob2-D_glob))

dD_dd_num=full((D_glob2-D_glob))/1e-6
