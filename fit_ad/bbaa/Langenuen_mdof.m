%%

c=0;
for idx1=2:3
for idx2=2:3
% idx1=2;
% idx2=2;

yt_r=[]; 
yt_i=[]; 
test_matrix=[]; 

yt_r(:,1)=data_K(idx1,idx2,:).^2.*data_AD_s(idx1,idx2,:);
yt_i(:,1)=data_K(idx1,idx2,:).^2.*data_AD_d(idx1,idx2,:);

test_matrix(:,1)=data_K(idx1,idx2,:);
test_matrix(:,2)=data_H(idx1,idx2,:);

%%

[d_opt,sigma_v_opt,sigma_opt,L_opt,alpha_opt]=ad_gpr_opt(test_matrix,yt_r,yt_i,model);

rng(0)

d_opt_all=[]; sigma_v_opt_all=[]; sigma_opt_all=[]; L_opt_all=[]; logL_opt_all=[];
for k=1:2

model.ini.d=rand(size(model.ini.d)).*(model.ub.d-model.lb.d)+model.lb.d;
model.ini.d=sort(model.ini.d);
% model.ini.d=[d1 d2 d3].'

model.ini.sigma=rand(size(model.ini.sigma)).*(10-0);
model.ini.L=rand(size(model.ini.L)).*(3-1)+1;


[logL_opt_all(:,k),d_opt_all(:,k),sigma_v_opt_all(:,k),sigma_opt_all(:,k),L_opt_all(:,k),~,~]=ad_gpr_opt(test_matrix,yt_r,yt_i,model);

%alpha_opt_all(:,k)
end

[logL_opt,idx_min]=min(logL_opt_all)

model.d=d_opt_all(:,idx_min);
model.sigma_v=sigma_v_opt_all(:,idx_min);
model.sigma=sigma_opt_all(:,idx_min);
model.L=L_opt_all(:,idx_min);
model.alpha=[];
model.abar=[];

%%

x_plot=[4.5:0.1:6.5].';
% x_plot=[5:0.1:7.5].';

K_plot=[0:0.05:1.5].';
plot_matrix=gridvec(K_plot,x_plot);


clc

[yp,yp_r,yp_i,std_yp_r,std_yp_i,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,plot_matrix,[yt_r;yt_i],model); 

% yp_r_plot(idx1,idx2,:)=yp_r;
% yp_i_plot(idx1,idx2,:)=yp_i;

c=c+1;
yp_r_plot{idx1-1,idx2-1}=yp_r;
yp_i_plot{idx1-1,idx2-1}=yp_i;

plot_matrix_all{idx1-1,idx2-1}=plot_matrix;
test_matrix_all{idx1-1,idx2-1}=test_matrix;

fr_point_all{idx1-1,idx2-1}=yt_r;
fi_point_all{idx1-1,idx2-1}=yt_i;

end
end

return
%%

close all

plotopt=struct();

plotopt.gap=[0.15 0.175];
plotopt.marg_h=[0.15 0.05];
plotopt.marg_w=[0.125 0.05];

plotopt.view=[120 60];

plotopt.xtick=[0 1];
plotopt.ytick=[4:0.5:7];
plotopt.xlabel='K';
plotopt.ylabel='D [m]';

plotopt.zlabel={'K^2 H_4^*' 'K^2 H_3^*' 'K^2 A_4^*' 'K^2 A_3^*'}

surfisomulti(plot_matrix_all,yp_r_plot,test_matrix_all,fr_point_all,plotopt)

% plotscriptmain('h',5,'w',8,'name','Langenuen_AD_stiff','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',8,'format',{'jpg'});



plotopt.zlabel={'K^2 H_1^*' 'K^2 H_2^*' 'K^2 A_1^*' 'K^2 A_2^*'}

surfisomulti(plot_matrix_all,yp_i_plot,test_matrix_all,fi_point_all,plotopt)

% plotscriptmain('h',5,'w',8,'name','Langenuen_AD_damp','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',8,'format',{'jpg'});

tilefigs
