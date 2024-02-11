

clc
clear all
close all

Generate_data2

%%

clc

nd=3;
na=nd+2;

model=modelopt('nd',nd,'kernel','se','noise','model2','basis','zero');

model.ini.d=linspace(0.2,1.5,nd).';
% model.ini.d=[d1 d2 d3].'

model.ini.sigma_v=0.1;
model.ini.sigma=10*ones(na,1);
model.ini.L=2*ones(na,1);

model.lb.d=0.1*ones(nd,1);
model.lb.sigma_v=1e-3;
model.lb.sigma=1e-1*ones(na,1);
model.lb.L=1e-1*ones(na,1);

model.ub.d=2.0*ones(nd,1);
model.ub.sigma_v=1e1;
model.ub.sigma=1e2*ones(na,1);
model.ub.L=1e1*ones(na,1);

% model.alpha0=1*ones(na,1);
% model.lb.alpha=1e-1*ones(na,1);
% model.ub.alpha=1e2*ones(na,1);

rng(0);
yt_r=real(F_test)+randn(size(F_test))*0.05;
yt_i=imag(F_test)+randn(size(F_test))*0.05;

%%
clc

[logL_opt,d_opt,sigma_v_opt,sigma_opt,L_opt,alpha_opt,abar_opt]=ad_gpr_opt(test_matrix,yt_r,yt_i,model);

% return
% model.d=d_opt;
% model.sigma_v=sigma_v_opt;
% model.sigma=sigma_opt;
% model.L=L_opt;
% model.alpha=alpha_opt;

rng(0)

for k=1:2

model.ini.d=rand(size(model.ini.d)).*(model.ub.d-model.lb.d)+model.lb.d;
model.ini.d=sort(model.ini.d);
% model.ini.d=[d1 d2 d3].'

model.ini.sigma=rand(size(model.ini.sigma)).*(10-0);
model.ini.L=rand(size(model.ini.L)).*(3-1)+1;

[logL_opt_all(:,k),d_opt_all(:,k),sigma_v_opt_all(:,k),sigma_opt_all(:,k),L_opt_all(:,k),~,~]=ad_gpr_opt(test_matrix,yt_r,yt_i,model);
%alpha_opt_all(:,k)
% abar_opt_all(:,k)
end

[~,idx_min]=min(logL_opt_all)

model.d=d_opt_all(:,idx_min);
model.sigma_v=sigma_v_opt_all(:,idx_min);
model.sigma=sigma_opt_all(:,idx_min);
model.L=L_opt_all(:,idx_min);
model.alpha=[];
model.abar=[];
% model.abar=abar_opt_all(:,idx_min);

%%

clc

[yp,yp_r,yp_i,std_yp_r,std_yp_i,std_yr_obs,std_yi_obs,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,plot_matrix,[yt_r;yt_i],model); 

close all

plotopt=struct();
plotopt.xlabel='K';
plotopt.ylabel='x';
plotopt.view=[-105 30];
plotopt.view=[-240 40];

plotopt.linestyle='-';
plotopt.linewidth=0.1;
plotopt.cbar=[0.5 1 0 0.1];
plotopt.cbar=NaN;
plotopt.xtick=[0 1];
plotopt.ytick=[2:2:8]

plotopt2=struct();
plotopt2.marker='o';
plotopt2.markersize=1;
plotopt2.color=[0 0 0];
plotopt2.linestyle='None';

figure(); 
ha=tight_subplot(2,4,[0.2 0.075],[0.15 0.15],[0.05 0.05]);

axes(ha(1)); hold on; grid on;
title('Truth','FontWeight','normal');
surfiso(plot_matrix,real(F_plot),plotopt,'zlabel','K^2 H_4^*');
plot3(test_matrix(:,1),test_matrix(:,2),yt_r,plotopt2);
% zlim([-2 4]);

axes(ha(1+4)); hold on; grid on;
surfiso(plot_matrix,imag(F_plot),plotopt,'zlabel','K^2 H_1^*');
plot3(test_matrix(:,1),test_matrix(:,2),yt_i,plotopt2);
% zlim([0 5]);

axes(ha(2)); hold on; grid on;
title('Prediction','FontWeight','normal');
surfiso(plot_matrix,yp_r,plotopt,'zlabel','K^2 H_4^*');
% zlim([-2 4]);

axes(ha(2+4)); hold on; grid on;
surfiso(plot_matrix,yp_i,plotopt,'zlabel','K^2 H_1^*');
% zlim([0 5]);

axes(ha(3)); hold on; grid on;
title('Error (abs.)','FontWeight','normal');
surfiso(plot_matrix,abs(real(F_plot)-yp_r),plotopt,'zlabel','K^2 H_4^*');
% zlim([0 2]);

axes(ha(3+4)); hold on; grid on;
surfiso(plot_matrix,abs(imag(F_plot)-yp_i),plotopt,'zlabel','K^2 H_1^*');
% zlim([0 2]);

axes(ha(4)); hold on; grid on;
title('Prediction uncertainty','FontWeight','normal');
surfiso(plot_matrix,std_yp_r,plotopt,'zlabel','K^2 H_4^*');
% zlim([0 2]);

axes(ha(4+4)); hold on; grid on;
surfiso(plot_matrix,std_yp_i,plotopt,'zlabel','K^2 H_1^*');
% zlim([0 2]);

tilefigs
% plotscriptmain('h',5,'w',15,'name','Numerical_ex2','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',8,'format',{'jpg'});

return
%%

% close all

K_pred=[0:0.05:1.5].';

plotopt=struct();
plotopt.title={'H=3.00 m' 'H=3.75 m' 'H=4.50 m' 'H=5.25 m' 'H=6.00 m' }
plotopt.gap=[0.1 0.05];

plotadsplit(test_matrix,yt_r,yt_i,model,K_pred,plotopt)

tilefigs

% x_pred=4 %5.25-2;
% K_pred=[0:0.05:1.5].';
% 
% pred_matrix=gridvec(K_pred,x_pred);
% 
% [y_pred,yr_pred,yi_pred,std_yp_r,std_yp_i]=ad_gpr_pred(test_matrix,pred_matrix,[yt_r;yt_i],model); 
% 
% F_plot2=F_fun(pred_matrix(:,1),pred_matrix(:,2));
% 
% % close all
% 
% figure(); hold on; grid on;
% plot(K_pred,real(F_plot2))
% plot(K_pred,yr_pred)
% h_shade=plotci(K_pred,yr_pred,std_yp_r,2);
% 
% figure(); hold on; grid on;
% plot(K_pred,imag(F_plot2))
% plot(K_pred,yi_pred)
% h_shade=plotci(K_pred,yi_pred,std_yp_i,2);
% 
% tilefigs

%%



