clc
clear all
% close all

importdata

%%

clc

nd=3;
na=nd+2;

d0=linspace(0.2,1.5,nd).';
d_lb=0.1*ones(nd,1);
d_ub=3.0*ones(nd,1);

model=modelopt('nd',nd,'kernel','se','noise','model3w','basis','zero','d0',d0,'d_lb',d_lb,'d_ub',d_ub,'sigma_v_ub',[10]);

% model.ini.d=[d1 d2 d3].'

% model.ini.sigma_v=0.1;
% model.ini.sigma=10*ones(na,1);
% model.ini.L=2*ones(na,1);
% 
% model.lb.d=0.1*ones(nd,1);
% model.lb.sigma_v=1e-3;
% model.lb.sigma=1e-2*ones(na,1);
% model.lb.L=1e-1*ones(na,1);
% 
% model.ub.d=3.0*ones(nd,1);
% model.ub.sigma_v=1e1;
% model.ub.sigma=1e2*ones(na,1);
% model.ub.L=1e2*ones(na,1);

% model.alpha0=1*ones(na,1);
% model.lb.alpha=1e-1*ones(na,1);
% model.ub.alpha=1e2*ones(na,1);

idx1=3;
idx2=3;

yt_r=[]; 
yt_i=[]; 
test_matrix=[]; 

yt_r(:,1)=data_K(idx1,idx2,:).^2.*data_AD_s(idx1,idx2,:);
yt_i(:,1)=data_K(idx1,idx2,:).^2.*data_AD_d(idx1,idx2,:);

test_matrix(:,1)=data_K(idx1,idx2,:);
test_matrix(:,2)=data_H(idx1,idx2,:);

%%

[logL_opt,d_opt,sigma_v_opt,sigma_opt,L_opt,alpha_opt]=ad_gpr_opt(test_matrix,yt_r,yt_i,model);

% return
% model.d=d_opt;
% model.sigma_v=sigma_v_opt;
% model.sigma=sigma_opt;
% model.L=L_opt;
% model.alpha=alpha_opt;

% rng(0)
% 
d_opt_all=[]; sigma_v_opt_all=[]; sigma_opt_all=[]; L_opt_all=[]; logL_opt_all=[];
for k=1:1

model.ini.d=rand(size(model.ini.d)).*(model.ub.d-model.lb.d)+model.lb.d;
model.ini.d=sort(model.ini.d);
% model.ini.d=[d1 d2 d3].'

model.ini.sigma=rand(size(model.ini.sigma)).*(10-0);
model.ini.L=rand(size(model.ini.L)).*(3-1)+1;


[logL_opt_all(:,k),d_opt_all(:,k),sigma_v_opt_all(:,k),sigma_opt_all(:,k),L_opt_all(:,k),~,~]=ad_gpr_opt(test_matrix,yt_r,yt_i,model);
% a_bar_opt_all(:,k)
%alpha_opt_all(:,k)
end

[logL_opt,idx_min]=min(logL_opt_all)

model.d=d_opt_all(:,idx_min);
model.sigma_v=sigma_v_opt_all(:,idx_min);
model.sigma=sigma_opt_all(:,idx_min);
model.L=L_opt_all(:,idx_min);
model.alpha=[];
model.abar=[];
% model.abar=a_bar_opt_all(:,idx_min);

return
%%

x_plot=[4.5:0.1:6.5].';
% x_plot=[5:0.1:7.5].';

K_plot=[0:0.05:1.5].';
plot_matrix=gridvec(K_plot,x_plot);

clc

[yp,yp_r,yp_i,std_yp_r,std_yp_i,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,plot_matrix,[yt_r;yt_i],model); 

% close all

plotopt=struct();
plotopt.xlabel='K';
plotopt.ylabel='x';
plotopt.view=[-105 30];
plotopt.view=[50 45];
plotopt.linestyle='-';
plotopt.linewidth=0.1;
plotopt.cbar=[0.5 1 0 0.1];
plotopt.cbar=NaN;
plotopt.xtick=[0 1];
plotopt.ytick=[2:2:8]

plotopt2=struct();
plotopt2.marker='o';
plotopt2.markersize=10;
plotopt2.color=[0 0 0];
plotopt2.linestyle='None';

figure(); 
ha=tight_subplot(2,2,[0.2 0.075],[0.15 0.15],[0.05 0.05]);

axes(ha(1)); hold on; grid on;
title('Prediction','FontWeight','normal');
surfiso(plot_matrix,yp_r,plotopt,'zlabel','K^2 H_1^*');
plot3(test_matrix(:,1),test_matrix(:,2),yt_r,plotopt2);

axes(ha(3)); hold on; grid on;
surfiso(plot_matrix,yp_i,plotopt,'zlabel','K^2 H_4^*');
plot3(test_matrix(:,1),test_matrix(:,2),yt_i,plotopt2);

axes(ha(2)); hold on; grid on;
title('Prediction uncertainty','FontWeight','normal');
surfiso(plot_matrix,std_yp_r,plotopt,'zlabel','K^2 H_1^*');
% zlim([0 0.05]);
% 
axes(ha(4)); hold on; grid on;
surfiso(plot_matrix,std_yp_i,plotopt,'zlabel','K^2 H_4^*');
% zlim([0 0.05]);

tilefigs
% plotscriptmain('h',5,'w',15,'name','Numerical_ex','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',8,'format',{'jpg'});

%%

% close all

K_pred=[0:0.05:1.3].';

plotopt=struct();
plotopt.title={ 'H=4.9 m' 'H=5.2 m' 'H=5.5 m' 'H=5.8 m' 'H=5.1 m' }
plotopt.gap=[0.1 0.05];

plotadsplit(test_matrix,yt_r,yt_i,model,K_pred,plotopt)


%%


[yp,yp_r,yp_i,std_yp_r,std_yp_i,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,test_matrix,[yt_r;yt_i],model); 


close all

figure(); hold on;
plot(yt_r,'ob');
plot(yp_r,'or');

figure(); hold on;
plot(yt_i,'ob');
plot(yp_i,'or');


figure(); hold on;
plot(yt_r-yp_r,'ob');
plot(yt_i-yp_i,'xb');

% figure(); hold on;
tilefigs






