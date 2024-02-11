
clc
clear all
close all

Numerical_sim_gendata

%%

clc

nd=3;
na=nd+2;

model=modelopt('nd',nd,'kernel','se','noise','model2','basis','zero');

model.ub.db=2;
model.lb.d=0.1;


% model.ini.d=linspace(0.2,1.5,nd).';
% % model.ini.d=[d1 d2 d3].'
% 
% model.ini.sigma_v=0.1;
% model.ini.sigma=10*ones(na,1);
% model.ini.L=2*ones(na,1);
% 
% model.lb.d=0.1*ones(nd,1);
% model.lb.sigma_v=1e-3;
% model.lb.sigma=1e-1*ones(na,1);
% model.lb.L=1e-1*ones(na,1);
% 
% model.ub.d=3.0*ones(nd,1);
% model.ub.sigma_v=1e1;
% model.ub.sigma=1e2*ones(na,1);
% model.ub.L=1e1*ones(na,1);

model=modelopt_fix(model);

rng(0);
yt_r=real(F_test)+randn(size(F_test))*0.1;
yt_i=imag(F_test)+randn(size(F_test))*0.2;

%% Optimize

clc

% [logL_opt,model]
[model,logL_opt,neg_logL_grad,neg_logL_grad_norm]=ad_gpr_opt(test_matrix,yt_r,yt_i,model,'globalsearch',true);

%% Predict
clc

[yp,yp_r,yp_i,std_yp_r,std_yp_i,std_yr_obs,std_yi_obs,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,plot_matrix,[yt_r;yt_i],model); 

close all

plotopt=struct();
plotopt.xlabel='K';
plotopt.ylabel='x';
% plotopt.view=[-105 30];
plotopt.view=[-240 40];
plotopt.view=[-260 40];

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
ha=tight_subplot(2,4,[0.2 0.1],[0.15 0.15],[0.05 0.05]);

axes(ha(1)); hold on; grid on;
title('Truth','FontWeight','normal');
surfiso(plot_matrix,real(F_plot),plotopt,'zlabel','K^2 H_4^*');
plot3(test_matrix(:,1),test_matrix(:,2),yt_r,plotopt2);
zlim([-4 0]);

axes(ha(1+4)); hold on; grid on;
surfiso(plot_matrix,imag(F_plot),plotopt,'zlabel','K^2 H_1^*');
plot3(test_matrix(:,1),test_matrix(:,2),yt_i,plotopt2);
zlim([-4 0]);

axes(ha(2)); hold on; grid on;
title('Prediction','FontWeight','normal');
surfiso(plot_matrix,yp_r,plotopt,'zlabel','K^2 H_4^*');
zlim([-4 0]);

axes(ha(2+4)); hold on; grid on;
surfiso(plot_matrix,yp_i,plotopt,'zlabel','K^2 H_1^*');
zlim([-4 0]);

axes(ha(3)); hold on; grid on;
title('Error (abs.)','FontWeight','normal');
surfiso(plot_matrix,abs(real(F_plot)-yp_r),plotopt,'zlabel','K^2 H_4^*');
zlim([0 2]);

axes(ha(3+4)); hold on; grid on;
surfiso(plot_matrix,abs(imag(F_plot)-yp_i),plotopt,'zlabel','K^2 H_1^*');
zlim([0 2]);

axes(ha(4)); hold on; grid on;
title('Prediction uncertainty','FontWeight','normal');
surfiso(plot_matrix,2*std_yp_r,plotopt,'zlabel','K^2 H_4^*');
zlim([0 2]);

axes(ha(4+4)); hold on; grid on;
surfiso(plot_matrix,2*std_yp_i,plotopt,'zlabel','K^2 H_1^*');
zlim([0 2]);

tilefigs

tight_subplot_letter(ha,6,[-0.025 -0.075],'()');

% plotscriptmain('h',5,'w',15,'name','Num_ex_surf','path','fig','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',8,'format',{'pdf' 'jpg'});
% plotscriptmain('h',5,'w',15,'name','Num_ex_surf','path','fig','renderer','opengl','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',8,'format',{'pdf' 'jpg'});
return

%%

close all

K_pred=[0:0.05:1.5].';

plotopt=struct();
for k=1:length(x_test)
    plotopt.title{k}=['x=' num2str(x_test(k),'%0.2f') ' m'];
end
plotopt.ylabel={'H_4^* K^2' 'H_1^* K^2'};
plotopt.markersize=3;
plotopt.gap=[0.15 0.05];
plotopt.marg_h=[0.15 0.2];
plotopt.marg_w=[0.05 0.05];

ha=plotadsplit(test_matrix,yt_r,yt_i,model,K_pred,plotopt);

plotadsplitref(ha,K_pred,real(F_fun(K_pred,x_test)),imag(F_fun(K_pred,x_test)),'displayname','Truth');

tilefigs

legendadjust('east',[0.06 0.075],3,20);

% plotscriptmain('h',5,'w',15,'name','Num_ex_split','path','fig','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',6,'format',{'pdf' 'jpg'});

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


