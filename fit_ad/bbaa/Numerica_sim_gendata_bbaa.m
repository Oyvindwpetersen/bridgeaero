

clc
clear all
close all

addPathTemp2('C:\Cloud\OD_OWP\Work\Matlab\Common\shared\GPR\Test_GPR')

%%

at1=@(x) 2*cos(x/2);
at2=@(x) exp(-(x-3).*(x-5)/2)+1

a1=@(x) 1+2*sin(x/2) 
a2=@(x) -0.5*sqrt(9+x.^2)
a3=@(x) -3*exp(-(x-4).*(x-5))+x+1;

d1=0.2
d2=0.4
d3=0.6

K_test=linspace(0.05,1,8);
x_test=linspace(3,6,4+1);

%%

F_fun= @(K,x) at1(x)+at2(x).*1i.*K+...
    a1(x).*1i.*K./(1i*K+d1)+...
    a2(x).*1i.*K./(1i*K+d2)+...
    a3(x).*1i.*K./(1i*K+d3)

x_plot=[1:0.1:8].';
K_plot=[0:0.05:1.5].';

test_matrix=gridvec(K_test,x_test);

plot_matrix=gridvec(K_plot,x_plot);

F_plot=F_fun(plot_matrix(:,1),plot_matrix(:,2));
F_test=F_fun(test_matrix(:,1),test_matrix(:,2));

% return
figure(); hold on; grid on;
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([plot_matrix],real(F_plot));
scatter3(test_matrix(:,1),test_matrix(:,2),real(F_test),'ok');
zlabel('real')
view([135 10]);
% 
figure(); hold on; grid on;
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([plot_matrix],imag(F_plot));
scatter3(test_matrix(:,1),test_matrix(:,2),imag(F_test),'ok');
zlabel('imag')
view([135 10]);


% tilefigs([2 2],'l')

% return
%%

% close all 
figure(); hold on; grid on; sizefig('bl')
plot(x_plot,at1(x_plot),'--','DisplayName','at1');
plot(x_plot,at2(x_plot),'--','DisplayName','at2');
plot(x_plot,a1(x_plot),'DisplayName','a1');
plot(x_plot,a2(x_plot),'DisplayName','a2');
plot(x_plot,a3(x_plot),'DisplayName','a3');
% plot(x_plot,a4(x_plot),'DisplayName','a4');
% plot(x_plot,a5(x_plot),'DisplayName','a5');
legend show

tilefigs([2 2],'l')
