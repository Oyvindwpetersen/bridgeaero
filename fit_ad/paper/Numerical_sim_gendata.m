

clc
clear all
close all

addPathTemp2('C:\Cloud\OD_OWP\Work\Matlab\Common\shared\GPR\Test_GPR')

%%

at1=@(x) exp(-(x-3).*(x-3)/10)-2
at2=@(x) -0.5*exp(-(x-3).*(x-5))-2

a1=@(x) -0.5*exp(-(x-2).*(x-5)/3)-1
a2=@(x) 0.1*exp(-x.*(x-7)/5)-1 %-x/20
a2=@(x) 0.5*exp((x-3).*(x-8)/5.*(x-4))-1
a3=@(x) 2*exp(-(x-4).*(x-6)/3)-x/10
% a3=@(x) exp(-(x-1).*(x-6).*(x-5)/10)-1



d1=0.3
d2=0.5
d3=1.2

x_test=[3 3.6 4.3 5.1 6.2]

K_test1=linspace(0.05,1.2,8);
K_test2=linspace(0.1,1,6);
K_test3=linspace(0.2,1,6);
K_test4=linspace(0.1,1.3,8);
K_test5=linspace(0.2,1,8);

test_matrix=[
    gridvec(K_test1,x_test(1)) ; 
    gridvec(K_test2,x_test(2)) ; 
    gridvec(K_test3,x_test(3)) ; 
    gridvec(K_test4,x_test(4)) ; 
    gridvec(K_test5,x_test(5)) ];

%%

F_fun= @(K,x) at1(x)+at2(x).*1i.*K+...
    a1(x).*1i.*K./(1i*K+d1)+...
    a2(x).*1i.*K./(1i*K+d2)+...
    a3(x).*1i.*K./(1i*K+d3)

x_plot=[1:0.1:8].';
K_plot=[0:0.05:1.4].';

plot_matrix=gridvec(K_plot,x_plot);

F_plot=F_fun(plot_matrix(:,1),plot_matrix(:,2));
F_test=F_fun(test_matrix(:,1),test_matrix(:,2));

return
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
