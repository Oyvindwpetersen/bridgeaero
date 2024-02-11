

%%


clc
clear all
close all



fun_obj= @(x) (x(:,1)-5).^2+x(:,2).*log(x(:,1)+15)+x(:,2).^2;


fun_grad= @(x) [ ...
                2*x(1)-10+x(2)./(x(1)+15) ;...
                log(x(1)+15)+2*x(2)];

fun_grad2= @(x) [2-x(2)./(x(1)+15).^2 ; 1./(x(1)+15)
                    1./(x(1)+15)   ; 2 ];


x_plot=linspace(-10,10,100)
y_plot=linspace(-10,10,100)


xy_plot=gridvec(x_plot,y_plot);

f_plot=fun_obj(xy_plot)


figure();hold on; grid on;
surfiso([xy_plot],(f_plot));

tilefigs([2 2]);




theta0=[-5 -5]
theta_lb=[-100 -100];
theta_ub=[100 100];

options = optimoptions('fmincon','Display','iter','TypicalX',theta0);
[theta_opt,f_opt]=fmincon(@testfun,theta0,[],[],[],[],theta_lb,theta_ub,[],options);




%%

clc

options = optimoptions('fmincon','Display','iter','TypicalX',theta0,'CheckGradients',true,...
    'Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'HessianFcn','objective');

[theta_opt,f_opt]=fmincon(@testfun,theta0,[],[],[],[],theta_lb,theta_ub,[],options);
 


