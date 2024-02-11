

clc
clear all

close all

addPathTemp2('C:\Cloud\OD_OWP\Work\Matlab\Common\shared\GPR\Test_GPR')

a1=@(x) cos(x);
a2=@(x) log(x+1);

a3=@(x) sqrt(1+x.^2)
a4=@(x) sin(x*3)
a5=@(x) -(1+x.^3).^0.25

d1=0.2
d2=0.4
d3=0.6

F_fun= @(K,x) a1(x)+a2(x).*1i.*K+a3(x).*1i.*K./(1i*K+d1) %+a4(x).*1i.*K./(1i*K+d2)+a5(x).*1i.*K./(1i*K+d2);

K_plot=linspace(0,1.2,100);
x_plot=linspace(3,7,100);

K_test=linspace(0.05,1,10);
x_test=linspace(4,6,5);


X_test=gridvec(K_test,x_test);
X_plot=gridvec(K_plot,x_plot);

F_plot=F_fun(X_plot(:,1),X_plot(:,2));
F_test=F_fun(X_test(:,1),X_test(:,2));

a1_plot=a1(X_plot(:,2));
a2_plot=a2(X_plot(:,2));
a3_plot=a3(X_plot(:,2));
a4_plot=a4(X_plot(:,2));
a5_plot=a5(X_plot(:,2));

figure(1); hold on; grid on;
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([X_plot],real(F_plot));
scatter3(X_test(:,1),X_test(:,2),real(F_test),'ok');


figure(2); hold on; grid on;
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([X_plot],imag(F_plot));
scatter3(X_test(:,1),X_test(:,2),imag(F_test),'ok');



%%

% 
% test_matrix=[K H]
% 
% close all
% 
% figure();
% scatter3(test_matrix(:,1),test_matrix(:,2),H1,'ob');
% 
% 
% figure();
% scatter3(test_matrix(:,1),test_matrix(:,2),H4,'ob');
% 
% 
% b_real=H1.*K.^2;
% b_imag=H4.*K.^2;
% 
% figure();
% scatter3(test_matrix(:,1),test_matrix(:,2),b_real,'ob');
% 
% 
% figure();
% scatter3(test_matrix(:,1),test_matrix(:,2),b_imag,'ob');
% 
% tilefigs

%%

%%%%%%

test_matrix=[X_test];
% y=[real(F_test) ; imag(F_test) ]

y=[real(F_test(1:10)) ; imag(F_test(1:10)) ;
    real(F_test(11:20)) ; imag(F_test(11:20)) ;
    real(F_test(21:30)) ; imag(F_test(21:30)) ;
    real(F_test(31:40)) ; imag(F_test(31:40)) ;
    real(F_test(41:50)) ; imag(F_test(41:50)) ]


b_real=real(F_test)
b_imag=imag(F_test)

La=ones(1,5);
sigma_a=ones(1,5)*10;
d=[d1 ]; %d2 d3
sigma_v=0.01

% La=ones(1,5);
% sigma_a=ones(1,5)*100;
% d=linspace(0.3,1,3)
% sigma_v=0.2
% y=[b_real;b_imag]

%%%%%%%%

na=2+length(d);

% test_matrix=
% [ K   x1  
%   .   .    
%   .   .    
%   .   .    
% ]

% Find unique x
[x_unique,~,ic]=unique(test_matrix(:,2));
nx=length(x_unique);
[Sa,Sa_inv]=restack_a(na,nx);

% Create bins of index for each unique x
bins=cell(1,nx);
for k=1:length(x_unique)
    bins{k}=find(ic==k);
end

% D-matrix for each unique x
D_blk=cell(1,nx);
for k=1:length(bins)
   K_this=test_matrix(bins{k},1);
   D_blk{k}=matrix(d,K_this);
end
D=blkdiag2(D_blk{:},'sparse');
Dr=D*Sa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_test=[a1(x_test(1)) a2(x_test(1)) a3(x_test(1)) ...
    a1(x_test(2)) a2(x_test(2)) a3(x_test(2)) ...
    a1(x_test(3)) a2(x_test(3)) a3(x_test(3)) ...
    a1(x_test(4)) a2(x_test(4)) a3(x_test(4)) ...
    a1(x_test(5)) a2(x_test(5)) a3(x_test(5)) ...
].'

F_test=D*a_test;
figure(1);
scatter3(X_test(:,1),X_test(:,2),F_test([[1:10] [21:30] [41:50] [61:70]  [81:90] ]),'xr');

figure(2);
scatter3(X_test(:,1),X_test(:,2),F_test([21:30] [41:50] [51:60] [71:80]  [91:100] ),'xr');

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% | b(x1,K1) | = | D(K1)   0   | | a(x1) |
% | b(x1,K1) |   | 0     D(K2) | | a(x2) |


% Noise matrix
N=sigma_v^2*speye(length(y));

% a=Sa*a_res;
% Kernel is created for restacked a

Ka_res_blk=cell(1,na);
for k=1:na
    Ka_res_blk{k}=kernel_se(x_unique,x_unique,sigma_a(k),La(k));
end
Ka_res=blkdiag2(Ka_res_blk{:},'sparse');

%%

% x_pred=linspace(0,12,100).';
x_pred=[5].';

Ka_res_blk2=cell(1,na);
for k=1:na
    Ka_res_blk2{k}=kernel_se(x_pred,x_unique,sigma_a(k),La(k));
end
Ka_res2=blkdiag2(Ka_res_blk2{:},'sparse');


Ka_res_star_blk2=cell(1,na);
for k=1:na
    Ka_res_star_blk2{k}=kernel_se(x_pred,x_pred,sigma_a(k),La(k));
end
Ka_res_star=blkdiag2(Ka_res_star_blk2{:},'sparse');

a_res_pred=Ka_res2*Dr.'/(Dr*Ka_res*Dr.'+N)*y

%%

alpha=Dr.'/(Dr*Ka_res*Dr.'+N)*y;


%%


b_real_plot=nan(length(X_plot),1);
b_imag_plot=nan(length(X_plot),1);

a1_pred=nan(length(X_plot),1);
a2_pred=nan(length(X_plot),1);
a3_pred=nan(length(X_plot),1);
a4_pred=nan(length(X_plot),1);
a5_pred=nan(length(X_plot),1);

for k1=1:length(X_plot)

    Ka_res_blk2=cell(1,na);
    for k=1:na
        Ka_res_blk2{k}=kernel_se(X_plot(k1,2),x_unique,sigma_a(k),La(k));
    end
    Ka_res2=blkdiag2(Ka_res_blk2{:},'sparse');
    
    a_res_pred=Ka_res2*alpha;
    
    Dp=matrix(d,X_plot(k1,1));

    y_pred=Dp*a_res_pred;


    b_real_plot(k1)=y_pred(1);
    b_imag_plot(k1)=y_pred(2);

    a1_pred(k1)=a_res_pred(1);
    a2_pred(k1)=a_res_pred(2);
    a3_pred(k1)=a_res_pred(3);
    % a4_pred(k1)=a_res_pred(4);
    % a5_pred(k1)=a_res_pred(5);
    % =[b_real_plot ; y_pred(1:101)];
    % b_imag_plot=[b_imag_plot ; y_pred(102:end)];


end
% a_pred=
close all


figure(); hold on; grid on;
scatter3(test_matrix(:,1),test_matrix(:,2),b_real,'ob');
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([X_plot],b_real_plot);


figure(); hold on; grid on;
scatter3(test_matrix(:,1),test_matrix(:,2),b_imag,'ob');
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([X_plot],b_imag_plot);


figure(); hold on; grid on; title('a1');
scatter3(test_matrix(:,1),test_matrix(:,2),a1(test_matrix(:,2)),'ob');
% scatter3(K_plot,x_plot,b_real_plot,'xr');
% surfiso([X_plot_full],a1_plot);
surfiso([X_plot],a1_pred);

tilefigs



return
%%

x_pred=[3.5:0.02:3.5].';
K_pred=[0:0.02:2].';


% x_pred=[4.5:0.01:6.5].';
% K_pred=[0:0.02:2].';

x_plot=[];
K_plot=[];
b_real_plot=[];
b_imag_plot=[];


for k1=1:length(x_pred)

    Ka_res_blk2=cell(1,na);
    for k=1:na
        Ka_res_blk2{k}=kernel_se(x_pred(k1),x_unique,sigma_a(k),La(k));
    end
    Ka_res2=blkdiag2(Ka_res_blk2{:},'sparse');
    
    a_res_pred=Ka_res2*Dr.'/(Dr*Ka_res*Dr.'+N)*y;
    
    x_plot=[x_plot ; x_pred(k1)*ones(length(K_pred),1)];
    K_plot=[K_plot ; K_pred];
    
    Dp=matrix(d,K_pred);

    y_pred=Dp*a_res_pred;


    b_real_plot=[b_real_plot ; y_pred(1:101)];
    b_imag_plot=[b_imag_plot ; y_pred(102:end)];


end
% a_pred=
close all


figure(); hold on; grid on;
scatter3(test_matrix(:,1),test_matrix(:,2),b_real,'ob');
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([K_plot x_plot],b_real_plot);


figure(); hold on; grid on;
scatter3(test_matrix(:,1),test_matrix(:,2),b_imag,'ob');
% scatter3(K_plot,x_plot,b_real_plot,'xr');
surfiso([K_plot x_plot],b_imag_plot);

% figure();
% scatter3(test_matrix(:,1),test_matrix(:,2),b_imag,'ob');

tilefigs


