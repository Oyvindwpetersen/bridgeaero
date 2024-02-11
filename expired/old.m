
% 
% BB_fun=@(K,x) tanh(-K/10+).*sin(2*pi*exp(-x/10))+tanh(-K/10+).*sin(2*pi*exp(-x/10));
% 
% 
% K_vec=linspace(0,10,100);
% x_vec=linspace(0,20,100);
% 
% BB_vec=
% 
% surfiso([K_vec.' x_vec.'],f)


%%
% 
% clc
% clear all
% close all
% 
% %%
% K1=[1 2 3];
% K2=[1.2 2.1 3.1];
% K3=[1.2 2.1 3.1 4];
% 
% x1=[5 5 5]
% x2=[7 7 7]
% x3=[10 10 10 10]
% 
% test_matrix=[ K1' x1' ; K2' x2' ; K3' x3' ]
% 
% test_matrix=[ K1' ones(size(x1))' ; K2' ones(size(x2))' ; K3' ones(size(x3))' ]
% 
% bb=tanh(test_matrix(:,1).^2/20+1i*test_matrix(:,2));
% 
% b_real=real(bb);
% b_imag=imag(bb);



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