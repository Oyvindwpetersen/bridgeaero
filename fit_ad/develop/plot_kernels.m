%%


L=1
alpha=1

x=linspace(0,5*L,100).';

[K1]=kernel_se(x,0,1,L);

[K2]=kernel_matern52(x,0,1,L);

[K3]=kernel_rq(x,0,1,L,alpha);

figure(); hold on; grid on;
plot(x,K1)
plot(x,K2);
plot(x,K3);
