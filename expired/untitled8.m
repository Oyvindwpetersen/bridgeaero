%%

rng(0)
X1=sort(rand(20,1)*10);
X2=X1;
sigma=10
L=[4].'
alpha=[3].'

[K,K_grad]=kernel_rq(X1,X2,sigma,L,alpha);



dsigma=1e-6;
[Kp,~]=kernel_rq(X1,X2,sigma+dsigma,L,alpha);

K_grad_num=(Kp-K)/dsigma

ratio=K_grad_num./K_grad{1}-1


dL=1e-6;
[Kp,~]=kernel_rq(X1,X2,sigma,L+dL,alpha);

K_grad_num=(Kp-K)/dL

ratio=K_grad_num./K_grad{2}-1


dalpha=1e-3;
[Kp,~]=kernel_rq(X1,X2,sigma,L,alpha+dalpha);

K_grad_num=(Kp-K)/dalpha

ratio=K_grad_num./K_grad{3}-1