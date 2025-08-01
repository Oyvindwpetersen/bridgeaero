function [N,N_grad]=noisecov(data_matrix,model)

%% Noise matrix
%
% Inputs:
% data_matrix: [M,2] matrix with K and x as columns
% model: struct with GPR model
%
% Outputs:
% N: [2M,2M] noise matrix
% N_grad: gradient of N wrt. sigma_v
%
% -------------------------------------
%
% Model 1: homogeneous noise on TF
%
% |y| = |K^2*ADs|+|vs|
%       |K^2*ADd|+|vd|
%
% Cov(V)=sigma^2*diag(I,I)
%
% -------------------------------------
%
% Model 2: separate noise on TF
%
% |y| = |K^2*ADs|+|vs|
%       |K^2*ADd|+|vd|
%
% Cov(V)=diag(sigma1^2*I,sigma2^2*I)
%
% -------------------------------------
%
% Model 3: separate noise, on stiffness and damping
%
% |y'| = |K^2*ADs|+|vs|
%        |K*ADd|  +|vd|
%
% |y|  = |K^2*ADs|+|vs|
%        |K^2*ADd|+|K*vd|
%
% Cov(V)=diag(sigma1^2*I,K^2*sigma2^2*I)
%
% -------------------------------------
%
% Model 3w: same as 3, but added white noise on last term
%
% |y'| = |K^2*ADs|+|vs|
%        |K*ADd|  +|vd|
%
% |y|  = |K^2*ADs|+|vs|
%        |K^2*ADd|+|K*vd|+|wd|
%
% Cov(V)=diag(sigma1^2*I,K^2*sigma2^2*I)+diag(0,sigma3^2*I)
%
% -------------------------------------
%
% Model 4: separate noise, on ADs
%
% |y'| = |ADs|  +|vs|
%        |ADd|  +|vd|
%
% |y|  = |K^2*ADs|+|K^2*vs|
%        |K^2*ADd|+|K^2*vd|
%
% Cov(V)=diag(K^4*sigma1^2*I,K^4*sigma2^2*I)
%
% -------------------------------------
%
% Model 4w: same as 4, but added white noise on both terms
%
% |y|  = |K^2*ADs|+|K^2*vs|+|ws|
%        |K^2*ADd|+|K^2*vd|+|wd|
%
% Cov(V)=diag(K^4*sigma1^2*I,K^4*sigma2^2*I)+diag(sigma3^2*I,sigma4^2*I)
%
% -------------------------------------
%
% Model 5: separate noise, on ADs*K
%
% |y'| = |ADs*K|  +|vs|
%        |ADd*K|  +|vd|
%
% |y|  = |K^2*ADs|+|K*vs|
%        |K^2*ADd|+|K*vd|
%
% Cov(V)=diag(K^2*sigma1^2*I,K^2*sigma2^2*I)
%
% -------------------------------------
%
% Model 5w: same as 5, but added white noise on both terms
%
% |y|  = |K^2*ADs|+|K^2*vs|+|ws|
%        |K^2*ADd|+|K^2*vd|+|wd|
%
% Cov(V)=diag(K^2*sigma1^2*I,K^2*sigma2^2*I)+diag(sigma3^2*I,sigma4^2*I)
%
% -------------------------------------
%

%%

n=size(data_matrix,1);
In=speye(n);
Zn=In*0;

K=data_matrix(:,1);

if strcmpi(model.noise,'model3') | strcmpi(model.noise,'model3w') | ...
        strcmpi(model.noise,'model4') | strcmpi(model.noise,'model4w') | ...
        strcmpi(model.noise,'model5') | strcmpi(model.noise,'model5w') 

    % K_diag=diag(K);
    K_diag=spdiags(K(:),0,numel(K),numel(K));
    K_squared_diag=K_diag.^2;
    K_cubic_diag=K_diag.^3;
    K_quad_diag=K_diag.^4;
end

%% Noise matrix

if strcmpi(model.noise,'model1')

    N=model.hyp.sigma_v^2*eye(2*n);

    N_grad{1}=2*model.hyp.sigma_v(1)*eye(2*n);

elseif strcmpi(model.noise,'model2')

    N=blkdiag2(model.hyp.sigma_v(1)^2*In,model.hyp.sigma_v(2)^2*In);

    N_grad{1}=blkdiag2(model.hyp.sigma_v(1)^2*In,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,model.hyp.sigma_v(2)^2*In,'sparse');

elseif strcmpi(model.noise,'model3')

    N=blkdiag2(model.hyp.sigma_v(1)^2*In,model.hyp.sigma_v(2)^2*K_squared_diag);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*In,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_squared_diag,'sparse');

elseif strcmpi(model.noise,'model3w')

    N=blkdiag2(model.hyp.sigma_v(1)^2*In,model.hyp.sigma_v(2)^2*K_squared_diag)...
      +blkdiag2(Zn,model.hyp.sigma_v(3)^2*In);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*In,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_squared_diag,'sparse');
    N_grad{3}=blkdiag2(Zn,2*model.hyp.sigma_v(3)*In,'sparse');

elseif strcmpi(model.noise,'model4')

    N=blkdiag2(model.hyp.sigma_v(1)^2*K_quad_diag,model.hyp.sigma_v(2)^2*K_quad_diag);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*K_quad_diag,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_quad_diag,'sparse');
    
elseif strcmpi(model.noise,'model4w')

    N=blkdiag2(model.hyp.sigma_v(1)^2*K_quad_diag,model.hyp.sigma_v(2)^2*K_quad_diag)...
     +blkdiag2(model.hyp.sigma_v(3)^2*In,model.hyp.sigma_v(4)^2*In);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*K_quad_diag,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_quad_diag,'sparse');
    N_grad{3}=blkdiag2(model.hyp.sigma_v(3)^2*In,Zn,'sparse');
    N_grad{4}=blkdiag2(Zn,model.hyp.sigma_v(4)^2*In,'sparse');
    
elseif strcmpi(model.noise,'model5')

    N=blkdiag2(model.hyp.sigma_v(1)^2*K_squared_diag,model.hyp.sigma_v(2)^2*K_squared_diag);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*K_squared_diag,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_squared_diag,'sparse');
    
elseif strcmpi(model.noise,'model5w')

    N=blkdiag2(model.hyp.sigma_v(1)^2*K_squared_diag,model.hyp.sigma_v(2)^2*K_squared_diag)...
     +blkdiag2(model.hyp.sigma_v(3)^2*In,model.hyp.sigma_v(4)^2*In);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*K_squared_diag,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_squared_diag,'sparse');
    N_grad{3}=blkdiag2(model.hyp.sigma_v(3)^2*In,Zn,'sparse');
    N_grad{4}=blkdiag2(Zn,model.hyp.sigma_v(4)^2*In,'sparse');
    
else

    error(['Noise option not permissible: ' model.noise]);
end
