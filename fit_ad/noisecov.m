function [N,N_grad]=noisecov(data_matrix,model)

%% Noise matrix
%
% Inputs:
% data_matrix: [n,2] matrix with inputs as columns
% model: struct with model valuesparameter
%
% Outputs:
% N: [2n,2n] noise matrix
% N_grad: gradient of N wrt. sigma_v
%

%%

n=size(data_matrix,1);
In=speye(n);
Zn=In*0;

K=data_matrix(:,1);

if strcmpi(model.noise,'model3') | strcmpi(model.noise,'model3w') | strcmpi(model.noise,'model4') | strcmpi(model.noise,'model4w')
    % K_diag=diag(K);
    K_diag=spdiags(K(:),0,numel(K),numel(K));
    K_squared_diag=K_diag.^2;
end

% if strcmpi(model.noise,'model4') | strcmpi(model.noise,'model4w')
%     K_squared_diag=diag(K.^2);
% end

% if strcmpi(model.noise,'model_a')
%     [D_glob,xt_uni,D_glob_grad]=rf_matrix_multi(data_matrix,model.hyp.d);
%     nx=length(xt_uni);
%     [S]=restack_a(model.na,length(xt_uni));
% end

%% Noise matrix

if strcmpi(model.noise,'model1')

    N=model.hyp.sigma_v^2*eye(n);

    N_grad{1}=2*model.hyp.sigma_v(1)*eye(n);

elseif strcmpi(model.noise,'model2')

    N=blkdiag2(model.hyp.sigma_v(1)^2*In,model.hyp.sigma_v(2)^2*In);

    N_grad{1}=blkdiag2(model.hyp.sigma_v(1)^2*In,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,model.hyp.sigma_v(2)^2*In,'sparse');

elseif strcmpi(model.noise,'model3')

    N=blkdiag2(model.hyp.sigma_v(1)^2*In,model.hyp.sigma_v(2)^2*K_diag);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*In,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_diag,'sparse');

elseif strcmpi(model.noise,'model3w')

    N=blkdiag2(model.hyp.sigma_v(1)^2*In,model.hyp.sigma_v(2)^2*K_diag)...
      +blkdiag2(Zn,model.hyp.sigma_v(3)^2*In);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*In,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_diag,'sparse');
    N_grad{3}=blkdiag2(Zn,2*model.hyp.sigma_v(3)*In,'sparse');

elseif strcmpi(model.noise,'model4')

    N=blkdiag2(model.hyp.sigma_v(1)^2*K_squared_diag,model.hyp.sigma_v(2)^2*K_squared_diag);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*K_squared_diag,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_squared_diag,'sparse');
    
elseif strcmpi(model.noise,'model4w')

    N=blkdiag2(model.hyp.sigma_v(1)^2*K_squared_diag,model.hyp.sigma_v(2)^2*K_squared_diag)...
     +blkdiag2(model.hyp.sigma_v(3)^2*In,model.hyp.sigma_v(4)^2*In);
    
    N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*K_squared_diag,Zn,'sparse');
    N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*K_squared_diag,'sparse');
    N_grad{3}=blkdiag2(model.hyp.sigma_v(3)^2*In,Zn,'sparse');
    N_grad{4}=blkdiag2(Zn,model.hyp.sigma_v(4)^2*In,'sparse');
    
% elseif strcmpi(model.noise,'model_a')

    % N=blkdiag2(model.hyp.sigma_v(1)^2*In,model.hyp.sigma_v(2)^2*In);
    %   +D_glob*S*model.hyp.sigma_v(3)^2*S.'*D_glob.';
    % 
    % N_grad{1}=blkdiag2(2*model.hyp.sigma_v(1)*In,Zn,'sparse');
    % N_grad{2}=blkdiag2(Zn,2*model.hyp.sigma_v(2)*In,'sparse');
    % N_grad{3}=D_glob*S*2*model.hyp.sigma_v(3)*S.'*D_glob.';
else
    error(['Noise option not permissible: ' model.noise]);
end
