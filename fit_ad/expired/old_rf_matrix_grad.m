% function [Dr_grad,Di_grad]=rf_matrix_grad(d,K)

%% Regression matrix in fit of aerodynamic derivatives using rational functions
%
% Inputs:
% d: vector with poles
% K: vector with frequencies
%
% Outputs:
% Dr_grad: gradient of real part regression matrix with respect to poles
% Di_grad: gradient of imaginary part regression matrix with respect to poles
%

%%

% Ensure column
K=K(:);

%%

Dr_grad=cell(length(d),1);
Di_grad=cell(length(d),1);

for k=1:length(d)

    val=-2*d(k)*K.^2./(d(k).^2+K.^2).^2;
    Dr_grad{k}=sparse(1:length(K),(2+k)*ones(length(K),1),val,length(K),2+length(d));

    val=K.*(K.^2-d(k).^2)./(d(k).^2+K.^2).^2;
    Di_grad{k}=sparse(1:length(K),(2+k)*ones(length(K),1),val,length(K),2+length(d));

end


