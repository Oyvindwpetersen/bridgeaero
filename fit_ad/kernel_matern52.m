function [K,K_grad]=kernel_matern52(X1,X2,sigma,L)
%% Matern 5/2 kernel
%
% Inputs:
% X1: [n1,1] matrix with inputs as columns
% X2: [n2,1] matrix with inputs as columns
% sigma: magnitude parameter
% L: length scale
%
% Outputs:
% K: covariance matrix
% K_grad: cell with gradients
%
% K=sigma^2*prod{ (1+sqrt(5)*r_i/L_i+5*r_i^2/(3*L_i))*exp(sqrt(5)*r_i/L_i) }
%
% dK/dsigma=2*sigma*prod{ (1+sqrt(5)*r_i/L_i+5*r_i^2/(3*L_i))*exp(sqrt(5)*r_i/L_i) }
%
% dK/dL_n=5/3*(r_n^2/L_n^4)*exp(sqrt(5)*r_n/L_n)*sigma^2*prod_not_n{ (1+sqrt(5)*r_i/L_i+5*r_i^2/(3*L_i))*exp(sqrt(5)*r_i/L_i) }
%

%% Check

if size(X1,2)~=size(X2,2)
    s1=size(X1,2)
    s2=size(X2,2)
    error('X1 and X2 must have equal amount of columns');
end

if size(X1,2)>1
    error('Multiple dimension x not supported yet');
end

%% Case single input

if size(X1,2)==1

    r=absdist(X1(:,1),X2(:,1));
    K=sigma^2.*(ones(size(r))+sqrt(5)*r/L+5/3*(r/L).^2).*exp(-sqrt(5)*r/L);

    if nargout>1
        K_grad{1}=2*sigma/sigma^2*K;
        K_grad{2}=1./L*(-5*sqrt(5)/3*(r/L).^2-5/3*(r/L)).*exp(-sqrt(5)*r/L);
        
    end

end