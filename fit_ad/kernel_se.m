function [K,K_grad]=kernel_se(X1,X2,sigma,L)
%% Squared exponential kernel
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
% K=sigma^2*prod{ exp(-0.5*r_i^2/L_i^2) }
%
% dK/dsigma=2*sigma*prod{ exp(-0.5*r_i^2/L_i^2 }
%
% dK/dL_n=(r_n^2/L_n^3)*sigma^2*prod{ exp(-0.5*r_i^2/L_i^2) }
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

    n=1;
    r_squared=absdist(X1(:,n),X2(:,n),2);

    K=sigma^2*exp(-0.5*r_squared/L(n)^2);

    if nargout>1
        K_grad{1}=2*sigma/sigma^2*K;
        K_grad{2}=r_squared./L(n)^3.*K;

    end
end
