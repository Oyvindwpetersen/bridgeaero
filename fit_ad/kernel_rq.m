function [K,K_grad]=kernel_rq(X1,X2,sigma,L,alpha)
%% Rational quadratic kernel
%
% Inputs:
% X1: [n1,1] matrix with inputs as columns
% X2: [n2,1] matrix with inputs as columns
% sigma: magnitude parameter
% L: length scale
% alpha: power
%
% Outputs:
% K: covariance matrix
% K_grad: cell with gradients
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

    j=1;
    r=absdist(X1(:,j),X2(:,j));
    K=sigma^2.*(ones(size(r))+r.^2./(2*alpha*L^2)).^(-alpha);

    if nargout>1
        K_grad{1}=2*sigma/sigma^2*K;
        K_grad{2}=r.^2/L^3.*sigma^2.*(ones(size(r))+r.^2./(2*alpha*L^2)).^(-alpha-1);

        tmp1=r.^2+(2*alpha*L^2);
        tmp2=ones(size(r))+r.^2./(2*alpha*L^2);

        K_grad{3}=sigma^2.*tmp2.^(-alpha).*(r.^2-tmp1.*log(tmp2))./tmp1;

    end

end
