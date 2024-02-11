function [K,K_grad]=kernel_rq(X1,X2,sigma,L,alpha)
%% Rational quadratic kernel
%
% Inputs:
% X1: [n1,m] matrix with inputs as columns
% X2: [n2,m] matrix with inputs as columns
% sigma: magnitude parameter
% L: [1,m] vector with length scales
% alpha: [1,m] vector with powers
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
%% Case X1 and X2 are matrices


%% Gradient

% if nargout>1
%
%     n=1;
%     K_grad{1}=2*sigma/sigma^2*K;
%
%     j=1;
%     r=absdist(X1(:,j),X2(:,j));
%
%
%
%
%
%
%
%     K_grad{2}=absdist(X1(:,n),X2(:,n),2)./L(n)^4.* (K./ ((ones(size(r))+sqrt(5)*r/L(j)+5/3*(r/L(j)).^2).*exp(sqrt(5)*r/L(j))) );
%
%     for n=2:length(L)
%
%         K_grad{n+2}=absdist(X1(:,n),X2(:,n),2)./L(n)^3.*K; %K.*exp(-0.5*absdist(X1(:,n),X2(:,n),2)/L(n)^2);
%
%     end
% end
%


%%

% for j=2:length(L)
% 
%     r=absdist(X1(:,j),X2(:,j));
% 
%     K=K.*(ones(size(r))+sqrt(5)*r/L(j)+5/3*(r/L(j)).^2).*exp(sqrt(5)*r/L(j))
% 
% end
