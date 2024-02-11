function [K,K_grad]=kernel_matern52(X1,X2,sigma,L)
%% Matern 5/2 kernel
%
% Inputs:
% X1: [n1,m] matrix with inputs as columns
% X2: [n2,m] matrix with inputs as columns
% sigma: magnitude parameter
% L: [1,m] vector with length scales
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
