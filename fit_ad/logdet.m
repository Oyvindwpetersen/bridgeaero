function [log_det_A,L]=logdet(A,method)

%% Logarithmic of matrix determinant
%
% For matrices that area large or close to singular, the det() function can give inaccurate results
% The log of the determinant can be calculated using the singular value or cholesky decomposition 
% The chol option requires a positive definite matrix
%
% Inputs:
% A: square matrix with inputs as columns
% method: 'det', 'svd', or 'chol'
%
% Outputs:
% log_det_A: log(det(A))
%

%%

[n1,n2]=size(A);

if n1~=n2
    error('Matrix must be square');
end

if isscalar(A)
    log_det_A=log(A);
    
    if isinf(log_det_A)
        warning('log(det()) is Inf');
    end

    return
end

if nargin==1
    method='svd';
end

if strcmpi(method,'det')
    log_det_A=log(det(A));

    if isinf(log_det_A)
        warning('log(det()) is Inf, consider using another option');
    end
end

if strcmpi(method,'svd')
    [u,s,v]=svd(A);
    log_det_A=sum(log(diag(s)));
end

L=[];
if strcmpi(method,'chol')
    L=chol(A,'lower');
    log_det_A=2*sum(log(diag(L)));
end