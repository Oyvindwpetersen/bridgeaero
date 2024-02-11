function b=expandcol(a,n)
%% Expand scalar into column vector with length n
%
% Inputs:
% a: scalar or vector
% n: desired length of vector
%
% Outputs:
% b: output vector

%%

b=a(:);

if length(b)==1
    b=b*ones(n,1);
end
