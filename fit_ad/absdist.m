function r=absdist(x1,x2,power)
%% Distance matrix in parameter space
%
% Inputs:
% x1: [n1*1] column vector
% x2: [n2*1] column vector
% power: exponent for distance measure 
%
% Outputs:
% r: [n1*n2] distance matrix
%
% r=[ |x1(1)-x2(1)|  |x1(1)-x2(2)|  ... |x1(1)-x2(n2)|  ]
%   [       .              .                  .         ]
%   [       .              .                  .         ]
%   [ |x1(n1)-x2(1)| |x1(n1)-x2(2)| ... |x1(n1)-x2(n2)| ]
%
%

%%

if ~isvector(x1) | ~isvector(x2)
    error('x1 and x2 must be vector');
end

if size(x1,2)>1
    error('Input must be column vector');
end

if size(x2,2)>1
    error('Input must be column vector');
end

r=abs(x1-x2.');

if nargin>2
    r=r.^power;
end
