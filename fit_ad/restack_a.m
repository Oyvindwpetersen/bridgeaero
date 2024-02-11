function [Sa,Sa_inv]=restack_a(na,nx)
%% Repopulation matrix to restack a coefficients
%
% Inputs:
% na: number of a coefficients
% nx: length of x
%
% Outputs:
% Sa: [na*nx,na*nx] repopulation matrix
% Sa_inv: [na*nx,na*nx] inverse repopulation matrix
%
% a=[a1(x1) a2(x1) ... aN(x1) , a1(x2) a2(x2) ... aN(x2) , a1(x3) a2(x3) ... aN(x3)];
%
% a_restack=[a1(x1) a1(x2) .... a1(xN) , a2(x1) a2(x2) ... a2(xN) , a3(x1) a3(x2) ... a3(xN) ];
%
% a=Sa*a_restack;
%

%%

val=ones(na*nx,1);

s=[1:length(val)];

col=reshape((reshape(s,[na nx])).',[],1);

row=[1:(na*nx)].';

Sa=sparse(col,row,val);

if nargout>1
    Sa_inv=inv(Sa);
end

%% Test data

% nx=3

% a_in=[ [1:10] [31:40] [101:110]].';

% [Sa,Sa_inv]=restack_a(na,nx)

% Sa*a_in