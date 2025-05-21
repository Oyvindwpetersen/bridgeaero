function [yr,yi]=rf2adks(a,d,K)
%% Convert rational function to AD*K^2
%
% Inputs:
% a: 3d-matrix with a-coefficients
% d: vector with poles
% K: vector with frequencies
%
% Outputs:
% yr: stiffness ADs times K^2
% yi: damping ADs times K^2
%
%%

[Dr,Di]=rf_matrix(d,K);

yr=zeros(size(a,1),size(a,2),length(K));
yi=zeros(size(a,1),size(a,2),length(K));

for k1=1:size(a,1)
    for k2=1:size(a,2)

        a_vec=squeeze(a(k1,k2,:));

        yr(k1,k2,:)=Dr*a_vec;
        yi(k1,k2,:)=Di*a_vec;

    end
end
