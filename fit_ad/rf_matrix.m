function [Dr,Di,Dr_grad,Di_grad]=rf_matrix(d,K)
%% Regression matrix in fit of aerodynamic derivatives using rational functions
%
% Inputs:
% d: vector with poles
% K: vector with frequencies
%
% Outputs:
% Dr: regression matrix for real part (K^2*AD_stiffness)
% Di: regression matrix for imaginary part (K^2*AD_damping)
%

%%

% Ensure column
K=K(:);

% Ensure row
d=d(:).';

D_real_qs=[ones(length(K),1) zeros(length(K),1)];
D_imag_qs=[zeros(length(K),1) K];

D_real_rf=K.^2./(d.^2+K.^2);
D_imag_rf=d.*K./(d.^2+K.^2);

Dr=[D_real_qs D_real_rf];
Di=[D_imag_qs D_imag_rf];

%% Gradient

if nargout>3

    Dr_grad=cell(length(d),1);
    Di_grad=cell(length(d),1);

    for k=1:length(d)

        row=[1:length(K)];
        col=(2+k)*ones(length(K),1);

        val=-2*d(k)*K.^2./(d(k).^2+K.^2).^2;
        Dr_grad{k}=sparse(row,col,val,length(K),2+length(d));

        val=K.*(K.^2-d(k).^2)./(d(k).^2+K.^2).^2;
        Di_grad{k}=sparse(row,col,val,length(K),2+length(d));

    end

end

