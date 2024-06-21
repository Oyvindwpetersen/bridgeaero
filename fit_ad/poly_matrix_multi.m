function [T_glob,x_uni]=poly_matrix_multi(data_matrix,p)
%% Regression matrix in fit of aerodynamic derivatives using rational functions for multiple x
%
% Inputs:
% data_matrix: [M,2] matrix with K and x as columns
% p: polynomial order
%
% Outputs:
% T_glob: [M*nx,p+1] regression matrix
% x_uni: [nx,1] vector with unique x
%
% data_matrix=
% [ alpha1   x1 ]
% [ alpha2   x2 ]
% [ alpha3   x3 ]
% [   .      .  ]
%
% T_glob=
% [T(alpha1,x1)                   ]
% [             T(alpha2,x2)      ]
% [                           ... ]
%

%%

% Find unique x
[x_uni,~,ic]=unique(data_matrix(:,2));
nx=length(x_uni);

% Create bins of index for each unique x
bins=cell(1,nx);
for k=1:length(x_uni)
    bins{k}=find(ic==k);
end

poly_matrix =@(alpha,p) [alpha(:)].^([0:p]);

% T-matrix for each unique x
T_blk=cell(nx,1);

for k=1:nx
    alpha_vec=data_matrix(bins{k},1);

    [T_blk{k}]=poly_matrix(alpha_vec,p);

end
T_glob=[blkdiag2(T_blk{:},'sparse')];

