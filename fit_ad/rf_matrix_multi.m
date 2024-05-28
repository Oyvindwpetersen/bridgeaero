function [D_glob,x_uni,D_glob_grad]=rf_matrix_multi(data_matrix,d)
%% Regression matrix in fit of aerodynamic derivatives using rational functions
% for multiple x
%
% Inputs:
% data_matrix: [M,2] matrix with K and x as columns
% d: [1,nd] vector with poles
%
% Outputs:
% D_glob: [2M*nx,nd+2] regression matrix
% x_uni: [nx,1] vector with unique x
% D_glob_grad: cell with gradients of regression matrix with respect to poles
%
% data_matrix=
% [ K1   x1 ]
% [ K2   x2 ]
% [ K3   x3 ]
% [ .    .  ]
%
% D_glob=
% [Dr(x1,K1)                ]
% [          Dr(x2,K2)      ]
% [                     ... ]
% [Di(x1,K1)                ]
% [          Di(x2,K2)      ]
% [                     ... ]
%

%%

if nargout>2
    return_grad=true;
else
    return_grad=false;
end

% Find unique x
[x_uni,~,ic]=unique(data_matrix(:,2));
nx=length(x_uni);

% Create bins of index for each unique x
bins=cell(1,nx);
for k=1:length(x_uni)
    bins{k}=find(ic==k);
end

% D-matrix for each unique x
Dr_blk=cell(nx,1);
Di_blk=cell(nx,1);
Dr_grad_blk=cell(nx,length(d));
Di_grad_blk=cell(nx,length(d));

for k=1:nx
    K_vec=data_matrix(bins{k},1);

    if return_grad
        [Dr_blk{k},Di_blk{k},Dr_grad_blk(k,:),Di_grad_blk(k,:)]=rf_matrix(d,K_vec);
    else
        [Dr_blk{k},Di_blk{k}]=rf_matrix(d,K_vec);
    end

end
D_glob=[blkdiag2(Dr_blk{:},'full') ; blkdiag2(Di_blk{:},'full')];

if return_grad
    D_glob_grad=cell(length(d),1);
    for j=1:length(d)
        
        D_glob_grad{j}=[blkdiag2(Dr_grad_blk{:,j},'full') ; blkdiag2(Di_grad_blk{:,j},'full')];

    end
end

