function [data_matrix_sorted,x_uni,T,T_inv]=sortmatrix(data_matrix)
%% Predict using GPR model
%
% Inputs:
% data_matrix: [N,2] matrix with K and x as columns
%
% Outputs:
% data_matrix_sorted: [N,2] matrix with K and x as columns, sorted
% x_uni: unique values for x
% T: transformation matrix, y=T*y_sorted
% T_inv: inverse transformation matrix, y_sorted=T_inv*y

%%

% Find unique x
[x_uni,~,ic]=unique(data_matrix(:,2));
nx=length(x_uni);

bins = accumarray(ic, (1:numel(ic))', [], @(v){v}).';

% Create transformation matrix
idx_row=[1:size(data_matrix,1)].';
idx_col=vertcat(bins{:});
val=ones(size(idx_row));
T_inv=sparse(idx_row,idx_col,val,size(data_matrix,1),size(data_matrix,1));

T=inv(T_inv);

data_matrix_sorted=T_inv*data_matrix;

