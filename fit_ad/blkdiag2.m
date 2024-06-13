function B=blkdiag2(varargin)
%% Fast block diagonal
%
% Code same as inbuildt blkdiag except that ...
% ... blkdiag has costly overhead checks
% ... blkdiag always gives out full matrices, even if inputs are sparse
%
% Inputs:
% X1,X2,X3,...: matrices to stack
% Add 'sparse' or 's' to end of input to preserve sparseness
%
% Outputs:
% B: block-diagonal matrix
%

%%

do_sparse=false;
n_blk=nargin;

if isstr(varargin{end})
    n_blk=n_blk-1;

    if strcmpi(varargin{end},'sparse') || strcmpi(varargin{end},'s') 
        do_sparse=true;
    end
end

% Call internal blkdiag

B=matlab.internal.math.blkdiag(varargin{1:n_blk});

if do_sparse==false

    % If sparse, make full
    if issparse(B)
        B=full(B);
    end

end
