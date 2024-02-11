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

% Assume supplied string at end implies sparse (if sparse input)
if isstr(varargin{end})
    do_sparse=true;
else
    do_sparse=false;
end

% Call internal blkdiag
if do_sparse

    B=matlab.internal.math.blkdiag(varargin{1:end-1});

else

    B=matlab.internal.math.blkdiag(varargin{:});

    % If sparse, make full
    if issparse(B)
        B=full(B);
    end

end
