function [yp,std_y_p,std_y_p_obs,ap,cov_ap]=sc_gpr_pred(test_matrix,pred_matrix,y,model)
%% Predict using GPR model
%
% Inputs:
% test_matrix: [M,2] matrix with alpha and x as columns
% pred_matrix: [M,2] matrix with alpha and x as columns
% y: [M,1] vector with outputs
% model: struct with GPR model
%
% Outputs:
% yp: [N,1] vector with K^2*AD_stiffness
% std_y_p: [N,1] vector with SD of yr_p
% std_y_p_obs: [N,1] vector with SD of yr_p (+noise)
% ap: predicted a-coefficients
% cov_ap: covariance of a-coefficients
%

%%

cell_check=[iscell(test_matrix) iscell(pred_matrix) iscell(y) iscell(model)];

if all(cell_check)
    % OK
    cell_out=true;
elseif ~any(cell_check)
    cell_out=false;
    test_matrix={test_matrix};
    pred_matrix={pred_matrix};
    y={y};
    model={model};
else
    error('All or none of inputs must be cell');
end

n1=size(test_matrix,2);

%% Predict for all

for idx1=1:n1

        [yp{idx1},std_y_p{idx1},std_y_p_obs{idx1},ap{idx1},cov_ap{idx1}]=...
            sc_gpr_pred_single(test_matrix{idx1},pred_matrix{idx1},y{idx1},model{idx1});

end

%%

if cell_out==false
    yp=yp{1,1};
    std_y_p=std_y_p{1,1};
    std_y_p_obs=std_y_p_obs{1,1};
    ap=ap{1,1};
    cov_ap=cov_ap{1,1};
end

