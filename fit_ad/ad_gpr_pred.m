function [yr_p,yi_p,std_yr_p,std_yi_p,std_yr_p_obs,std_yi_p_obs,ap,cov_ap]=ad_gpr_pred(test_matrix,pred_matrix,yr_t,yr_i,model)
%% Predict using GPR model
%
% Inputs:
% test_matrix: [M,2] matrix with K and x as columns
% pred_matrix: [N,2] matrix with K and x as columns
% yr_t: [M,1] vector with K^2*AD_stiffness
% yr_i: [M,1] vector with K^2*AD_damping
% model: struct with GPR model
%
% Outputs:
% yr_p: [N,1] vector with K^2*AD_stiffness
% yi_p: [N,1] vector with K^2*AD_damping
% std_yr_p: [N,1] vector with SD of yr_p
% std_yi_p: [N,1] vector with SD of yi_p
% std_yr_p_obs: [N,1] vector with SD of yr_p (+noise)
% std_yi_p_obs: [N,1] vector with SD of yr_p (+noise)
% ap: predicted a coefficients
%
%

%%

cell_check=[iscell(test_matrix) iscell(pred_matrix) iscell(yr_t) iscell(yr_i) iscell(model)];

if all(cell_check)
    % OK
    cell_out=true;
elseif ~any(cell_check)
    cell_out=false;
    test_matrix={test_matrix};
    pred_matrix={pred_matrix};
    yr_t={yr_t};
    yr_i={yr_i};
    model={model};
else
    error('All or none of inputs must be cell');
end

n1=size(test_matrix,1);
n2=size(test_matrix,2);

%% Predict ADs for a

for idx1=1:n1
    for idx2=1:n2

        [yr_p_loc,yi_p_loc,std_yr_p_loc,std_yi_p_loc,std_yr_p_obs_loc,std_yi_p_obs_loc,ap_loc,cov_ap_loc]=...
            ad_gpr_pred_single(test_matrix{idx1,idx2},pred_matrix{idx1,idx2},yr_t{idx1,idx2},yr_i{idx1,idx2},model{idx1,idx2});

        yr_p{idx1,idx2}=yr_p_loc;
        yi_p{idx1,idx2}=yi_p_loc;
        std_yr_p{idx1,idx2}=std_yr_p_loc;
        std_yi_p{idx1,idx2}=std_yi_p_loc;
        std_yr_p_obs{idx1,idx2}=std_yr_p_obs_loc;
        std_yi_p_obs{idx1,idx2}=std_yi_p_obs_loc;
        ap{idx1,idx2}=ap_loc;
        cov_ap{idx1,idx2}=cov_ap_loc;

    end
end

%%

if cell_out==false
    yr_p=yr_p{1,1};
    yi_p=yi_p{1,1};
    std_yr_p=std_yr_p{1,1};
    std_yi_p=std_yi_p{1,1};
    std_yr_p_obs=std_yr_p_obs{1,1};
    std_yi_p_obs=std_yi_p_obs{1,1};
    ap=ap{1,1};
    cov_ap=cov_ap{1,1};
end
