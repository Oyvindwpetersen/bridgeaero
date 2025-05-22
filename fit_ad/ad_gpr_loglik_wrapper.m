function [neg_logL,neg_logL_grad,neg_logL_hess,neg_logL_grad_modelfit,neg_logL_grad_complexity]=ad_gpr_loglik_wrapper(test_matrix,yr,yi,model,theta)

%% Calculation of log likelihood and gradient (wrapper)
%
% Inputs:
% test_matrix: [M,2] matrix with K and x as columns
% yr: [M,1] vector with K^2*AD_stiffness
% yi: [M,1] vector with K^2*AD_damping
% model: struct with GPR model
% theta: vector with hyperparameters
%
% Outputs:
% neg_logL: -log(L)
% neg_logL_grad: gradient of -log(L) wrt. hyperparameters
%

%%

cell_check=[iscell(test_matrix) iscell(yr) iscell(yi) iscell(model)];

if all(cell_check)
    % OK
    cell_out=true;
elseif ~any(cell_check)
    cell_out=false;
    test_matrix={test_matrix};
    yi={yi};
    yi={yi};
    model={model};
else
    error('All or none of inputs must be cell');
end

n1=size(test_matrix,1);
n2=size(test_matrix,2);

%% Assign hyperparameters to model

for idx1=1:size(model,1)
    for idx2=1:size(model,2)

        theta_loc=theta(model{idx1,idx2}.idx.glob);
        model{idx1,idx2}=hyp2model(model{idx1,idx2},theta_loc);

    end
end

%% Log likelihood for each DOF

neg_logL_grad=zeros(length(theta),1);
neg_logL=zeros(size(model,1),size(model,2));

neg_logL_grad_modelfit=zeros(length(theta),1);
neg_logL_grad_complexity=zeros(length(theta),1);

for idx1=1:size(model,1)
    for idx2=1:size(model,2)

        %
        [neg_logL_loc,neg_logL_grad_loc,neg_logL_grad_modelfit_loc,neg_logL_grad_complexity_loc]=ad_gpr_loglik(test_matrix{idx1,idx2},yr{idx1,idx2},yi{idx1,idx2},model{idx1,idx2});

        % Gradient contribution assignment
        neg_logL_grad(model{idx1,idx2}.idx.glob,1)=neg_logL_grad(model{idx1,idx2}.idx.glob,1)+neg_logL_grad_loc;

        % Gradient contribution assignment
        neg_logL_grad_modelfit(model{idx1,idx2}.idx.glob,1)=neg_logL_grad_modelfit(model{idx1,idx2}.idx.glob,1)+neg_logL_grad_modelfit_loc;
        neg_logL_grad_complexity(model{idx1,idx2}.idx.glob,1)=neg_logL_grad_complexity(model{idx1,idx2}.idx.glob,1)+neg_logL_grad_complexity_loc;

        % Contribution for each DOF
        neg_logL(idx1,idx2)=neg_logL_loc;

    end
end

% Sum all contributions
neg_logL=sum(sum(neg_logL));

%%

if nargout>2

    for idx=1:length(theta)

        theta_perturb=theta;

        theta_perturb(idx)=theta_perturb(idx)*1.01;

        [~,neg_logL_grad_perturb]=ad_gpr_loglik_wrapper(test_matrix,yr,yi,model,theta_perturb);

        neg_logL_hess(:,idx)=(neg_logL_grad-neg_logL_grad_perturb)./theta_perturb(idx)*0.01;

    end

else
    neg_logL_hess=[];
end


