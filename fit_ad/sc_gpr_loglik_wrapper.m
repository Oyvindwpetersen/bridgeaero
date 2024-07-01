function [neg_logL,neg_logL_grad]=sc_gpr_loglik_wrapper(test_matrix,y,model,theta)

%% Calculation of log likelihood and gradient (wrapper)
%
% Inputs:
% test_matrix: [M,2] matrix with alpha and x as columns
% y: [M,1] vector with outputs
% model: struct with GPR model
% theta: vector with hyperparameters
%
% Outputs:
% neg_logL: -log(L)
% neg_logL_grad: gradient of -log(L) wrt. hyperparameters
%

%% Assign hyperparameters to model

model=hyp2model_poly(model,theta);

[neg_logL,neg_logL_grad]=sc_gpr_loglik(test_matrix,y,model);