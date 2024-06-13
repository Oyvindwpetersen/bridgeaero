function model=hyp2model(model,theta)
%% Assign hyperparameters to model
%
% Inputs:
% model: struct with GPR model
% theta: vector with hyperparameters
%
% Outputs:
% model: struct with GPR model (updated)
%

%%

theta=theta(:);

model.hyp.d=theta(model.idx.d);
model.hyp.sigma_v=theta(model.idx.sigma_v);
model.hyp.sigma=theta(model.idx.sigma);
model.hyp.L=theta(model.idx.L);
model.hyp.alpha=theta(model.idx.alpha);
model.hyp.abar=theta(model.idx.abar);

