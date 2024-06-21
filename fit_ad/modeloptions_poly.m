function model=modeloptions_poly(varargin)
%% Create options for model: kernel, noise, hyperparameters, bounds, ...
%
% Inputs:
% given below
%
% Outputs:
% model: struct with GPR model
%

%% Parse inputs

p=inputParser;
addParameter(p,'kernel','se',@ischar)

addParameter(p,'p',[2],@isnumeric)

addParameter(p,'sigma_v0',[],@isnumeric)
addParameter(p,'sigma0',[],@isnumeric)
addParameter(p,'L0',[],@isnumeric)

addParameter(p,'sigma_v_lb',[],@isnumeric)
addParameter(p,'sigma_lb',[],@isnumeric)
addParameter(p,'L_lb',[],@isnumeric)

addParameter(p,'sigma_v_ub',[],@isnumeric)
addParameter(p,'sigma_ub',[],@isnumeric)
addParameter(p,'L_ub',[],@isnumeric)


parse(p,varargin{:})

%% Assign to model

model.p=p.Results.p;

model.kernel=p.Results.kernel;

model.ini.sigma_v=p.Results.sigma_v0;
model.ini.sigma=p.Results.sigma0;
model.ini.L=p.Results.L0;

model.lb.sigma_v=p.Results.sigma_v_lb;
model.lb.sigma=p.Results.sigma_lb;
model.lb.L=p.Results.L_lb;

model.ub.sigma_v=p.Results.sigma_v_ub;
model.ub.sigma=p.Results.sigma_ub;
model.ub.L=p.Results.L_ub;

%%

model=modeloptions_poly_fix(model);




