function model=modeloptions(varargin)

%% Create options for model: kernel, noise, hyperparameters, bounds, ...
%
% Inputs:
% variable inputs 
%
% Outputs:
% model: struct
%

%% Parse inputs

p=inputParser;
addParameter(p,'nd',3,@isnumeric)
addParameter(p,'kernel','se',@ischar)
addParameter(p,'noise','model1',@ischar)
addParameter(p,'basis','zero',@ischar)

addParameter(p,'d0',[],@isnumeric)
addParameter(p,'sigma_v0',[],@isnumeric)
addParameter(p,'sigma0',[],@isnumeric)
addParameter(p,'L0',[],@isnumeric)
addParameter(p,'alpha0',[],@isnumeric)
addParameter(p,'abar0',[],@isnumeric)

addParameter(p,'d_lb',[0.1],@isnumeric)
addParameter(p,'sigma_v_lb',[],@isnumeric)
addParameter(p,'sigma_lb',[],@isnumeric)
addParameter(p,'L_lb',[],@isnumeric)
addParameter(p,'alpha_lb',[],@isnumeric)
addParameter(p,'abar_lb',[],@isnumeric)

addParameter(p,'d_ub',[3],@isnumeric)
addParameter(p,'sigma_v_ub',[],@isnumeric)
addParameter(p,'sigma_ub',[],@isnumeric)
addParameter(p,'L_ub',[],@isnumeric)
addParameter(p,'alpha_ub',[],@isnumeric)
addParameter(p,'abar_ub',[],@isnumeric)

addParameter(p,'delta_d',0.05,@isnumeric)

parse(p,varargin{:})

%% Assign to model

model.nd=p.Results.nd;
model.kernel=p.Results.kernel;
model.noise=p.Results.noise;
model.basis=p.Results.basis;

model.ini.d=p.Results.d0;
model.ini.sigma_v=p.Results.sigma_v0;
model.ini.sigma=p.Results.sigma0;
model.ini.L=p.Results.L0;
model.ini.alpha0=p.Results.alpha0;
model.ini.abar0=p.Results.abar0;

model.lb.d=p.Results.d_lb;
model.lb.sigma_v=p.Results.sigma_v_lb;
model.lb.sigma=p.Results.sigma_lb;
model.lb.L=p.Results.L_lb;
model.lb.alpha=p.Results.alpha_lb;
model.lb.abar=p.Results.abar_lb;

model.ub.d=p.Results.d_ub;
model.ub.sigma_v=p.Results.sigma_v_ub;
model.ub.sigma=p.Results.sigma_ub;
model.ub.L=p.Results.L_ub;
model.ub.alpha=p.Results.alpha_ub;
model.ub.abar=p.Results.abar_ub;

model.delta_d=p.Results.delta_d;

%%

model=modeloptions_fix(model);




