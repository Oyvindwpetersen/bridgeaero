function model=modelopt_poly_fix(model)
%% Fix and check options
%
% Inputs:
% model: struct with GPR model
%
% Outputs:
% model: struct with GPR model
%

%%

if ~permopt(model.kernel,{'se'})
    error(['Kernel option not permissible: ' model.kernel]);
end

%%

model.nv=1;

% If scalars given, expand to vector
model.ini.sigma_v=expandcol(model.ini.sigma_v,model.nv);
model.lb.sigma_v=expandcol(model.lb.sigma_v,model.nv);
model.ub.sigma_v=expandcol(model.ub.sigma_v,model.nv);

model.ini.sigma=expandcol(model.ini.sigma,model.p+1);
model.lb.sigma=expandcol(model.lb.sigma,model.p+1);
model.ub.sigma=expandcol(model.ub.sigma,model.p+1);

model.ini.L=expandcol(model.ini.L,model.p+1);
model.lb.L=expandcol(model.lb.L,model.p+1);
model.ub.L=expandcol(model.ub.L,model.p+1);

%% Index

model.idx.sigma_v=[1:model.nv];
model.idx.sigma=model.idx.sigma_v(end)+[1:(model.p+1)];
model.idx.L=model.idx.sigma(end)+[1:(model.p+1)];

% model.idx.alpha=[];
% model.idx.abar=[];

%%

model=orderstruct(model);
