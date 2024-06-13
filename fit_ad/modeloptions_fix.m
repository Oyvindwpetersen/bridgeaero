function model=modelopt_fix(model)
%% Fix and check options
%
% Inputs:
% model: struct with GPR model
%
% Outputs:
% model: struct with GPR model
%

%%
if ~permopt(model.kernel,{'se' 'matern52' 'rq'})
    error(['Kernel option not permissible: ' model.kernel]);
end

if ~permopt(model.noise,{'model1' 'model2' 'model3' 'model3w' 'model4' 'model4w'})
    error(['Noise option not permissible: ' model.noise]);
end

if ~permopt(model.basis,{'zero' 'constant'})
    error(['Basis option not permissible: ' model.basis]);
end

%%

if strcmpi(model.noise,'model1')
    model.nv=1;
elseif strcmpi(model.noise,'model2')
    model.nv=2;
elseif strcmpi(model.noise,'model3')
    model.nv=2;
elseif strcmpi(model.noise,'model3w')
    model.nv=3;
elseif strcmpi(model.noise,'model4')
    model.nv=2;
elseif strcmpi(model.noise,'model4w')
    model.nv=4;
elseif strcmpi(model.noise,'model_a')
    model.nv=NaN;
end

model.na=model.nd+2;

% 
if isempty(model.ini.d)
    model.ini.d=linspace(min(model.lb.d),max(model.ub.d),model.nd+2);
    model.ini.d=model.ini.d(2:end-1);
end

model.ini.d=model.ini.d(:);
model.lb.d=expandcol(model.lb.d,model.nd);
model.ub.d=expandcol(model.ub.d,model.nd);

% If scalars given, expand to vector
model.ini.sigma_v=expandcol(model.ini.sigma_v,model.nv);
model.lb.sigma_v=expandcol(model.lb.sigma_v,model.nv);
model.ub.sigma_v=expandcol(model.ub.sigma_v,model.nv);

model.ini.sigma=expandcol(model.ini.sigma,model.na);
model.lb.sigma=expandcol(model.lb.sigma,model.na);
model.ub.sigma=expandcol(model.ub.sigma,model.na);

model.ini.L=expandcol(model.ini.L,model.na);
model.lb.L=expandcol(model.lb.L,model.na);
model.ub.L=expandcol(model.ub.L,model.na);

model.ini.alpha0=expandcol(model.ini.alpha0,model.na);
model.lb.alpha=expandcol(model.lb.alpha,model.na);
model.ub.alpha=expandcol(model.ub.alpha,model.na);

model.ini.abar0=expandcol(model.ini.abar0,model.na);
model.lb.abar=expandcol(model.lb.abar,model.na);
model.ub.abar=expandcol(model.ub.abar,model.na);

if ~strcmpi(model.kernel,'rq')
    model.ini.alpha0=[];
    model.lb.alpha=[];
    model.ub.alpha=[];
end

if strcmpi(model.basis,'zero')
    model.ini.abar0=[];
    model.lb.abar=[];
    model.ub.abar=[];
end

%% Index

model.idx.d=[1:model.nd];
model.idx.sigma_v=model.idx.d(end)+[1:model.nv];
model.idx.sigma=model.idx.sigma_v(end)+[1:model.na];
model.idx.L=model.idx.sigma(end)+[1:model.na];

if strcmpi(model.kernel,'rq')
    model.idx.alpha=model.idx.L(end)+[1:model.na];
else
    model.idx.alpha=[];
end

if strcmpi(model.basis,'constant')
    model.idx.abar=model.idx.L(end)+[1:model.na];

    if strcmpi(model.kernel,'rq')
        model.idx.abar=model.idx.abar+model.na;
    end
else
    model.idx.abar=[];
end

%%

model=orderstruct(model);
