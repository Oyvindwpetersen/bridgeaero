function model=hypdefault2model_poly(test_matrix,y,model)
%% Fix and check options
%
% Inputs:
% model: struct with GPR model
%
% Outputs:
% model: struct with GPR model
%
%%

y_std=std([y]);
L_char=max(test_matrix(:,2))-min(test_matrix(:,2));
% y_mean=mean([yr_t;yr_i]);

if isempty(model.ini.sigma_v); model.ini.sigma_v=expandcol(y_std*1e-1,model.nv); end
if isempty(model.lb.sigma_v); model.lb.sigma_v=expandcol(y_std*1e-3,model.nv); end
if isempty(model.ub.sigma_v); model.ub.sigma_v=expandcol(y_std*1e1,model.nv); end

if isempty(model.ini.sigma); model.ini.sigma=expandcol(1*1e-1,model.p+1); end
if isempty(model.lb.sigma); model.lb.sigma=expandcol(1*1e-6,model.p+1); end
if isempty(model.ub.sigma); model.ub.sigma=expandcol(1*1e2,model.p+1); end

if isempty(model.ini.L); model.ini.L=expandcol(L_char*1e0,model.p+1); end
if isempty(model.lb.L); model.lb.L=expandcol(L_char*1e-1,model.p+1); end
if isempty(model.ub.L); model.ub.L=expandcol(L_char*1e1,model.p+1); end

