function model=hypdefault2model(test_matrix,yr_t,yr_i,model)
%% Fix and check options
%
% Inputs:
% model: struct with GPR model
%
% Outputs:
% model: struct with GPR model
%
%%

y_std=std([yr_t;yr_i]);
L_char=max(test_matrix(:,2))-min(test_matrix(:,2));
y_mean=mean([yr_t;yr_i]);

if isempty(model.ini.sigma_v); model.ini.sigma_v=expandcol(y_std*1e-1,model.nv); end
if isempty(model.lb.sigma_v); model.lb.sigma_v=expandcol(y_std*1e-3,model.nv); end
if isempty(model.ub.sigma_v); model.ub.sigma_v=expandcol(y_std*1e0,model.nv); end

if isempty(model.ini.sigma); model.ini.sigma=expandcol(y_std*1e0,model.na); end
if isempty(model.lb.sigma); model.lb.sigma=expandcol(y_std*1e-2,model.na); end
if isempty(model.ub.sigma); model.ub.sigma=expandcol(y_std*1e2,model.na); end

if isempty(model.ini.L); model.ini.L=expandcol(L_char*1e0,model.na); end
if isempty(model.lb.L); model.lb.L=expandcol(L_char*1e-1,model.na); end
if isempty(model.ub.L); model.ub.L=expandcol(L_char*1e1,model.na); end

if ~isempty(model.idx.alpha)
    if isempty(model.ini.alpha); model.ini.alpha=expandcol(1e0,model.na); end
    if isempty(model.lb.alpha); model.lb.alpha=expandcol(1e-1,model.na); end
    if isempty(model.ub.alpha); model.ub.alpha=expandcol(1e1,model.na); end
end

if ~isempty(model.idx.abar)
    if isempty(model.ini.abar); model.ini.abar=expandcol(y_mean*1e0,model.na); end
    if isempty(model.lb.abar); model.lb.abar=expandcol(y_mean*1e-1,model.na); end
    if isempty(model.ub.abar); model.ub.abar=expandcol(y_mean*1e1,model.na); end
end