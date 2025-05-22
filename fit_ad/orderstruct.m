function model_out=orderstruct(model)
%% Reorder field of the model struct
%
% Inputs:
% model: struct
%
% Outputs:
% model_out: same struct but with fields reordered
%

%%

% Fields in preferable order
preferred = {'basis','kernel','noise',...
    'nd','na','nv',...
    'delta_d',...
    'ini','lb','ub',...
    'idx','hyp',...
};

fn=fieldnames(model);

new_ord=[];
for k=1:length(preferred)
    new_ord(k)=NaN;
    for j=1:length(fn)
        if strcmpi(fn{j},preferred{k})
            new_ord(k)=j;
            break
        end
    end
end

% Index of preferred fields not found in the struct
idx_not_found=find(isnan(new_ord));

% Remove them
new_ord=new_ord(~isnan(new_ord));

% If the struct has additional fields, add them to the end
extra=setdiff([1:length(fn)],new_ord);
new_ord=[new_ord extra];

% Reorder
model_out=orderfields(model,new_ord);
