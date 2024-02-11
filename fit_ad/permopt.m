function accept=permopt(str_check,str_perm)
%% Check if string input is among permissible string arguments
%
% Inputs:
% str_check: string to check
% str_perm: cell with strings of permissible string arguments
%
% Outputs:
% matched: true/false logic
%

%%

if ~iscell(str_perm)
    str_perm={str_perm};
end

for k=1:length(str_perm)
    str_match(k)=strcmpi(str_check,str_perm{k});
end

if any(str_match)
    accept=true;
else
    accept=false;
end