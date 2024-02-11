function check_options(model)

%%

required_fields={
    'nd' , 'num' ;
    'na' , 'num' ;
    'nv' , 'num' ;
    'kernel' , 'str' ;
    'noise' , 'str' ;
    'basis' ,  'str'
    };

%% Identify culprits

idx_missing=[];
idx_wrong=[];

for k=1:length(required_fields)

    if isfield(model,required_fields{k,1})

        tmp=getfield(model,required_fields{k,1});

        if strcmpi(required_fields{k,2},'num') & ~isnumeric(tmp)
            idx_wrong(end+1)=k;
        elseif strcmpi(required_fields{k,2},'str') & ~isstring(tmp)
            idx_wrong(end+1)=k;
        end

    else
        idx_missing(end+1)=k;
    end
end

%% Error

do_error=false;

if ~isempty(idx_missing)
    for k=1:length(idx_missing)
        warning(['options is missing the following field: ' required_fields{idx_missing(k),1}]);
    end
    do_error=true;
end

if ~isempty(idx_wrong)
    for k=1:length(idx_missing)
        warning(['options field is the wrong type: ' required_fields{idx_wrong(k),1}]);
    end
    do_error=true;
end

if do_error
    error('Missing or wrong options, cannot continue');
end