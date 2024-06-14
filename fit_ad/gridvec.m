function x=gridvec(varargin)
%% Create grid based on a number of vectors
%
% Inputs:
% x1,x2,x3...: vectors
%
% Outputs:
% x: [n1*n2*n3...,N] combination matrix 
%

%% Check if x is a matrix of vectors

% e.g. x_in=[1:10 ; 1:10];

input_cell={};

% Cell with vectors
if nargin==1 & iscell(varargin{1})
    input_cell_temp=varargin{1};
    n_input=length(input_cell_temp);
    for k=1:n_input
        input_cell{k}(:,1)=input_cell_temp{k};
    end
end

% Multiple vectors
if nargin>1 & isnumeric(varargin{1})
    n_input=nargin;
    for k=1:n_input
        input_cell{k}(:,1)=varargin{k};
    end
end


%%

% input_cell={[1:10] [50 60 70] [0.1 0.2]}
for i=1:length(input_cell)
    siz(i)=numel(input_cell{i});
end

out=zeros(prod(siz),length(input_cell));

nout=length(input_cell);

for i=1:length(input_cell)
        x=input_cell{i};
        s=ones(1,nout); 
        s(i)=numel(x);
        x=reshape(x,s);
        s=siz;
        s(i)=1;
        tmp=repmat(x,s);
        out(:,i)=tmp(:);
end

x=out;
