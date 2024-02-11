function plotadsplitref(ha,K_ref,yr_ref,yi_ref,varargin)

%% Plot separate
%
% Inputs:
% ha: vect0r with handles to axes
% K_ref: vector with frequencies
% yr_ref: matrix with real part of transfer function (column for each x)
% yi_ref: matrix with imaginary part of transfer function (column for each x)
%
% Outputs:
% 

%%

p=inputParser;
addParameter(p,'color',[1 0.5 0],@isnumeric)
addParameter(p,'linestyle','--',@ischar)
addParameter(p,'linewidth',[1],@isnumeric)
addParameter(p,'displayname','Ref',@ischar)

parse(p,varargin{:})

color=p.Results.color;
linestyle=p.Results.linestyle;
linewidth=p.Results.linewidth;
displayname=p.Results.displayname;

%%

plotopt=struct();
plotopt.color=color;
plotopt.linestyle=linestyle;
plotopt.linewidth=linewidth;
plotopt.displayname=displayname;

for k=1:length(ha)/2

    axes(ha(k));
    h1(k)=plot(K_ref,yr_ref(:,k),plotopt);

    axes(ha(k)+length(ha)/2);
    h2(k)=plot(K_ref,yi_ref(:,k),plotopt);
    
    % Move down, just above patches (uncertainties)
    n_move=stackdownpatch(ha(k));
    for j=1:n_move; uistack(h1(k),'down'); end

    n_move=stackdownpatch(ha(k)+length(ha)/2);
    for j=1:n_move; uistack(h2(k),'down'); end

end

end

%%

function n_move=stackdownpatch(hax)

    n_move=[];

    hc=get(hax,'Children');
    for j=1:length(hc)
        if strcmpi(hc(j).Type,'patch')
            n_move=j-2;
            break
        end
    end

end

