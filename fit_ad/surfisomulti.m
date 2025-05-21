function ha=surfisomulti(plot_matrix,f,test_matrix,f_point,varargin)

%% Plot mulitple 2D surface with isolines
%
% Inputs:
% plot_matrix: [n1,n2] cell with input data matrices (columns)
% f: [n1,n2] cell with surface function data (columns)
% test_matrix: [n1,n2] cell with test data matrices (columns)
% f_point: [n1,n2] cell with discrete test data (columns)
%
% Outputs:
%

%%

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'gap',[0.2 0.075])
addParameter(p,'marg_h',[0.15 0.15])
addParameter(p,'marg_w',[0.05 0.05])

addParameter(p,'xlabel','x')
addParameter(p,'ylabel','y')
addParameter(p,'zlabel','z')
% addParameter(p,'facealpha',0.8,@isnumeric)
% addParameter(p,'displayname','',@ischar)
% addParameter(p,'isolines',[20],@isnumeric) % Number of isolines between [f_min,f_max]
% addParameter(p,'linestyle','-',@ischar)
% addParameter(p,'linewidth',0.3,@isnumeric)
addParameter(p,'view',[40 15],@isnumeric)
% addParameter(p,'cbar',[0.5 1 0 0],@isnumeric) % Set to NaN to turn off colorbar
addParameter(p,'xtick',[],@isnumeric)
addParameter(p,'ytick',[],@isnumeric)

parse(p,varargin{1:end});

gap=p.Results.gap;
marg_h=p.Results.marg_h;
marg_w=p.Results.marg_w;

xlab=p.Results.xlabel;
ylab=p.Results.ylabel;
zlab=p.Results.zlabel;
% facealpha=p.Results.facealpha;
% displayname=p.Results.displayname;
% isolines=p.Results.isolines;
% linestyle=p.Results.linestyle;
% linewidth=p.Results.linewidth;
viewvec=p.Results.view;
% cbar=p.Results.cbar;
xtick=p.Results.xtick;
ytick=p.Results.ytick;

%%

plotopt=struct();
plotopt.xlabel=xlab;
plotopt.ylabel=ylab;
plotopt.view=viewvec;
plotopt.linestyle='-';
plotopt.linewidth=0.1;
plotopt.cbar=[0.5 1 0 0.1];
plotopt.cbar=NaN;
plotopt.xtick=xtick; %[0 1];
plotopt.ytick=ytick; %[2:2:8]

plotopt2=struct();
plotopt2.marker='o';
plotopt2.markersize=1;
plotopt2.color=[0 0 0];
plotopt2.linestyle='None';

%%

if ~iscell(plot_matrix)
    plot_matrix={plot_matrix};
end

if ~iscell(f)
    f={f};
end

if ~iscell(test_matrix)
    test_matrix={test_matrix};
end

if ~iscell(f_point)
    f_point={f_point};
end

[n1,n2]=size(f);

for idx1=1:n1
    for idx2=1:n2

        if ~isnumeric(plot_matrix{idx1,idx2})
            error('Must be a matrix');
        end

        if ~isnumeric(f{idx1,idx2})
            error('Must be a matrix');
        end

        if ~isnumeric(test_matrix{idx1,idx2})
            error('Must be a matrix');
        end

        if ~isnumeric(f_point{idx1,idx2})
            error('Must be a matrix');
        end

    end
end

%%

figure();
ha=tight_subplot(n1,n2,gap,marg_h,marg_w);

c=0;
for idx1=1:n1
    for idx2=1:n2

        c=c+1;

        zl=zlab{idx1,idx2}; 

        axes(ha(c)); hold on; grid on;
        surfiso(plot_matrix{idx1,idx2},f{idx1,idx2},plotopt,'zlabel',zl);
        plot3(test_matrix{idx1,idx2}(:,1),test_matrix{idx1,idx2}(:,2),f_point{idx1,idx2},plotopt2);

    end
end


% for k=1:length()
Link = linkprop([ha],{'View','XLim', 'YLim'});
setappdata(gcf, 'StoreTheLink', Link);


