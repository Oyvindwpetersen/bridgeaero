function plot_sens(model,neg_logL_grad,varargin)

%% Create options for model: kernel, noise, hyperparameters, bounds, ...
%
% Inputs:
% model: struct with GPR model
% neg_logL_grad: gradient of -log(L) wrt. hyperparameters
%
% Outputs:
%

%%

p=inputParser;
p.KeepUnmatched=true;
addParameter(p,'markersize',3,@isnumeric)
addParameter(p,'normalize',true,@islogical)
addParameter(p,'vertline',true,@islogical)
addParameter(p,'titles',{},@iscell)

addParameter(p,'gap',[0.2 0.125],@isnumeric)
addParameter(p,'marg_h',[0.15 0.15],@isnumeric)
addParameter(p,'marg_w',[0.05 0.05],@isnumeric)

parse(p,varargin{:})

markersize=p.Results.markersize;
normalize=p.Results.normalize;
vertline=p.Results.vertline;
titles=p.Results.titles;

gap=p.Results.gap;
marg_h=p.Results.marg_h;
marg_w=p.Results.marg_w;

%%

plotopt=struct();
plotopt.markersize=markersize;

figure(); 
ha=tight_subplot(1,1,gap,marg_h,marg_w);

hold on; grid on;
ylog;

% i1=model{1,1}.idx.glob(model{1,1}.idx.d);
% plot(i1,abs(neg_logL_grad(i1)),'ks',plotopt);
% sens=abs(neg_logL_grad(i2));
% if normalize==true; sens=sens./model{idx1,idx2}.hyp.sigma_v; end
% plot(i2,sens,'ob',plotopt);

i1=model{1,1}.idx.glob(model{1,1}.idx.d);
sens=abs(neg_logL_grad(i1));
if normalize==true; sens=sens./model{1,1}.hyp.d; end
plot(i1,sens,'sk',plotopt);

vertline_x=0;
vertline_x(end+1)=i1(end);

dollar='$';

for idx1=1:size(model,1)
    for idx2=1:size(model,2)

        i2=model{idx1,idx2}.idx.glob(model{idx1,idx2}.idx.sigma_v);
        sens=abs(neg_logL_grad(i2));
        if normalize==true; sens=sens./model{idx1,idx2}.hyp.sigma_v; end
        plot(i2,sens,'ob',plotopt);

        i3=model{idx1,idx2}.idx.glob(model{idx1,idx2}.idx.sigma);
        sens=abs(neg_logL_grad(i3));
        if normalize==true; sens=sens./model{idx1,idx2}.hyp.sigma; end
        plot(i3,sens,'xr',plotopt);

        i4=model{idx1,idx2}.idx.glob(model{idx1,idx2}.idx.L);
        sens=abs(neg_logL_grad(i4));
        if normalize==true; sens=sens./model{idx1,idx2}.hyp.L; end
        plot(i4,sens,'dm',plotopt);

        for k=1:length(i1)
            xl{i1(k)}=[dollar 'd_' num2str(k) dollar];
        end
        for k=1:length(i2)
            xl{i2(k)}=[dollar '\sigma_{v,' num2str(k) '}' dollar];
        end
        for k=1:length(i3)
            if k==1 | k==2
                sub_k=k;
                a_str='a';
            else
                sub_k=k+1;
                a_str='a';
            end
            
            xl{i3(k)}=[dollar '\sigma_{ ' a_str ',' num2str(sub_k) '}' dollar];
        end

        for k=1:length(i4)
            if k==1 | k==2
                sub_k=k;
                a_str='a';
            else
                sub_k=k+1;
                a_str='a';
            end
            
            xl{i4(k)}=[dollar 'L_{' a_str ',' num2str(sub_k) '}' dollar];
        end

        vertline_x(end+1)=i4(end);

    end
end

xlabel('Parameter');

if normalize==true; ylabel({'Normalized' 'sensitivity'});
else; ylabel({'Sensitivity'});
end

xticks([1:length(xl)]);
xticklabels(xl);
set(gca,'TickLabelInterpreter','latex')

axistight(gca,[0.05 0.05],'x','ylog');
xlim([0 i4(end)+1]);

if vertline
    h_line=linevertical(gca,vertline_x(2:end-1)+0.5,'--',0.5,[0 0 0],'bottom',true);
end

yl=ylim(); yl=yl(2); y_text=10.^(1.2*log10(yl));

if vertline & ~isempty(titles)

    x_titles=(vertline_x(2:end)+vertline_x(1:end-1))/2+0.5;

    for k=1:length(titles)

        interpreter='latex';

        text(x_titles(k),y_text,titles{k},'HorizontalAlignment','center','interpreter',interpreter,'BackgroundColor',[1 1 1],'Margin',1,...
            'Fontsize',8);
    end

end

setlogtick(gca,5,'y');

