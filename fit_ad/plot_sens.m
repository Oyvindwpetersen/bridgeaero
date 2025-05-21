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
addParameter(p,'markersize',5,@isnumeric)
addParameter(p,'normalize',true,@islogical)

parse(p,varargin{:})

markersize=p.Results.markersize;
normalize=p.Results.normalize;

%%

plotopt=struct();
plotopt.markersize=markersize;

figure(); hold on; grid on;

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
            yl{i1(k)}=['d_' num2str(k)];
        end
        for k=1:length(i2)
            yl{i2(k)}=['\sigma_{v,' num2str(k) '}'];
        end
        for k=1:length(i3)
            yl{i3(k)}=['\sigma_' num2str(k)];
        end

        for k=1:length(i4)
            yl{i4(k)}=['L_' num2str(k)];
        end
        
    end
end

xlabel('Parameter');

if normalize==true; ylabel({'Normalized' 'sensitivity'});
else; ylabel({'Sensitivity'});
end

xticks([1:length(yl)]);
xticklabels(yl);