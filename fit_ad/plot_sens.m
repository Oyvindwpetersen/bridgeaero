function plot_sens(model,neg_logL_grad)
%% Create options for model: kernel, noise, hyperparameters, bounds, ...
%
% Inputs:
% model: struct with GPR model
% neg_logL_grad: gradient of -log(L) wrt. hyperparameters
%
% Outputs:
%

%%

figure(); hold on; grid on;

for idx1=1:size(model,1)
    for idx2=1:size(model,2)

        i1=model{idx1,idx2}.idx.glob(model{idx1,idx2}.idx.d);
        plot(i1,abs(neg_logL_grad(i1)),'ks');

        i2=model{idx1,idx2}.idx.glob(model{idx1,idx2}.idx.sigma_v);
        plot(i2,abs(neg_logL_grad(i2)),'ob');

        i3=model{idx1,idx2}.idx.glob(model{idx1,idx2}.idx.sigma);
        plot(i3,abs(neg_logL_grad(i3)),'xr');

        i4=model{idx1,idx2}.idx.glob(model{idx1,idx2}.idx.L);
        plot(i4,abs(neg_logL_grad(i4)),'md');

        for k=1:length(i1)
            yl{i1(k)}=['d' num2str(k)];
        end
        for k=1:length(i2)
            yl{i2(k)}=['\sigma_v' num2str(k)];
        end
        for k=1:length(i3)
            yl{i3(k)}=['\sigma' num2str(k)];
        end

        for k=1:length(i4)
            yl{i4(k)}=['L' num2str(k)];
        end
        
    end
end

ylog;

xlabel('Parameter');
ylabel('Sensitivity');

xticks([1:length(yl)]);
xticklabels(yl);