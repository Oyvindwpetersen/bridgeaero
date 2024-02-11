function ha=plotadsplit(test_matrix,yt_r,yt_i,model,K_pred,varargin)

%% Plot separate
%
% Inputs:
% test_matrix: [M,2] matrix with K and x as columns
% yt_r: [M,1] vector with real part of transfer function
% yt_i: [M,1] vector with imaginary part of transfer function
% model: struct with model valuesparameter options
% K_pred: vector with frequencies to predict
%
% Outputs:
% 
%

%%

p=inputParser;
addParameter(p,'uncertainty',true,@islogical)
addParameter(p,'gap',[0.2 0.125],@isnumeric)
addParameter(p,'marg_h',[0.15 0.15],@isnumeric)
addParameter(p,'marg_w',[0.05 0.05],@isnumeric)
addParameter(p,'marker','o',@ichar)
addParameter(p,'markersize',3,@isnumeric)
addParameter(p,'color',[0 0 0 ; 0 0 1],@isnumeric)
addParameter(p,'n_sd',[2],@isnumeric)
addParameter(p,'xlabel','K',@ischar)
addParameter(p,'ylabel',{'K^2 AD_{stiffness}' 'K^2 AD_{damping}'},@iscell)
addParameter(p,'title',{''},@iscell)

parse(p,varargin{:})

uncertainty=p.Results.uncertainty;
gap=p.Results.gap;
marg_h=p.Results.marg_h;
marg_w=p.Results.marg_w;
marker=p.Results.marker;
markersize=p.Results.markersize;
color=p.Results.color;
n_sd=p.Results.n_sd;
xlab=p.Results.xlabel;
ylab=p.Results.ylabel;
tit=p.Results.title;

%%

xt_uni=uniquetol(test_matrix(:,2));

n=length(xt_uni);

figure();
ha=tight_subplot(2,n,gap,marg_h,marg_w);
get(gca,'YLabel');

plotopt_train=struct();
plotopt_train.linestyle='None';
plotopt_train.marker=marker;
plotopt_train.markersize=markersize;
plotopt_train.color=color(1,:);
plotopt_train.displayname='Data';

plotopt_pred=struct();
plotopt_pred.color=color(2,:);
plotopt_pred.displayname='Prediction';

for k=1:n

    idx_h=find(ismembertol(test_matrix(:,2),xt_uni(k)));

    pred_matrix=gridvec(K_pred,xt_uni(k));

    [y_pred,yr_pred,yi_pred,std_yr_pred,std_yi_pred,std_yr_obs,std_yi_obs,a_pred,cov_a_pred]=ad_gpr_pred(test_matrix,pred_matrix,[yt_r;yt_i],model);

    % Stiffness AD
    axes(ha(k)); hold on; grid on; 
    plot(K_pred,yr_pred,plotopt_pred);
    plot(test_matrix(idx_h,1),yt_r(idx_h),plotopt_train);

    if uncertainty
        displayname1=['\pm' num2str(n_sd) '\sigma '];
        displayname2=['\pm' num2str(n_sd) '\sigma '];

        h_shade=plotci(K_pred,yr_pred,std_yr_pred,n_sd,'Color',[0 0 1],'displayname',displayname1);
        h_shade=plotci(K_pred,yr_pred,std_yr_obs,n_sd,'Color',[0.25 0.25 0.25],'displayname',displayname2);
    end
    
    xlabel(xlab);
    if k==1; ylabel(ylab{1}); end
    axistight(gca,[0 0.05],'x','y');

    if ~isempty(tit)
        title(tit{k},'FontSize',8,'FontWeight','normal');
    end
    
    % Damping AD
    axes(ha(k+n)); hold on; grid on;
    plot(K_pred,yi_pred,plotopt_pred);
    plot(test_matrix(idx_h,1),yt_i(idx_h),plotopt_train);

    if uncertainty
        h_shade=plotci(K_pred,yi_pred,std_yi_pred,n_sd,'Color',[0 0 1]);
        h_shade=plotci(K_pred,yi_pred,std_yi_obs,n_sd,'Color',[0.25 0.25 0.25]);
    end
    
    xlabel(xlab);
    if k==1; ylabel(ylab{2}); end
    axistight(gca,[0 0.05],'x','y');

end

axes(ha(1)); legend show;
