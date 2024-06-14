function [ad_label_s,ad_label_d]=ad_label(idx,type)
%% Labels for aerodynamic derivatives
%
% Inputs:
% idx: input slice of which labels to output, e.g. [2 3] for vertical and torsion only
% type: 1 for regular AD, 2 for K^2*AD
%
% Outputs:
% ad_label_s: cell with labels for stiffness ADs
% ad_label_d: cell with labels for damping ADs
%

%%

ad_label_d{1,1}='P_1'; ad_label_d{1,2}='P_5'; ad_label_d{1,3}='P_2';
ad_label_d{2,1}='H_5'; ad_label_d{2,2}='H_1'; ad_label_d{2,3}='H_2';
ad_label_d{3,1}='A_5'; ad_label_d{3,2}='A_1'; ad_label_d{3,3}='A_2';

ad_label_s{1,1}='P_4'; ad_label_s{1,2}='P_6'; ad_label_s{1,3}='P_3';
ad_label_s{2,1}='H_6'; ad_label_s{2,2}='H_4'; ad_label_s{2,3}='H_3';
ad_label_s{3,1}='A_6'; ad_label_s{3,2}='A_4'; ad_label_s{3,3}='A_3';

dollar='$';

if type==1
    k_squared='';
elseif type==2
    k_squared=' K^2';
end

for idx1=1:3
    for idx2=1:3

        ad_label_d{idx1,idx2}=[dollar ad_label_d{idx1,idx2} k_squared dollar];
        ad_label_s{idx1,idx2}=[dollar ad_label_s{idx1,idx2} k_squared dollar];
    end
end

if any(idx)>3 | any(idx)<1
    error('idx must contain 1,2 and 3 only');
end

ad_label_d=ad_label_d(idx,idx);
ad_label_s=ad_label_s(idx,idx);
