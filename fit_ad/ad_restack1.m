function AD=ad_restack1(AD_stiff,AD_damp)
%% Convert ADs: from 3*3 (stiffness and damping) to single vector
%
% Inputs
% AD_stiff: 3*3*N matrix
% AD_damp: 3*3*N matrix
%
% Outputs:
% AD: 18*N matrix
%

%%
%
% AD_label_vec=...
%     { 'P1' 'P2' 'P3' 'P4' 'P5' 'P6' ...
%     'H1' 'H2' 'H3' 'H4' 'H5' 'H6' ...
%     'A1' 'A2' 'A3' 'A4' 'A5' 'A6' };
% 
% AD_d{1,1}='P1'; AD_d{1,2}='P5'; AD_d{1,3}='P2';
% AD_d{2,1}='H5'; AD_d{2,2}='H1'; AD_d{2,3}='H2';
% AD_d{3,1}='A5'; AD_d{3,2}='A1'; AD_d{3,3}='A2';
% 
% AD_s{1,1}='P4'; AD_s{1,2}='P6'; AD_s{1,3}='P3';
% AD_s{2,1}='H6'; AD_s{2,2}='H4'; AD_s{2,3}='H3';
% AD_s{3,1}='A6'; AD_s{3,2}='A4'; AD_s{3,3}='A3';

idx_s=[...
    4   6   3 ;
    12  10   9 ;
    18  16  15];

idx_d=[...
    1    5    2
    11    7    8
    17   13   14];

AD=zeros(18,size(AD_stiff,3));
for k1=1:3
    for k2=1:3

        AD(idx_s(k1,k2),:)=squeeze(AD_stiff(k1,k2,:));
        AD(idx_d(k1,k2),:)=squeeze(AD_damp(k1,k2,:));

    end
end
