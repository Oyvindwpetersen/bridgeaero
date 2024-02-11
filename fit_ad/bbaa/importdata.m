
clc
clear all
% close all

%%

FolderAD=['C:\Cloud\OD_OWP\Work\Projects\Long_span_bridges\Langenuen\WindTunnelData\Results\'];

% H_axis=[5.5 5.8 6.1 6.4 6.7 7.0]; Year=2020;
H_axis=[4.9 5.2 5.5 5.8 6.1]; Year=2021;

for k=1:length(H_axis)
    
    SectionNameCell{k}=['LN' num2str(Year-2000) '_' num2str(H_axis(k)*1000)];
    
    Section_data=load([FolderAD '\' SectionNameCell{k} '\' 'AD_' SectionNameCell{k} '.mat']);
    
    Section_AD{1,k}=Section_data.AD;
    Section_K{1,k}=1./Section_data.Vred;
    Section_H{1,k}=H_axis(k)*ones(size(Section_K{1,k}));

end

Section_AD=cell2mat(Section_AD).';
Section_K=cell2mat(Section_K).';
Section_H=cell2mat(Section_H).';

%%

AD_sort=...
{ 'P1' 'P2' 'P3' 'P4' 'P5' 'P6' ...
    'H1' 'H2' 'H3' 'H4' 'H5' 'H6' ...
    'A1' 'A2' 'A3' 'A4' 'A5' 'A6' };

AD_mat_d{1,1}='P1';
AD_mat_d{1,2}='P5';
AD_mat_d{1,3}='P2';
AD_mat_d{2,1}='H5';
AD_mat_d{2,2}='H1';
AD_mat_d{2,3}='H2';
AD_mat_d{3,1}='A5';
AD_mat_d{3,2}='A1';
AD_mat_d{3,3}='A2';


AD_mat_s{1,1}='P4';
AD_mat_s{1,2}='P6';
AD_mat_s{1,3}='P3';
AD_mat_s{2,1}='H6';
AD_mat_s{2,2}='H4';
AD_mat_s{2,3}='H3';
AD_mat_s{3,1}='A6';
AD_mat_s{3,2}='A4';
AD_mat_s{3,3}='A3';


data_AD_s=[];
data_AD_d=[];
for k1=1:3
for k2=1:3

    idx_d=cellsubindex(AD_mat_d{k1,k2},AD_sort);
    idx_s=cellsubindex(AD_mat_s{k1,k2},AD_sort);    

    data_AD_s(k1,k2,:)=Section_AD(:,idx_d);
    data_AD_d(k1,k2,:)=Section_AD(:,idx_d);

    data_K(k1,k2,:)=Section_K(:,idx_d);
    data_H(k1,k2,:)=Section_H(:,idx_d);

end
end

% col_h1=7;
% col_h4=10;
% 
% H1=Section_AD(col_h1,:).';
% H4=Section_AD(col_h4,:).';
% 
% H=Section_H(col_h1,:).';
% K=Section_K(col_h1,:).';



