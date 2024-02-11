%%


clc
clear all
close all


mod=importbeam(100,210e9,0.25/100,1000,100);

x_axis=mod.nodecoord(:,2);

f=mod.f;
phi=mod.phi;
doflabel=mod.doflabel;

%%%%


K=mod.K;
for k=1:2:size(K,1)
    K(k,k)=K(k,k)+1e5;
end

K=mod.S_red.'*K*mod.S_red;

[v,d]=eig(K,mod.M_red);

dd=diag(d);
[~,i_sort]=sort(dd);
oms=dd(i_sort);
phi=v(:,i_sort);

phi=mod.S_red*phi

f=sqrt((oms))/2/pi;

figure();
for k=1:5
    hold on;
    plot(phi(1:2:end,k));
end

tilefigs






% clear mod


%%



dt=1/(1000*10);
% dt=0.01;
% dt=0.05;

modes=[1:10];


[Sd,Sa,Sp]=DofSelection({'50_U'},{},{'10_U' '11_U'},doflabel)


Omega=diag(f(modes)*2*pi);
% Xi=diag([0.01 0.01 0.01 0.1*ones(1,7)])
Xi=diag([0.01 ]*ones(1,length(modes)));
Gamma=2*Omega.*Xi;

[A B G J Ac Bc Gc Jc F]=ssmod_modal(phi(:,modes),Omega,Gamma,Sa,Sd,[],dt);


% V_all=[10:10:100];
V_all=[1 20 40 80];

phi_u=phi(1:2:end,modes);

for j=1:length(V_all)

    T=max(x_axis)/V_all(j);

    t=[0:dt:(T)*1.0];

    x=zeros(size(A,1),length(t));
    y=zeros(1,length(t));

    for k=1:length(t)

        x_current=V_all(j)*t(k);

        [idx,weight]=findpos(x_axis,x_current);

        p_modal=phi_u(idx,:).'*weight';

        x(:,k+1)=A*x(:,k)+B*p_modal;
        y(:,k)=G*x(:,k)+J*p_modal;

    end

    x=x(:,1:length(t));
    
    t_all{j}=t;
    x_all{j}=x;
    y_all{j}=y;

end


%%


% close all

plottime(...
    t_all{1}*V_all(1),y_all{1},...
    t_all{2}*V_all(2),y_all{2},...
    t_all{3}*V_all(3),y_all{3},...
    t_all{4}*V_all(4),y_all{4},...
'linestyle',{'-' '--' '--' '-'}...
    )


% plottime(...
%     t_all{1},y_all{1}...
%     )



plotfreq(...
    t_all{1},y_all{1},...
    t_all{2},y_all{2},...
    t_all{3},y_all{3},...
    t_all{4},y_all{4},...
'linestyle',{'-' '--' '--' '-'},...
    'xlim',[0 500]);


