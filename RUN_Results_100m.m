%% RUN Results
tic
clear all;
close all;
clc;

%% Parameters
lin_stp=3;        % combo
no_type=5;       % SF type
no_vt=3;        % Vesseltrack
dz=0.5;           % for Inversion Seafloor projection of PSO results

count=0;
tab_1_all=[];
tab_2_all=[];
%% Main Loop
for iii=1:no_type
  for kkk=1:no_vt
        
type=iii;
vt=kkk;
p1=load(sprintf('p1_vt_%d',vt));
p2=load(sprintf('p2_vt_%d',vt));

% Load SVP-Particles from PSO
results=load(sprintf('results_SF_%d_VT_%d',type,vt));

% Find best misfit
[val idn]=min(results(end-2,:)); 
best_result=results(1:end-3,idn);
no_reset=max(results(end,:));
no_reset_gb=results(end,idn);

l=(length(best_result(:,1))-1)/2; %Number of Layers
m=(length(best_result(:,1))-1)/2+1;   %Number of velocities

% Create SVP from Particles (Combo mode)
[v_neu,z_neu]=fun_SVP_100m(best_result(1:end),dz,lin_stp);

v{kkk}=v_neu;
z{kkk}=z_neu;
%% Load SF Projection from Forward
xxq_for_proj_1=load(sprintf('SF_%d_VT_%d_X_proj',type,vt)); %linear
yyq_for_proj_1=load(sprintf('SF_%d_VT_%d_Y_proj',type,vt));
zzq_for_proj_1=load(sprintf('SF_%d_VT_%d_Z_proj',type,vt));

xxq_for_proj_prof_1=xxq_for_proj_1(:,length(xxq_for_proj_1)/2); 
yyq_for_proj_prof_1=yyq_for_proj_1(:,length(yyq_for_proj_1)/2);
zzq_for_proj_prof_1=zzq_for_proj_1(:,length(zzq_for_proj_1)/2);

top_bound=-(zzq_for_proj_prof_1+(zzq_for_proj_prof_1*0.005));
bottom_bound=-(zzq_for_proj_prof_1-(zzq_for_proj_prof_1*0.005));

%% Load actual SF from Forward
Xq=load(sprintf('SF_%d_mod_X',type));
Yq=load(sprintf('SF_%d_mod_Y',type));
Zq=load(sprintf('SF_%d_mod_Z',type));

x_prof_SF=Xq(:,length(xxq_for_proj_prof_1)/2); 
y_prof_SF=Yq(:,length(xxq_for_proj_prof_1)/2);
z_prof_SF=Zq(:,length(xxq_for_proj_prof_1)/2);

%% Load Traveltimes from Forwrd
P1_data=load(sprintf('P1_Data_SF_%d_VT_%d',type,vt));
P2_data=load(sprintf('P2_Data_SF_%d_VT_%d',type,vt));

%% Load angle information and Number of Pings
load('theta_info');
theta=theta_info(2:end-1);
q=theta_info(end);
phi1=P1_data(1:length(P1_data)/q,3);
phi2=P2_data(1:length(P2_data)/q,3);

%% Load initial SVP of PSO and errorbars
load('SVP_init');

v_init=SVP_init(((length(SVP_init)+1)/2):end);   % Geschwindigkeit
m_init=SVP_init(1:((length(SVP_init)+1)/2)-1);

err_ver=[0,load('SVP_init_M_err')];
err_hor=load('SVP_init_v_err');

 z_init=nan(length(m_init)+1,1);
 z_init(1)=0;
 z_init(2)=m_init(1);

 for k=3:length(m_init)+1
     z_init(k)=z_init(k-1)+m_init(k-1);
 end
%% Load Data from Forward
SVP_mod=load('SVP_mod');

v_mod=SVP_mod(((length(SVP_mod)+1)/2):end);
ma_mod=SVP_mod(1:((length(SVP_mod)+1)/2)-1);

z_mod=nan(length(m)+1,1);
z_mod(1)=0;
z_mod(2)=ma_mod(1);

for k=3:length(ma_mod)+1
    z_mod(k)=z_mod(k-1)+ma_mod(k-1);
end

% Inversion Seafloor projection of PSO results                                                                                                      
[x_PSO_prof_1,y_PSO_prof_1,z_PSO_prof_1]=fun_resultprojection_100m(Xq,Yq,P1_data,P2_data,best_result,dz,lin_stp,theta,q,phi1,phi2,p1,p2);

% mean difference of SVP
[v_modp,z_modp]=fun_SVP_100m([ma_mod';v_mod'],dz,1);

if length(v_modp)<length(v_neu)
    std_dev_svp=sqrt(sum((abs(v_neu(1:length(v_modp))-v_modp)).^2)/length(v_modp)-1)
    mean_diff_svp=sum(v_neu(1:length(v_modp)))/length(v_modp)-sum(v_modp)/length(v_modp)
else
    std_dev_svp=sqrt(sum((abs(v_neu-v_modp(1:length(v_neu)))).^2)/length(v_neu)-1)
    mean_diff_svp=sum(v_neu)/length(v_neu)-sum(v_modp(1:length(v_neu)))/length(v_neu)
end

dlmwrite(sprintf('diff_SVP_%d_VT_%d',type,vt),std_dev_svp);
dlmwrite(sprintf('mean_diff_SVP_%d_VT_%d',type,vt),mean_diff_svp);

% mean difference on profile
diff=sum(abs(z_PSO_prof_1{1}(2:end-1,1)-(-zzq_for_proj_prof_1(2:end-1,1))));
diff=diff/(length(z_PSO_prof_1{1})-2); %procent
diff=(diff/100)*100;
dlmwrite(sprintf('diff_SF_%d_VT_%d',type,vt),diff);

% Plots
no_fig=1+count;

figure(no_fig)
hold on
axis ij
grid on
title(sprintf('Seafloor:%d VT:%d',type,vt));
h1=plot(yyq_for_proj_prof_1,-zzq_for_proj_prof_1,'g','LineWidth',4); %SF Projection forward linear
h2=plot(yyq_for_proj_prof_1,top_bound,'k'); %topbound 0.5%
h3=plot(yyq_for_proj_prof_1,bottom_bound,'k'); %bottom bound 0.5%
h4=plot(y_PSO_prof_1{1}(:,1),z_PSO_prof_1{1}(:,1),'r'); %SF Projections from PSO linear
xlabel('Y [m]')
ylabel('depth [m]')
legend([h1 h2 h4],{'Forward Projection','Projection accuracy','PSO Projection'});
% 
 count=count+1;

 
 tab_1(kkk,:)=[type vt idn no_reset no_reset_gb]; % ];
 tab_2(kkk,:)=[(val)  diff];
  end
 
 tab_1_all=cat(1,tab_1_all,tab_1);
 tab_2_all=cat(1,tab_2_all,tab_2);

 
col=hsv(no_vt);
no_fig=1+count;
figure(no_fig)       %SVP from PSO linear
hold on
grid on
axis ij
title(sprintf('Seafloor:%d',type));
h1=plot(v_mod,z_mod,'Color',[0 0 0]+0.5);                      %SVP forward
h2=plot(v_init,z_init,'k'); %initial svp for pso
errorbar(v_init,z_init,err_ver,'k')
errorbar(v_init,z_init,err_hor,'horizontal','k')

for ppp=1:no_vt
    h3(ppp)=plot(v{ppp},z{ppp},'color',col(ppp,:));
end
ylabel('depth [m]')
xlabel('velocity [m/s]')
legend([h1,h2,h3(1),h3(2),h3(3)],'model','inital','GB VT1','GB VT2','GB VT3');
%axis([1400 1600 0 z_init(end)+10]) 
 count=count+1;
end

toc