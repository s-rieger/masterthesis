%% RUN FORWARD

clear all;
close all;
clc;
tic

%% Input Parameters
v1=1500;   
v2=1420;
v3=1420;
v4=1440;
v5=1525;
v6=1560;

m1=8;    
m2=5;
m3=15;
m4=44;
m5=28;   

v=[v1,v2,v3,v4,v5,v6];   
m=[m1,m2,m3,m4,m5];
SVP_mod=[m,v];

z_SF=sum(m);              %For creation of Seafloor

lin_stp=1;                %SVP type:1=linear, 2=stepped 
dz=sum(m)*0.005;           %Incremention of depth for SVP/Raytrace  1% of depth

dlmwrite('SVP_mod',SVP_mod);

dtheta=1;
theta=(0:dtheta:55);           %Opening angle and inclimation for Raytrace (more then 45 deg requires larger SF dimension)

no_type=5;                   %Seafloor type:1=flat, 2=tilted, 3=wavy, 4=random, 5=rastrigin
no_vt=3;                     %Track type: 1=angled, 2=parallel, 3=sinously

q=10;                     %Number of Pings on each track
dis_track=1.5*z_SF;            %Distance of the vessel tracks [m]
dis_ping=0.2;             %Distance bewtween each Ping [m]
head=10;                 %Angle of P1 to P2 [deg]
start_off=0.5;            %Offset of P1 for more intersections [m]

mk_sq=dis_track/(dis_ping*q); %Ratio for squaring grids to better find intersection with SF
dlmwrite('mk_sq',mk_sq);



for iii=1:no_type
  for jjj=1:no_vt
  
vt=jjj;
type=iii;

fprintf('Seafloor:%d Vesseltrack:%d \n',type,vt);
%% SVP Model
[v,z]=fun_SVP_100m(SVP_mod,dz,lin_stp);  %[v(z),z]

clear fun_SVP
%% Vessel Tracks
[p1,p2,alpha1,alpha2,phi1,phi2]=fun_vesseltrack_100m(q,dis_track,dis_ping,head,start_off,vt); 

dlmwrite(sprintf('p1_vt_%d',vt),p1);
dlmwrite(sprintf('p2_vt_%d',vt),p2);

p1_cell=num2cell(p1,1);
p2_cell=num2cell(p2,1);

clear fun_vesseltrack
%% Creation of Seafloor Model
x=linspace(0,(start_off+q*dis_ping),10);       
y=linspace(0,max(p2(2,:)),10);
[X,Y]=meshgrid(x,y);

[Xq,Yq,Zq,SF]=fun_seafloor_100m(type,X,Y,x,y,z_SF);  %writes SF_mod_X/Y/Z

SF_sq=[SF(:,1)*mk_sq, SF(:,2),SF(:,3)];  %Interpolate the Grid to square for better intersection

% figure(iii)  % Seafloor
% %title('Seafloor')
% surf(Xq,Yq,Zq);
% view(-60,60);
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('depth [m]');


clear fun_seafloor
clear SF;
%% Raytrace
z_neu=-z;
z_neu(end)=[];

[T_ray,X_ray,p1_ray_xyz,p2_ray_xyz,p1_ray_xyz_sq,p2_ray_xyz_sq]=fun_raytrace_100m(q,theta,phi1,phi2,p1_cell,p2_cell,mk_sq,v,z,z_neu,dz);

clear fun_raytrace
%% Calculating intersection of Seafloor with Rays and find Traveltimes
P1_Data=[];
P2_Data=[];
T_Boden_P1=cell(q-2,length(theta));
T_Boden_P2=cell(q-2,length(theta));

[SF_P1,SF_P2]=fun_intersec_SF_100m(q,theta,z_neu,p1_ray_xyz_sq,p2_ray_xyz_sq,SF_sq);  %SF_P1/P2=depth to seafloor 

for i=2:q-1       
    for k=1:length(theta)          
        T_Boden_P1{i-1,k}=[p1_cell{i}(1,1),p1_cell{i}(2,1),phi1(i-1),theta(k),T_ray{k}(SF_P1(i-1,k),1)]';  %[x,y,phi,theta,T_SF]
        T_Boden_P2{i-1,k}=[p2_cell{i}(1,1),p2_cell{i}(2,1),phi2(i-1),theta(k),T_ray{k}(SF_P2(i-1,k),1)]';  
        
        P1_T_Boden=cell2mat(T_Boden_P1(i-1,:))';
        P2_T_Boden=cell2mat(T_Boden_P2(i-1,:))';               
    end

    P1_Data=cat(1,P1_Data,P1_T_Boden);
    P2_Data=cat(1,P2_Data,P2_T_Boden);              
end

clear fun_intersec
clear T_Boden_P1
clear T_Boden_P2
clear P1_T_Boden
clear P2_T_Boden
%% Save Data for PSO-Inversion
Info=[(length(theta)),theta,(length(p1)-2)];               
dlmwrite('theta_info',Info);                                %[No. of theta,theta_1,....,theta_n,No. of pings]
dlmwrite(sprintf('P1_Data_SF_%d_VT_%d',type,vt),P1_Data);   %[x,y,phi,theta,T_SF]
dlmwrite(sprintf('P2_Data_SF_%d_VT_%d',type,vt),P2_Data);

%% Projection of Seafloor
for i=1:q-2
   for k=1:length(theta)
             
xx1(i,k)=[p1_ray_xyz{i,k}(SF_P1(i,k),1)];
xx2(i,k)=[p2_ray_xyz{i,k}(SF_P2(i,k),1)];

yy1(i,k)=[p1_ray_xyz{i,k}(SF_P1(i,k),2)];
yy2(i,k)=[p2_ray_xyz{i,k}(SF_P2(i,k),2)];

zz1(i,k)=[p1_ray_xyz{i,k}(SF_P1(i,k),3)];
zz2(i,k)=[p2_ray_xyz{i,k}(SF_P2(i,k),3)];

   end 
end

xx=cat(1,xx1,xx2);
xx=xx(:)*mk_sq;     %mk_squar for interpolation of scattered data to square grid
yy=cat(1,yy1,yy2);
yy=yy(:);
zz=cat(1,zz1,zz2);
zz=zz(:);

xx_mesh=linspace(min(Xq(1,:)*mk_sq),max(Xq(1,:)*mk_sq),100);  %ursprünglich min(xx)/min(yy) aber wegen profile auf Xq um gleiche x_position im vergleich zu bekommen
yy_mesh=linspace(min(Yq(:,1)),max(Yq(:,1)),100);

[xxq,yyq]=meshgrid(xx_mesh,yy_mesh);
zzq=griddata(xx,yy,zz,xxq,yyq,'linear');

xxq=xxq/mk_sq;

dlmwrite(sprintf('SF_%d_VT_%d_X_proj',type,vt),xxq);
dlmwrite(sprintf('SF_%d_VT_%d_Y_proj',type,vt),yyq);
dlmwrite(sprintf('SF_%d_VT_%d_Z_proj',type,vt),zzq);


%% Profile from Forward

x_prof=xxq(:,length(xxq)/2);
y_prof=yyq(:,length(yyq)/2);
z_prof=zzq(:,length(zzq)/2);

%% Profile from SF-Model

x_prof_SF=Xq(:,length(xxq)/2);
y_prof_SF=Yq(:,length(yyq)/2);
z_prof_SF=Zq(:,length(zzq)/2);
 
 

end
end



%% Pli Pla Plottz

% figure(1)  %SVP
% plot(v,z)
% grid on
% ylabel('depth [m]')
% xlabel('velocity [m/s]')
% axis ij
% axis ([1400 1600 0 140])

% figure(2)  %Vessel Tracks
% hold on
% %title('Vessel Trails')
% plot(p1(1,:),p1(2,:),'go')
% plot(p2(1,:),p2(2,:),'mo')
% grid on
% xlabel('X [m]')
% ylabel('Y [m]')
% axis([0 2 -10 110]);
%drawArrow(pxa2,pya2);%'color','m');
%drawArrow(pxa1,pya1);%,'color','g')

% figure(3)  % Seafloor
% %title('Seafloor')
% surf(Xq,Yq,Zq);
% view(-60,60);
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('depth [m]');

 figure(4) % Overview figure of the Raytrace
 hold on
 %title('Forward Model')
 for i = 1:q-2
     for k=1:length(theta)
         plot3(p1_ray_xyz{i,k}(1:SF_P1(i,k),1),p1_ray_xyz{i,k}(1:SF_P1(i,k),2),p1_ray_xyz{i,k}(1:SF_P1(i,k),3),'g');
         plot3(p2_ray_xyz{i,k}(1:SF_P2(i,k),1),p2_ray_xyz{i,k}(1:SF_P2(i,k),2),p2_ray_xyz{i,k}(1:SF_P2(i,k),3),'m'); 
     end
 end
 plot3(p1(1,:)',p1(2,:)',p1(3,:)','go',p2(1,:)',p2(2,:)',p2(3,:)','mo');
 surf(Xq,Yq,Zq);
 %alpha 0.5
 xlabel('X [m]')
 ylabel('Y [m]')
 zlabel('depth [m]')
 axis ij
 grid on;
 view(-60,25);
% 
% figure(5) % Seafloor projection
% surf(xxq,yyq,zzq)
% hold on
% for i=1:q-2
%     for k=1:length(theta)
%         plot3(p1_ray_xyz{i,k}(SF_P1(i,k),1),p1_ray_xyz{i,k}(SF_P1(i,k),2),p1_ray_xyz{i,k}(SF_P1(i,k),3),'g+')
%         plot3(p2_ray_xyz{i,k}(SF_P2(i,k),1),p2_ray_xyz{i,k}(SF_P2(i,k),2),p2_ray_xyz{i,k}(SF_P2(i,k),3),'m+')
% 
%     end
% end
% view(-60,5);
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('depth [m]');
% set(gca,'Ydir','reverse')
% %axis([0 2 0 120 -90.5 -89.5]);
% 
% figure(6) % Profiles
% hold on
% surf(xxq,yyq,zzq)
% colormap(gray)
% hold on
% for i=1:q-2
%     for k=1:length(theta)
%         plot3(p1_ray_xyz{i,k}(SF_P1(i,k),1),p1_ray_xyz{i,k}(SF_P1(i,k),2),p1_ray_xyz{i,k}(SF_P1(i,k),3),'g+')
%         plot3(p2_ray_xyz{i,k}(SF_P2(i,k),1),p2_ray_xyz{i,k}(SF_P2(i,k),2),p2_ray_xyz{i,k}(SF_P2(i,k),3),'m+')
% 
%     end
% end
% plot3(x_prof,y_prof,z_prof,'ro','Markersize',10,'MarkerFaceColor','r')
% %plot3(xxq,yyq,zzq,'k')  % nur die in y-richtung liegenden profile
% view(135,15);
% grid on
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('depth [m]')

% figure(7)
% hold on
% plot(y_prof,z_prof,'r')
% plot(y_prof_SF,z_prof_SF,'g')
% grid on
% xlabel('Y [m]')
% ylabel('depth [m]')
% legend('forward projection','modelled SF')
toc
  