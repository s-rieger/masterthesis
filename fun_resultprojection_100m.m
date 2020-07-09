function [x_PSO_prof,y_PSO_prof,z_PSO_prof]=fun_resultprojection_100m(Xq,Yq,P1_data,P2_data,result_1,dz,lin_stp,theta,q,phi1,phi2,p1,p2)
%% Create Cell for each Ping
kk=1:length(theta):length(P1_data); 

for i=1:length(kk)  
    P1_data_cell{i}=num2cell(P1_data([kk(i):(kk(i)+length(theta)-1)],:),5);
    P2_data_cell{i}=num2cell(P2_data([kk(i):(kk(i)+length(theta)-1)],:),5);
end

%% Load Coordinates of each Ping
% p1=load('p1_vtrack');
% p2=load('p2_vtrack');
p1_cell=num2cell(p1,1);
p2_cell=num2cell(p2,1);

%% Load ratio for Interpolation before finding intersections
load('mk_sq');

%% Define particle atributes
% l=length((:,1))/2-1; %Number of Layers
% m=length((:,1))/2;     %Number of velocities
l=(length(result_1(:,1))-1)/2; %Number of Layers
m=(length(result_1(:,1))-1)/2+1;   %Number of velocities

%% Main Loop for results
for w=1:length(result_1(1,:))

%SVP    
[v,z]=fun_SVP_100m(result_1(:,w),dz,lin_stp);

z_neu=z;  
z_neu(end)=[];

%Raytrace
[T_ray,X_ray,p1_ray_xyz,p2_ray_xyz,p1_ray_xyz_sq,p2_ray_xyz_sq]=fun_raytrace_100m((q+1),theta,phi1,phi2,p1_cell,p2_cell,mk_sq,v,z,z_neu,dz);
%T_ray=num2cell(t_ray,2);

%Find depth of Seafloor
dt_P1_t_Boden=cell(1,length(theta));
dt_P2_t_Boden=cell(1,length(theta));
SF_depth_P1=nan(1,q);
SF_depth_P2=nan(1,q);


for i=1:q
    
    for k=1:length(theta)
        
              dt_P1_t_Boden{k}=abs(T_ray{k}(:,1)-cell2mat(P1_data_cell{i}(k,5)));
              [val_P1,idx_P1] = min([dt_P1_t_Boden{k}]);
              %diff_P1{k}(i)=[val_P1];   %Anstand Laufzeit_Modelliert und Raytrace
              SF_depth_P1(i,k)=idx_P1;
                       
              dt_P2_t_Boden{k}=abs(T_ray{k}(:,1)-cell2mat(P2_data_cell{i}(k,5)));
              [val_P2,idx_P2] = min([dt_P2_t_Boden{k}]);
              %diff_P2{k}(i)=[val_P2];   %Anstand Laufzeit_Modelliert und Raytrace
              SF_depth_P2(i,k)=idx_P2;
                             
    end
end

%Projection of Seafloor
for i=1:q-2
   for k=1:length(theta)
             
        xx1(i,k)=[p1_ray_xyz{i,k}(SF_depth_P1(i,k),1)];
        xx2(i,k)=[p2_ray_xyz{i,k}(SF_depth_P2(i,k),1)];

        yy1(i,k)=[p1_ray_xyz{i,k}(SF_depth_P1(i,k),2)];
        yy2(i,k)=[p2_ray_xyz{i,k}(SF_depth_P2(i,k),2)];

        zz1(i,k)=[p1_ray_xyz{i,k}(SF_depth_P1(i,k),3)];
        zz2(i,k)=[p2_ray_xyz{i,k}(SF_depth_P2(i,k),3)];

   end 
end


xx=cat(1,xx1,xx2);
xx=xx(:)*mk_sq;     
yy=cat(1,yy1,yy2);
yy=yy(:);
zz=cat(1,zz1,zz2);
zz=zz(:);

xx_mesh=linspace(min(Xq(1,:)*mk_sq),max(Xq(1,:)*mk_sq),100);  %ursprünglich min(xx)/min(yy) aber wegen profile auf Xq um gleiche x_position im vergleich zu bekommen
yy_mesh=linspace(min(Yq(:,1)),max(Yq(:,1)),100);

[xxq,yyq]=meshgrid(xx_mesh,yy_mesh);
zzq=griddata(xx,yy,zz,xxq,yyq,'linear');

xxq=xxq/mk_sq;

%% Profile from PSO

x_PSO_prof{w}=xxq(:,length(xxq)/2);
y_PSO_prof{w}=yyq(:,length(yyq)/2);
z_PSO_prof{w}=zzq(:,length(zzq)/2);

end
