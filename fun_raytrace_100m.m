function [T_ray,X_ray,p1_ray_xyz,p2_ray_xyz,p1_ray_xyz_sq,p2_ray_xyz_sq]=fun_raytrace_100m(q,theta,phi1,phi2,p1_cell,p2_cell,mk_sq,v,z,z_neu,dz)

p=sind(theta)./v(1);   % Strahlparameter definieren
t_ray=nan(length(theta),length(z_neu));
x_ray=nan(length(theta),length(z_neu));

for i=1:length(theta)     
    for k=1:length(z)-1   
             
        t_ray(i,k)=0;
        for j=1:int32(z(k)/dz)
            t_ray(i,k)=t_ray(i,j)+(dz/(v(k)*sqrt(1-p(i)^2*v(k)^2)));  
        end
      
        x_ray(i,k)=0;
        for j=1:int32(z(k)/dz)
            x_ray(i,k)=x_ray(i,j)+(p(i)*v(k))*dz/sqrt(1-p(i)^2*v(k)^2);  
        end                    
    end
end

T_ray=cell(1,length(theta));
X_ray=cell(1,length(theta));

for i=1:length(theta)
  
    T_ray{i}=[t_ray(i,:)',z_neu']; %[t z]
    X_ray{i}=[x_ray(i,:)',zeros(length(z_neu),1),z_neu']; %[x 0 z] 
    
end

%Apply raytrace to all pings
p1_ray_xyz=cell(q-2,length(theta));
p2_ray_xyz=cell(q-2,length(theta));
p1_ray_xyz_sq=cell(q-2,length(theta));
p2_ray_xyz_sq=cell(q-2,length(theta));


for i=2:q-1
    for k=1:length(theta)                       
        p1_ray_xyz{i-1,k}=[X_ray{k}(:,1).*cosd(phi1(i-1))+p1_cell{i}(1,1),X_ray{k}(:,1)*sind(phi1(i-1))+p1_cell{i}(2,1),z_neu'+p1_cell{i}(3,1)];  % [x,y,z]
        p2_ray_xyz{i-1,k}=[X_ray{k}(:,1).*cosd(phi2(i-1))+p2_cell{i}(1,1),X_ray{k}(:,1)*sind(phi2(i-1))+p2_cell{i}(2,1),z_neu'+p2_cell{i}(3,1)];
        
        p1_ray_xyz_sq{i-1,k}=[p1_ray_xyz{i-1,k}(:,1)*mk_sq , p1_ray_xyz{i-1,k}(:,2),p1_ray_xyz{i-1,k}(:,3)];  % X-Value times mk_sq for square grid for better location of intersection
        p2_ray_xyz_sq{i-1,k}=[p2_ray_xyz{i-1,k}(:,1)*mk_sq , p2_ray_xyz{i-1,k}(:,2),p2_ray_xyz{i-1,k}(:,3)];        
    end
end
end