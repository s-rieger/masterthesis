function [SF_P1,SF_P2]=fun_intersec_SF_100m(q,theta,z_neu,p1_ray_xyz_sq,p2_ray_xyz_sq,SF_sq,mk_sq)
diff_P1=cell(q-2,length(theta),length(z_neu));
diff_P2=cell(q-2,length(theta),length(z_neu));
C_P1=nan(1,length(z_neu));
C_P2=nan(1,length(z_neu));
SF_P1=nan(q-2,length(theta));
SF_P2=nan(q-2,length(theta));

for i=2:q-1        
    for k=1:length(theta)
        for l=1:length(z_neu)  
            
            %diff_P1{i,k,l}=abs([p1_ray_xyz{i,k}(l,1)*mk_sq,p1_ray_xyz{i,k}(l,2),p1_ray_xyz{i,k}(l,3)]-[SF(:,1)*mk_sq,SF(:,2),SF(:,3)]);
            %diff_P2{i,k,l}=abs([p2_ray_xyz{i,k}(l,1)*mk_sq,p2_ray_xyz{i,k}(l,2),p2_ray_xyz{i,k}(l,3)]-[SF(:,1)*mk_sq,SF(:,2),SF(:,3)]);
            
            diff_P1{i-1,k,l}=abs(p1_ray_xyz_sq{i-1,k}(l,:)-SF_sq);
            diff_P2{i-1,k,l}=abs(p2_ray_xyz_sq{i-1,k}(l,:)-SF_sq);           
                                   
            sum1=sum(diff_P1{i-1,k,l}(:,:),2);
            sum2=sum(diff_P2{i-1,k,l}(:,:),2);
            
            [val_P1,idx]=min(sum1);
            [val_P2,idx]=min(sum2);
            
            C_P1(l)=val_P1;
            C_P2(l)=val_P2; 
            
            clear diff_P1
            clear diff_P2
        end
            [vl_1,ix_P1]=min(C_P1);
            [vl_2,ix_P2]=min(C_P2);
        
            SF_P1(i-1,k)=[ix_P1];
            SF_P2(i-1,k)=[ix_P2];
            
%             clear C_P1
%             clear C_P2
%             clear ix_P1
%             clear ix_P2          
% 
%             clear sum1
%             clear sum2
%             clear val_P1
%             clear val_P2
%             clear idx

    end

end