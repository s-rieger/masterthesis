function [v_neu,z_neu]=fun_SVP_100m(svp,dz,lin_stp)

m=svp(1:((length(svp)+1)/2)-1);
v=svp(((length(svp)+1)/2):end);


 z=nan(length(m)+1,1);
 z(1)=0;

 for k=2:length(m)+1
     z(k)=z(k-1)+m(k-1);
 end

z(end+1)=z(end)+(z(end)/4); % +25% for making sure the raytrace is "long" enough   
v(end+1)=v(end)+(v(end-1)-v(end))/(z(end-1)-z(end))*(z(end)-z(end-1)); %+1 with slope of previous points  
 
z_neu=z(1):dz:z(end);  
        
v_neu=nan(1,length(z_neu));
k=2;

%linear
if lin_stp==1  
    
   
        for i=1:length(z_neu)                 
            if z_neu(i)<=z(k)
              v_neu(i)=v(k-1)+(z_neu(i)-z(k-1))*((v(k)-v(k-1))/(z(k)-z(k-1)));
                                                                                  
            elseif z_neu(i)>z(k)
              k=k+1;
              v_neu(i)=v(k-1)+(z_neu(i)-z(k-1))*((v(k)-v(k-1))/(z(k)-z(k-1)));   
                             
            end
        end  
end

% stepped
if lin_stp==2  
     
        for i=1:length(z_neu)                 
            if z_neu(i)>=z(k)
                k=k+1;
                v_neu(i)=v(k-1);   
                                
            elseif z_neu(i)<z(k)
                v_neu(i)=v(k-1);
                                     
            end
        end  
end

% combo
if lin_stp==3  %last 5 layers linear
    
    for i=1:length(z_neu)      
            if z_neu(i)>=z(end-1)
                v_neu(i)=v(end)+(z_neu(i)-z(end))*((v(end)-v(end-1))/(z(end)-z(end-1)));         
            elseif z_neu(i)>=z(end-2)
                v_neu(i)=v(end-1)+(z_neu(i)-z(end-1))*((v(end-1)-v(end-2))/(z(end-1)-z(end-2)));
            elseif z_neu(i)>=z(end-3)
                v_neu(i)=v(end-2)+(z_neu(i)-z(end-2))*((v(end-2)-v(end-3))/(z(end-2)-z(end-3)));                
            elseif z_neu(i)>=z(end-4)
                v_neu(i)=v(end-3)+(z_neu(i)-z(end-3))*((v(end-3)-v(end-4))/(z(end-3)-z(end-4)));
            elseif z_neu(i)>=z(end-5)
                v_neu(i)=v(end-4)+(z_neu(i)-z(end-4))*((v(end-4)-v(end-5))/(z(end-4)-z(end-5)));
                       
            elseif z_neu(i)>=z(k)
                k=k+1;
                v_neu(i)=v(k-1);                                   
            elseif z_neu(i)<z(k)
                v_neu(i)=v(k-1);                                    
            end
        end  

end  