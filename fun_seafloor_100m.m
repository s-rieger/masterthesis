function [Xq,Yq,Zq,SF]=fun_seafloor_100m(type,X,Y,x,y,z_SF)
 if type==1 %flat seafloor  
    Z=zeros(length(X),length(Y))-z_SF+15;
    
 elseif type==2 %tilted seafloor   
    Z=ones(length(X),length(Y))-z_SF+10; 
    slope=0;
    for i=1:length(x)
        Z(i,:)=Z(i,:)-slope;
        slope=slope+1;
    end
    
 elseif type==3 %wavey seafloor    
    for i=1:length(x)
        for k=1:length(y)
            Z(k,i)=[10*-sin(40*x(k))-z_SF+30];
        end
    end
    
 elseif type==4 %random seafloor    
    Z=randi([0 10],length(X),length(Y))-z_SF+10;

 elseif type==5 %Rastrigin Funktion    
    Z=nan(length(x),length(y)); 
    for i=1:length(x)
        for j=1:length(y)
            Z(i,j)=((x(i)).^2+((y(j)+25)/30).^2-25*(0.5*cos(5*pi*x(i))+0.5*cos(50*pi*y(j))))-z_SF+20;
        end
    end
 end

%Interpolate the seafloor
xq=linspace(min(x),max(x),100);
yq=linspace(min(y),max(y),100);
[Xq,Yq]= meshgrid(xq,yq);
Zq = interp2(X,Y,Z,Xq,Yq,'spline');

SF=cell(length(Xq),3);

for i=1:length(Xq)
    for j=1:length(Yq)        
        SF{i}(j,:)=[xq(i) yq(j) Zq(j,i)];      
    end    
end

dlmwrite(sprintf('SF_%d_mod_X',type),Xq);
dlmwrite(sprintf('SF_%d_mod_Y',type),Yq);
dlmwrite(sprintf('SF_%d_mod_Z',type),Zq);

SF=cat(1,SF{:});  %[x y z] of Seafloor
end