function [Delta_Z]=fun_delta(particle,P1_data,P2_data,eee,mk_sq,lin_stp,dz,vt);

%% SVP Model
[v,z]=fun_SVP_100m(particle,dz,lin_stp);  %[v(z),z]

z_neu=z;  
z_neu(end)=[];
%% %Winkel und Strahlparameter für Integral

theta=P1_data{1}(:,4); % Öffnungswinkel aus Daten nehmen
p=sind(theta)./v(1);   % Strahlparameter definieren

%% % Numerische Integration Zeit und X

t_ray=nan(length(theta),1);
x_ray=nan(length(theta),1);

for i=1:length(theta)      
    for k=1:length(z)-1  % Letzter überspringen wegen index
       
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

T_ray=num2cell(t_ray,2);

%% Finde im Integral die Tiefe für die Bodenlaufzeit aus Daten

dt_P1_t_Boden=cell(1,length(theta));
dt_P2_t_Boden=cell(1,length(theta));
Tiefe_P1=cell(1,length(P1_data));
Tiefe_P2=cell(1,length(P1_data));
X_P1_data_Boden=cell(1,length(P1_data));
X_P2_data_Boden=cell(1,length(P1_data));

for k=1:length(P1_data)
    
    for i=1:length(theta)
        
              dt_P1_t_Boden{i}=abs(T_ray{i}(1,:)-P1_data{k}(i,5));
              [val_P1,idx_P1] = min([dt_P1_t_Boden{i}]);
              %diff_P1{k}(i)=[val_P1];   %Anstand Laufzeit_Modelliert und Raytrace
              Tiefe_P1{k}(i)=idx_P1;
              X_P1_data_Boden{k}(i)=x_ray(i,Tiefe_P1{k}(i));    % Einsetzten der Tiefe für X-Wert
              
          
              dt_P2_t_Boden{i}=abs(T_ray{i}(1,:)-P2_data{k}(i,5));
              [val_P2,idx_P2] = min([dt_P2_t_Boden{i}]);
              %diff_P2{k}(i)=[val_P2];   %Anstand Laufzeit_Modelliert und Raytrace
              Tiefe_P2{k}(i)=idx_P2;
              X_P2_data_Boden{k}(i)=x_ray(i,Tiefe_P2{k}(i));    % Einsetzten der Tiefe für X-Wert
                             
    end
end

for k=1:length(P1_data)
    Tiefe_P1{k}=Tiefe_P1{k}*dz; % reset back to normal depth in m
    Tiefe_P2{k}=Tiefe_P2{k}*dz;
end
%% Interpolieren der Daten zwischen den Theta

Z_mod_P1=cell(1,length(P1_data));
Z_mod_P2=cell(1,length(P1_data));
X_intp_P1=cell(1,length(P1_data));
X_intp_P2=cell(1,length(P1_data));
Z_intp_P1=cell(1,length(P1_data));
Z_intp_P2=cell(1,length(P1_data));

for i=1:length(P1_data)  % extend interpolation linspace to more than 100!!!
   
   %Profil 1 
   Z_mod_P1{i}=Tiefe_P1{i};  
   X_intp_P1{i}=linspace(min(X_P1_data_Boden{i}),max(X_P1_data_Boden{i}),100);  % Interpoliere linear(!!!) auf 100 segmente für ein Profil 7
   %Z_intp{i}=interp1(X_mod{i},Z_mod{i},X_intp{i},'spline');
   
   %POLYFIT
   po_1=polyfit(X_P1_data_Boden{i},Z_mod_P1{i},5); %1=linear, 2=quadratisch ...
   Z_intp_P1{i}=polyval(po_1,X_intp_P1{i});
   
   %Profil 2
   Z_mod_P2{i}=Tiefe_P2{i};
   X_intp_P2{i}=linspace(min(X_P2_data_Boden{i}),max(X_P2_data_Boden{i}),100);
   
   po_2=polyfit(X_P2_data_Boden{i},Z_mod_P2{i},5);
   Z_intp_P2{i}=polyval(po_2,X_intp_P2{i});
   
   % Dont display the polynomial bad conditioned warining
   [msg, id] = lastwarn;
 
 warning('off', id) 
end

P1_Boden_Inv=cell(1,length(P1_data));
P2_Boden_Inv=cell(1,length(P1_data));


for i=1:length(P1_data)
          
        P1_Boden_Inv{i}=[(X_intp_P1{i}.*cosd(P1_data{i}(1,3))+P1_data{i}(1,1)); ...  % [x,y,z]
                         X_intp_P1{i}.*sind(P1_data{i}(1,3))+P1_data{i}(1,2); ...
                         Z_intp_P1{i}];
                      
        P2_Boden_Inv{i}=[(X_intp_P2{i}.*cosd(P2_data{i}(1,3))+P2_data{i}(1,1)); ...  % [x,y,z]
                         X_intp_P2{i}.*sind(P2_data{i}(1,3))+P2_data{i}(1,2); ...
                         Z_intp_P2{i}];  
        
        % bottom part is normalized through mk_sq to square grid for better
        % intersection calculation
        P1_Boden_Inv_mesh{i}=[P1_Boden_Inv{i}(1,:)*mk_sq;P1_Boden_Inv{i}(2,:);P1_Boden_Inv{i}(3,:)]; % mk_sq ist um quadratisches grid zur besseren intersectionsuche zu erzeugen
        P2_Boden_Inv_mesh{i}=[P2_Boden_Inv{i}(1,:)*mk_sq;P2_Boden_Inv{i}(2,:);P2_Boden_Inv{i}(3,:)]; % geht ja nur zum index                

        
end
%% Finde die Koordinaten für den einen Punkt der am nächsten der Überschneidung kommt

for u=1:length(P2_data)
 for i=1:length(P1_data)  % jeder Messpunkt
  for k=1:length(X_intp_P1{1});% jeder Eintrag (Interpolierte Winkel) % wäre für POLYFIT(XZ_P1_Intep)  
    for m=1:length(X_intp_P1{1})
    
    diff(k,m)=sum(abs((P2_Boden_Inv_mesh{u}([1 2],k)-P1_Boden_Inv_mesh{i}([1 2],m))));  % Für jeden Wert von P1 zu einem Profil eines Messpunktes von P2
            
    end    
  end
Cross_min{u,i}=min(min(diff));   % gibt das Minimum (wo P1 und P2 sich am nächsten kommen)

[K,M] = find(diff<eee);
row{u,i}=[K];   % K ist indizie von P2
col{u,i}=[M];   % M ist indizie von P1

if length(K)>50
    sort_diff=sort(diff,1);
    sort_diff=sort(sort_diff,2);
    [K,M]= find(diff<=sort_diff(1,10)); 
    %[K,M] = find(diff<=Cross_min{u,i}*2);  
    row{u,i}=[K];   % K ist indizie von P2
    col{u,i}=[M];   % M ist indizie von P1
end
    
    if K~=0 & M~=0
        for zz=1:length(row{u,i})
            delta_z{u,i}(zz)=P2_Boden_Inv{u}([3],row{u,i}(zz))-P1_Boden_Inv{i}([3],col{u,i}(zz));  % gibt das Delta Z
        end        
            Delta_Z{u,i}=sum(delta_z{u,i})/length(row{u,i});
    else
            Delta_Z{u,i}=nan;  % Damit bei PSO nicht 0 als misfit rauskommt
    end      
      
 end   
end

Delta_Z=cell2mat(Delta_Z);
Delta_Z=Delta_Z(~isnan(Delta_Z)); %delete Nan
Delta_Z=sum(abs(Delta_Z))/length(Delta_Z);


