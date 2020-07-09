function [p1,p2,alpha1,alpha2,phi1,phi2]=fun_vesseltrack_100m(q,dis_track,dis_pings,head,start_off,vt)

if vt==1;  %Angled profiles

%Profil 1   (green)
px1=start_off:dis_pings:(q*dis_pings+start_off-dis_pings);
py1=ones(1,q);
pz1=zeros(1,q);
p1=[px1;py1;pz1];


%Profil 2    (magenta)
px2=0:dis_pings:(q*dis_pings-dis_pings);
py2=px2*tand(head)+dis_track;
pz2=zeros(1,q);

p2=[px2;py2;pz2];
p2=fliplr(p2);      %track 2 heads in opposite direction

elseif vt==2;  %Parallel Profiles

%Profil 1   (green)
px1=start_off:dis_pings:(q*dis_pings+start_off-dis_pings);
py1=ones(1,q);
pz1=zeros(1,q);
p1=[px1;py1;pz1];

%Profil 2    (magenta)
px2=start_off:dis_pings:(q*dis_pings+start_off-dis_pings);
py2=ones(1,q)+dis_track;
pz2=zeros(1,q);

p2=[px2;py2;pz2];
p2=fliplr(p2);      %track 2 heads in opposite direction

elseif vt==3;   %Sinously Profiles
    
px1=start_off:dis_pings:(q*dis_pings+start_off-dis_pings);
py1=0.05*sin(px1*pi);;    
pz1=zeros(1,q);
p1=[px1;py1;pz1];    

px2=start_off:dis_pings:(q*dis_pings+start_off-dis_pings);
py2=0.05*sin(px2*pi+2)+dis_track;    
pz2=zeros(1,q);
p2=[px2;py2;pz2]; 

    
end

%Calculation of heading and beam angle from vessel track positions
alpha1=nan(1,q-2);
alpha2=nan(1,q-2);
phi1=nan(1,q-2);
phi2=nan(1,q-2);


for i=2:length(p1)-1   
  
    alpha1(i-1)=atan(py1(i+1)-py1(i-1))/(px1(i+1)-px1(i-1)); %Angle of heading
    alpha2(i-1)=atan(py2(i+1)-py2(i-1))/(px2(i+1)-px2(i-1));
    
    phi1(i-1)=alpha1(i-1)+90;   %perpendicular to heading
    phi2(i-1)=alpha2(i-1)-90;   %minus for examiniation of area in between
  
end

end