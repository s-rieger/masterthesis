%% RUN PSO

clear all;
close all;
clc;
tic

%% Initial Model
v_ref=[1.5000 1.4440 1.4407 1.4590 1.4790 1.4958 1.5125 1.5292 1.5363 1.5430 1.5497]*1e3;   %Velocity
M_ref=[10 10 10 10 10 10 10 10 10 10];                %Thickness
d_v_ref=[0 30 30 30 30 30 30 30 30 30 30];            %Velocity errors
d_M_ref=[1 1 1 1 1 1 1 1 1 1];              %Thickness errors

SVP_init=[M_ref,v_ref];

dlmwrite('SVP_init',SVP_init);
dlmwrite('SVP_init_v_err',d_v_ref);
dlmwrite('SVP_init_M_err',d_M_ref);

%% Data from Forward
mk_sq=load('mk_sq');                      %Ratio of X-/Y-Axis for sqauring grid for intersection calculations
info=load('theta_info');                  %Format of Data

%% PSO calibration
n=15;         %Swarmsize
omega=0.5;	 %Weight of initial momentum
c1=0.9;		 %Weight of individual best
c2=0.5;		 %Weight of global best
nmax=200;     %maximum nr. of iterations 100/10
Rmax=10;      %Do reset of swarm after Rmax unsuccessful iterations 10
ee=6.0e-4;	 %An iteration is called unsuccessful if the change in misfit is smaller than ee
seednr=5000; %Seed of random number generator 1500 (aussen)

%% Inversion Parameters
eee=1*mk_sq/10;       %Parameter to determine intersections between rays [x y] (squared Grid! acutally [*1/mk_sq])
dz=sum(M_ref)*0.005;  %Incremention of depth for SVP/Raytrace 0.5%

%% Model Parameters
lin_stp=3; %combo 
no_vt=3;
no_type=5;


 for jjj=1:no_vt;  %1:angled 2:parallel 3:sinous
  for iii=1:no_type;

        
        vt=jjj;
        type=iii;              
        
        fprintf('Seafloor:%d VT:%d \n',iii,jjj);

        P1_Data=load(sprintf('P1_Data_SF_%d_VT_%d',type,vt));          % Eigentlich muss ja nur 1 SVP model aus vorwärts reichen, da man das ja als "richitg" annimmt! die PSO von dem forward modell dann aber sowohl lin als auch stepped
        P2_Data=load(sprintf('P2_Data_SF_%d_VT_%d',type,vt)); 

        L=ones([1,info(length(info))])*info(1);  %For splitting Data into Pings
        P1_data=mat2cell(P1_Data,L);
        P2_data=mat2cell(P2_Data,L);


%% Run PSO-Inversion
results=fun_PSO_100m(n,omega,c1,c2,nmax,Rmax,ee,seednr,v_ref,M_ref,d_v_ref,d_M_ref,P1_data,P2_data,eee,mk_sq,type,lin_stp,dz,vt);
dlmwrite(sprintf('results_SF_%d_VT_%d',type,vt),results);
clear P1_Data P2_Data

  end
end
toc
