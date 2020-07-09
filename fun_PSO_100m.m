function [results]=fun_PSO_100m(n,omega,c1,c2,nmax,Rmax,ee,seednr,v_ref,M_ref,d_v_ref,d_M_ref,P1_data,P2_data,eee,mk_sq,type,lin_stp,dz,vt)
%% PSO-I. Initialize PSO:------------------------------------------------------------------------------------------
%#####1. Create Initial Particles:
m=length(v_ref)+length(M_ref); %problem dimension
l=length(M_ref); %nr of layers above halfspace


misfit=zeros(n,1);ibmisfit=zeros(n,1);particle=zeros(m,n);ibparticle=zeros(m,n);
gbmisfit=10000.0;
gbmisfit_old=10000.0;
gbparticle=zeros(m,1);
particlev=zeros(m,n);
rng(seednr,'twister');% initialize mersenne twister MATLAB
%rand("state",seednr); % initialize mersenne twister OCTAVE

for j=1:n
%randomize:
particle(1:l,j)=M_ref+(2.*rand(1,l)-1).*d_M_ref;
particle(l+1:m,j)=v_ref+(2.*rand(1,m-l)-1).*d_v_ref;

%calculate misfit:
misfit(j)=fun_delta_100m((particle(:,j)'),P1_data,P2_data,eee,mk_sq,lin_stp,dz,vt);

% Dont display the polynomial bad conditioned warining
[msg, id] = lastwarn;
 
warning('off', id)   
end
%#####2. Set individual best:
ibparticle=particle;
ibmisfit=misfit;
%#####3. Find global best:
[gbmisfit,gbindex]=min(misfit);   
%[gbmisfit,gbindex]=find(misfit==0); % Hier nicht suchen wo 0 ist?!
gbparticle=particle(:,gbindex);

%% PSO-II. Main Loop:----------------------------------------------------------------------------------------------
printcount=0;
Reset=0;
resetcount=0;
no_reset=0;
for main=1:nmax %MAINLOOP 
printcount=printcount+1;
fprintf('DOING ITERATION %d\n',main);
%#####1. Alternate Particles randomly:
if Reset==1
 %reinitialze particles and individual best and global best, printout so far global best, set
 printcount=0;
 disp('RESETTING SWARM')
 %disp('So far all best:')
 %disp(gbparticle) %ACHTUNG ABSPEICHERN mit misfit
 no_reset=no_reset+1;
 for j=1:n %PARTICLES
   %randomize:
   particle(1:l,j)=M_ref+(2.*rand(l,1)-1)'.*d_M_ref;
   particle(l+1:m,j)=v_ref+(2.*rand(m-l,1)-1)'.*d_v_ref;
   
      sumoflayer=sum(particle(1:l,j));
   while sumoflayer<sum(M_ref)-M_ref(1) ;%& sumoflayer>sum(M_ref)-M_ref(1);   %to prevent layers of getting too small
        r1=rand(m,1);
        r2=rand(m,1);
        particlev(:,j)=omega*particlev(:,j) + c1*(r1.*(ibparticle(:,j)-particle(:,j)))+c2*(r2.*(gbparticle-particle(:,j)));
        particle(:,j)=particle(:,j)+particlev(:,j);
        sumoflayer=sum(particle(1:l,j));
   end
   %calculate misfit:
   misfit(j)=fun_delta_100m(particle(:,j),P1_data,P2_data,eee,mk_sq,lin_stp,dz,vt);
 end
 
  %% ...Set individual best...

  if misfit(j)<ibmisfit(j)
      ibparticle(:,j)=particle(:,j);
      ibmisfit(j)=misfit(j);
  end
  
  %% ...Find global best...
  [gbmisfit,gbindex]=min(misfit);
  gbparticle=particle(:,gbindex);
  Reset=0;
 
else
 
  for j=1:n %PARTICLES update
  %fprintf('updating particle %d\n',j);
  r1=rand(m,1);
  r2=rand(m,1);
  particlev(:,j)=omega*particlev(:,j) + c1*(r1.*(ibparticle(:,j)-particle(:,j)))+c2*(r2.*(gbparticle-particle(:,j)));
  particle(:,j)=particle(:,j)+particlev(:,j);
 
  sumoflayer=sum(particle(1:l,j));
  sumticker=0;
  while sumoflayer<sum(M_ref)-M_ref(1) ;%& sumoflayer>sum(M_ref)-M_ref(1);   %to prevent layers of getting too small
        r1=rand(m,1);
        r2=rand(m,1);
        particlev(:,j)=omega*particlev(:,j) + c1*(r1.*(ibparticle(:,j)-particle(:,j)))+c2*(r2.*(gbparticle-particle(:,j)));
        particle(:,j)=particle(:,j)+particlev(:,j);
        sumoflayer=sum(particle(1:l,j));
        sumticker=sumticker+1;
        if sumticker>1000;
           particle(1:l,j)=M_ref+(2.*rand(l,1)-1)'.*d_M_ref;
           particle(l+1:m,j)=v_ref+(2.*rand(m-l,1)-1)'.*d_v_ref; 
           reset=1;
           break
        end
  end
  
    if find(particle(1:l,j)<=0)
      particle(1:l,j)=M_ref+(2.*rand(1,l)-1).*d_M_ref;
      reset=1;
	end
  %% ...calculate misfit...
  print=0;
  plotF=0;
  misfit(j)=fun_delta_100m(particle(:,j),P1_data,P2_data,eee,mk_sq,lin_stp,dz,vt);
    
  end
   
end

%#####2. Determine individual best:
for j=1:n
    if misfit(j)<ibmisfit(j)
        ibparticle(:,j)=particle(:,j);
        ibmisfit(j)=misfit(j);
    end
end
%#####3. Find global best:
[gbmisfit_ac,gbindex]=min(misfit);

%######4. ...GCPSO...
  rho=[2,0.5,1];
  gbparticle_gcpso=nan(length(particle(:,1)),1);
  
for b=1:length(rho)
    
    r1=rand(m,1);    
    gbparticle_gcpso=gbparticle+omega*particlev(:,gbindex)+(1-2*r1)*rho(b);  
    sumoflayer=sum(gbparticle_gcpso(1:l));
    sumticker=0;
    
    while sumoflayer<sum(M_ref)-M_ref(1) ;%& sumoflayer>sum(M_ref)-M_ref(1);   %to prevent layers of getting too small
        
        r1=rand(m,1);
        gbparticle_gcpso=gbparticle+omega*particlev(:,gbindex)+(1-2*r1)*rho(b);
        sumoflayer=sum(gbparticle_gcpso(1:l));
        sumticker=sumticker+1;
        if sumticker>1000;
            gbparticle_gcpso=M_ref+(2.*rand(1,l)-1).*d_M_ref;
           
            break
        end     
    end
    
    if find(gbparticle_gcpso(1:l)<=0)
      gbparticle_gcpso(1:l)=M_ref+(2.*rand(1,l)-1).*d_M_ref;
    end
    
    misfit_gcpso=fun_delta_100m(gbparticle_gcpso,P1_data,P2_data,eee,mk_sq,lin_stp,dz,vt);
    
    if misfit_gcpso<gbmisfit
        gbparticle=gbparticle_gcpso;
        gbmisfit=misfit_gcpso;
    end
end
  
%%%% conditions
if gbmisfit_ac<gbmisfit
 if abs((gbmisfit_ac-gbmisfit_old)/gbmisfit_ac) <= ee
  resetcount=resetcount+1;
  fprintf('Similar minimum. Resetcount=%d \n',resetcount)
 else
  resetcount=0;
 end
 if resetcount==Rmax 
  Reset=1;
  resetcount=0;
 end
 fprintf('Found better misfit %f with:\n',gbmisfit_ac)
 gbmisfit_old=gbmisfit_ac;
 gbparticle=particle(:,gbindex);
 gbmisfit=misfit(gbindex);
 

else
%#####4. Reset?: (if yes, set Reset=1)
 resetcount=resetcount+1;
 fprintf('Unsuccessful step. Resetcount=%d \n',resetcount)
 gbmisfit_old=gbmisfit;
end


gbparticles_all{main}=[gbparticle;gbmisfit_ac;main;no_reset]; %misfit is at bottom of colom
  
if resetcount==Rmax  
Reset=1;
resetcount=0;
%save all the best global particles with misfit

end

end%MAINLOOP

%results=[particle;misfit];
 results=gbparticles_all(~cellfun('isempty',gbparticles_all)); % Löschen aller leeren Zellen
 results=cell2mat(results);
end
