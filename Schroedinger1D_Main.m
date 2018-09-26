%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% last update 26Sept2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
warning off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Model activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

FE_Method=1;            % Diagonalization of the Hamiltonian (FEM)
Euler_Method=0;         % Scanning in Energy method (Euler)
PWE_Method=0;           % Plane Wave Expansion (PWE)
TM_Method=1;Ns=3;       % Transfer Matrix Method (TMM), Ns=number of segment per layer !!!WORK ONLY with Pot_MultiLayers!!!

saveV=0;
savePSI=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz=2E-10;               % resolution of the grid [m]
n=4;                   % number of solution asked 

Mass = 0.067;           % effective mass, constant over all the structure...
F0=0;%-1e7;             % Electric field [V/m]
ScF=0.1;                % scaling factor to plot the wave function [Without Dimension]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Potential definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two vectors of the same length must be defined
% z [meter]
% V0 [eV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose your between the next 3 potentials or build your own!

Pot_MultiLayers
%Pot_Parabolic
%Pot_Sinus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0=(F0*z)+V0;           % adding the electric field to the potential
V0=V0-min(V0);

if TM_Method==1
  VVt=(F0*ZZ)+VVt;      % adding the electric field to the potential
  VVt=VVt-min(VVt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1=[];E2=[];E3=[];E4=[];
display('=======================================')

if FE_Method==1
    tic
    [E1,psi1] = Schroed1D_FEM_f(z,V0,Mass,n);  % m=cste
    display(strcat('-> Finite Elements method =',num2str(toc),'sec'))
end

if Euler_Method==1
    dE=0.005;         % minimum step [eV] between 2 eigen values (bigger => faster)
    precision=1e-7;   % precision [eV] of the eigen values (bigger => faster)
    
    tic
    [E2,psi2] = Schroed1D_Euler_f(z,V0,Mass,n,dE,precision);
    display(strcat('-> Euler method =',num2str(toc),'sec'))
end

if PWE_Method==1
    Nz = 256 ;        % number of points on the z grid % has to be a power of 2 (32,64,128,256,512,...) (smaller => faster)
    NG = Nz/2-1  ;    % number of harmonics % has to be at least 2 times -1 smaller than Nz (smaller => faster)
    
    tic
    [E3,psi3] = Schroed1D_PWE_f(z,V0,Mass,n,Nz,NG);
    display(strcat('-> PWE method =',num2str(toc),'sec'))
    if F0~=0
      display('Warning: The PWE method is not appropriate for non-periodic potential')
    end
end

if TM_Method==1
    
    dE=0.001;          % minimum step [eV] between 2 eigen values (bigger => faster)
    precision=1e-16;   % precision [eV] of the eigen values (bigger => faster)
    
    tic
    [E4,psi4] = Schroed1D_TMM_f(ZZ,zv,VVt,Mass,n,dE,precision);
    display(strcat('-> TMM method =',num2str(toc),'sec'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=nan(n,4);
E(1:length(E1),1)=E1;
E(1:length(E2),2)=E2;
E(1:length(E3),3)=E3;
E(1:length(E4),4)=E4;

display('=======================================')
display('Results:')
display('=======================================')

display(strcat('E(eV)='))
display(strcat(num2str(E)))

%break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FE_Method==1
  for i=1:length(E1)
      psi1(:,i)=abs(psi1(:,i)).^2/max(abs(psi1(:,i)).^2)*ScF + E1(i); % normalisation for the plotting
  end
end

if Euler_Method==1
  for i=1:length(E2)
      psi2(:,i)=abs(psi2(:,i)).^2/max(abs(psi2(:,i)).^2)*ScF + E2(i); % normalisation for the plotting
  end
end

if PWE_Method==1
  for i=1:length(E3)
      psi3(:,i)=abs(psi3(:,i)).^2/max(abs(psi3(:,i)).^2)*ScF + E3(i); % normalisation for the plotting
  end
end

if TM_Method==1
  for i=1:length(E4)
      psi4(:,i)=abs(psi4(:,i)).^2/max(abs(psi4(:,i)).^2)*ScF + E4(i); % normalisation for the plotting
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[10 100 1000 700]);
LW=2;
FS=15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;%grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(z*1e9,V0, 'b-','linewidth',LW)

s{1}='\fontsize{20}\color{blue}Potential';
s{2}='\fontsize{20}\color{blue}----------';

if FE_Method==1
  for i=1:length(E1)
      plot(z*1e9,psi1(:,i),'r-','linewidth',LW)
  end
  s{end+1}='\fontsize{20}\color{red}FE-Method';
end

if Euler_Method==1
  for i=1:length(E2)
      plot(z*1e9,psi2(:,i),'k--','linewidth',LW)
  end
  s{end+1}='\fontsize{20}\color{black}Euler-Method';
end

if PWE_Method==1
  for i=1:length(E3)
      plot(z*1e9,psi3(:,i),'m--','linewidth',LW)
  end
  s{end+1}='\fontsize{20}\color{magenta}PWE-Method';
end

if TM_Method==1
  for i=1:length(E4)
      plot(z*1e9,psi4(:,i),'g--','linewidth',LW)
  end
  s{end+1}='\fontsize{20}\color{green}TM-Method';
end

xlabel('z (nm)');
ylabel('Energy (eV)');

ylim([min(V0)-0.05 max(V0)+0.1])
title('Multi layers potential')

xl=get(gca,'xlim');
yl=get(gca,'ylim');

text(0.02*xl(2),yl(2),s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveV==1;
    MV = [z' V0' ];
    save('data_V.txt','MV','-ascii')
end

if savePSI==1;
    
  if FE_Method==1
    M1 = [z' psi1];
    save('data_PSI_FEM.txt','M1','-ascii')
  end
  if Euler_Method==1
    M2 = [z' psi2];
    save('data_PSI_Euler.txt','M2','-ascii')
  end
  if PWE_Method==1
    M3 = [z' psi3];
    save('data_PSI_PWE.txt','M3','-ascii')
  end
  if TM_Method==1
    M4 = [z' psi4];
    save('data_PSI_TMM.txt','M4','-ascii')
  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%