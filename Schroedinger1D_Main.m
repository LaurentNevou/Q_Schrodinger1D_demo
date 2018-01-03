clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 23Dec2017, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% model activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

FE_Method=1;            % Diagonalization of the Hamiltonian (FEM)
Euler_Method=1;         % Scanning in Energy method (Euler)
PWE_Method=0;           % Plane Wave Expansion (PWE)

saveV=0;
savePSI=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz=1E-10;               % resolution of the grid [m]
n=10;                   % number of solution asked 

Mass = 0.067;           % effective mass, constant over all the structure...
F0=0%;-1e7;%-20         % Electric field Volt/meter
ScF=0.1;                % scaling factor to plot the wave function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Potential definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two vectors of the same length must be defined
% z [meter]
% V0 [eV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose your between the next 3 potentials or build your own!

Pot_MultiLayers
%Pot_Parabolic
%Pot_Sinus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V0=(F0*z)+V0;           % adding the electric field to the potential
V0=V0-min(V0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1=[];E2=[];E3=[];
display('=======================================')

if FE_Method==1
    tic
    [E1,psi1] = Schroed1D_FEM_f(z,V0,Mass,n);  % m=cste
    display(strcat('-> Finite Elements method =',num2str(toc),'sec'))
end

if Euler_Method==1
    dE=0.005;         % minimum step [eV] between 2 eigen values (bigger => faster)
    precision=1e-5;   % precision [eV] of the eigen values (bigger => faster)
    
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=nan(n,3);
E(1:length(E1),1)=E1;
E(1:length(E2),2)=E2;
E(1:length(E3),3)=E3;

display('=======================================')
display('Results:')
display('=======================================')

display(strcat('E(eV)='))
display(strcat(num2str(E)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 1000 700]);
subplot(1,1,1,'fontsize',15)
hold on;%grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(z*1e9,V0, 'b.-')

if FE_Method==1
  for i=1:length(E1)
      plot(z*1e9,psi1(:,i),'r-')
  end
end

if Euler_Method==1
  for i=1:length(E2)
      plot(z*1e9,psi2(:,i),'k--')
  end
end

if PWE_Method==1
  for i=1:length(E3)
      plot(z*1e9,psi3(:,i),'m--')
  end
end

xlabel('z (nm)');
ylabel('Energy (eV)');

ylim([min(V0)-0.05 max(V0)+0.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Data save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveV==1;
    SS1 = [z' V0' ];
    save('data_V.txt','SS1','-ascii')
end

if savePSI==1;
    
  if FE_Method==1
    SS2 = [z' psi1];
    save('data_PSI_FEM.txt','SS2','-ascii')
  end
  if Euler_Method==1
    SS2 = [z' psi2];
    save('data_PSI_Euler.txt','SS2','-ascii')
  end
  if PWE_Method==1
    SS2 = [z' psi3];
    save('data_PSI_PWE.txt','SS2','-ascii')
  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%