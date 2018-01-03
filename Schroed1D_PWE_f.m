function[E,psi]=Schroed1D_PWE_f(z,V0,Mass,n,Nz,NG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [C]
me=9.10938188E-31;              %% electron mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Creation of the inex Geometrie %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NG = 2*floor(NG/2);       % round to lower even number

zz=linspace(z(1),z(end),Nz);
V=interp1(z,V0,zz);

dz=z(2)-z(1);
dzz=zz(2)-zz(1);
Ltot=zz(end)-zz(1);

Gamma=1./(Mass+0*zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Building of the potential in Fourier space %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vk = fft(V)*dzz;
Vk=fftshift(Vk);
%size(Vk)
Vk =Vk(Nz/2-NG+1:Nz/2+NG+1);
%size(Vk)

Gammak = fft(Gamma)*dzz;
Gammak=fftshift(Gammak);
Gammak = Gammak(Nz/2-NG+1:Nz/2+NG+1);
%size(Gammak)

%V_back = invFFT1D_Ver2(Vk,Nz)/dzz ;   % fft to have back the compressed potential
%
%figure
%hold on;grid on;
%plot(zz,V,'b.-')
%plot(zz,V_back,'r.-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Initialisation of Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = (-NG/2:NG/2)*2*pi/Ltot;	        % reciprocal lattice vectors
NG=length(G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Building Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build H matrix
% gammaK and VK are spectra of real valued functions, thus they
% are symmetric VK(n) == VK(-n) == VK(abs(n))
%
% we have to compute a matrix such that
% H(i,j) = (G(i) + k)*(G(j) + k)*gammaK(abs(i - j)+1) + VK(i - j + N)
%  

idx = ones(NG,1)*(1:NG);    % matrix with elements equal to the column index idx(i,j) = j
idx = idx - idx';           % matrix with elements equal to row - column idx(i,j) = i - j;
idx = idx + NG;             % idx(i,j) = i - j + N
        
HV = Vk(idx)*e/(Ltot*hbar^2/(2*me)); % H(i,j) =  VK(i - j + N)
% idx = abs(idx - N) + 1; % idx(i,j) = abs(i - j) + 1
Gk = G'*G; % Gk(i,j) = (G(i) + k)*(G(j) + k)
H = HV +  Gk.*Gammak(idx)/Ltot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Solving Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = sparse(H);
[psik, Ek] = eigs(H,n,'SM');
E = diag(Ek) * hbar^2/(2*me) / e;
E=abs(E);
%E=reel(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Scaling of the waves functions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:n
    psizz(:,j) = invFFT1D_Ver2(psik(:,j)',Nz)/dzz ;
    psi(:,j)   = interp1(zz,psizz(:,j),z);
    psi(:,j)   = psi(:,j)/sqrt(sum(abs(psi(:,j)).^2)*dz) ;  % normalisation at 1
end

psi=psi(:,end:-1:1);
E=E(end:-1:1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vz] = invFFT1D_Ver2(Vk,Nz)


Nk=length(Vk);

Nxx=Nz/2-ceil(Nk/2);

Vk1=[zeros(1,Nxx+1)  Vk  zeros(1,Nxx)];

Vk2=ifftshift(Vk1);

Vz=ifft(Vk2);

end