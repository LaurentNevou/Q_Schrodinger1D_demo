function[E,psi]=Schroed1D_TMM_f(zz,zv,Vt,Mass,n,dE,precision)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Igor A. Sukhoivanov and Igor V. Guryev
% Photonic Cristal
% Chap3: Fundamentals of Computation of Photonic Crystal Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
q=1.602176487E-19;              %% electron charge [C]
m0=9.10938188E-31;              %% electron mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEBUG=0;

e=min(Vt+dE);
Emax=max(Vt);

C=0; N=0; E=[];
b0=0; b0_old=b0;
aN=0; aN_old=aN;
PSI=[];

while e<Emax && length(E)<n
    
    C = C+1;
    
    b0_old = b0;
    aN_old = aN;
    kt=sqrt( 2*Mass*q*m0*(e-Vt) ) / hbar;
    [A,B]=Schroed1D_TMM_Eval(zz,kt,Mass);
    b0 = real(B(1));
    aN = real(A(end));
        
    if ((sign( b0 ) ~= sign( b0_old ) ) && C>1) ||  ((sign( aN ) ~= sign(aN_old ) ) && C>1)
       % here, I catch a quantum state because the last point of psi change sign
    
        N=N+1; de=dE;
        
        while abs(de)>precision
            if (sign( b0 ) ~= sign( b0_old )) || (sign( aN ) ~= sign( aN_old )) 
                 de = -de/2;
            end
            e=e+de;
            if DEBUG==1
              ee(end+1)=e;
            end
            b0_old=b0;
            aN_old=aN;
            
            kt=sqrt( 2*Mass*q*m0*(e-Vt) ) / hbar;
            [A,B]=Schroed1D_TMM_Eval(zz,kt,Mass);
            b0=real(B(1));
            aN = real(A(end));
            
        end
        
        psi=[];
        for j=1:length(Vt)
          psi= [ psi  A(j+1)*exp(1i*kt(j)*zv{j}) + B(j+1)*exp(-1i*kt(j)*zv{j}) ];
        end
        E(N,:)=e; PSI(:,N)= psi;
        C=0;
    end	
    e=e+dE;
    if DEBUG==1
      ee(end+1)=e;
    end
end

psi=PSI;

if DEBUG==1
  figure
  hold on;grid on;
  plot(ee,'b.-')
  xlabel('steps')
  ylabel('Energy (eV)')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[A,B]=Schroed1D_TMM_Eval(zz,kt,Mass)

AmplitudeInput=1;
kL=kt(1);
kR=kt(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left bondary condition

M(1,1:3) = [    1 -1 -1 ];
M(2,1:3) = [ -kL  -kt(1)  kt(1) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filling the matrix

for j=1:length(kt)-1
    M(j*2+1,2*j:2*j+3) = [ exp(1i*kt(j)*zz(j)) exp(-1i*kt(j)*zz(j)) -exp(1i*kt(j+1)*zz(j)) -exp(-1i*kt(j+1)*zz(j))];
    M(2*j+2,2*j:2*j+1) =  kt(j)    *[ +exp(1i*kt(j)  *zz(j))  -exp(-1i*kt(j)  *zz(j))  ];
    M(2*j+2,2*j+2:2*j+3) = -kt(j+1)  *[ +exp(1i*kt(j+1)*zz(j))  -exp(-1i*kt(j+1)*zz(j))  ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right bondary condition

M(length(kt)*2+1,2*length(kt):2*length(kt)+2) = [exp(1i*kt(end)*zz(end)) exp(-1i*kt(end)*zz(end)) -exp(1i*kR*zz(end)) ];
M(length(kt)*2+2,2*length(kt):2*length(kt)+2) = [ kt(end)*exp(1i*kt(end)*zz(end))  -kt(end)*exp(-1i*kt(end)*zz(end)) -kR*exp(1i*kR*zz(end)) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=sparse(M);

D=zeros(length(M),1);
D(1)=-sqrt(AmplitudeInput);
D(2)=-sqrt(AmplitudeInput)*kL;

AB=inv(M)*D;
A=[1 ; AB(2:2:end)];
B=[AB(1:2:end-1) ; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end