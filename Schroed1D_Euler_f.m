function[E,psi]=Schroed1D_Euler_f(z,V0,Mass,n,dE,precision)

method=2;         % method "2" is much more acurate than "1" but take also more time...
e=min(V0);
Emax=max(V0)+0.1;

C=0; N=0; E=[];
psi=z*0+1; psi_old=psi;


while e<Emax && length(E)<n
    
    C = C+1;
    
    psi_old = psi;
    psi=Schroed1D_Euler_Eval(z,V0,Mass,e,method);
    
    if (sign( psi(end) ) ~= sign( psi_old(end) ) ) && C>1 % here, I catch a quantum state because the last point of psi change sign
    
        N=N+1; de=dE;
        
        while abs(de)>precision
            if sign( psi(end) ) ~= sign( psi_old(end) )
                 de = -de/2;
            end
            e=e+de;
            psi_old=psi;
            psi=Schroed1D_Euler_Eval(z,V0,Mass,e,method);
        end

        E(N,:)=e; PSI(:,N)= psi;
        C=0;
    end	
    
    e=e+dE;
end

psi=PSI;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[psi]=Schroed1D_Euler_Eval(z,V0,Mass,ee,method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [C]
me=9.10938188E-31;              %% electron mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz=z(2)-z(1);
y = [0 1];

if method==1
    
    for i = 1:length(z)                                  %% number of z steps to take
        dy(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * y(1) ; %% Equation for dv/dz
        dy(1) =  Mass*me*y(2);                           %% Equation for dx/dz
        y = y + dz*dy ;                                  %% integrate both equations with Euler
        psi(i)= y(1);
    end
    
elseif method==2
    
    for i = 1:length(z)                                  %% number of z steps to take
    
        dy(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * y(1) ; %% Equation for dv/dz
        dy(1) =  Mass*me*y(2);                           %% Equation for dx/dz
        
        K = y + 0.5*dz*dy;
        dK(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * K(1) ; %% Equation for dv/dz
        dK(1) =  Mass*me*K(2);                           %% Equation for dx/dz
        
        y =  y + dz*dK;
        psi(i)= y(1);
    end
end

end