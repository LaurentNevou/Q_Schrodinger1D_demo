function[E,psi]=Schroed1D_Euler_f(z,V0,Mass,n,dE,precision)

method=2;         % method 2 is more acurate than the 1
Emin=min(V0);
Emax=max(V0)+0.1;
e=Emin;

Evec=Emin:dE:Emax;

CCC=0;
NNN=0;
E=[];

psi=z*0;
psi(end)=1;
psi_old=psi;


while e<Emax && length(E)<n

    CCC = CCC+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx=[0 0];
    epsi=1;
    while length(idx)>1
        idx=find( abs( (Evec-e )) < epsi ) ;
      if isempty(idx)
          idx=IDX(1);
        else
          IDX=idx;
          epsi=epsi/2;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    psi_old = psi;
    psi=Schrod_Euler_Eval_v03(z,V0,Mass,e,method);
    
    if (sign( psi(end) ) ~= sign( psi_old(end) ) ) && CCC>1
    
        NNN=NNN+1; de=dE;
        
        while abs(de)>precision
            if sign( psi(end) ) ~= sign( psi_old(end) )
                 de = -de/2;
            end
            e=e+de;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            idx=[0 0];
            epsi=1;
            while length(idx)>1
                idx=find( abs( (Evec-e )) < epsi ) ;
              if isempty(idx)
                  idx=IDX(1);
                else
                  IDX=idx;
                  epsi=epsi/2;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            psi_old=psi;
            psi=Schrod_Euler_Eval_v03(z,V0,Mass,e,method);
        end

        E(NNN,:)=e;
        PSI(:,NNN)= psi;
        CCC=0;
    end	
    
    e=e+dE;
end

psi=PSI;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[psi]=Schrod_Euler_Eval_v03(z,V0,Mass,ee,method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [C]
me=9.10938188E-31;              %% electron mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz=z(2)-z(1);
psi_E=z*0;
y = [0 1];

if method==1

    for i = 1:length(z)                                 %% number of time steps to take
        dy(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * y(1) ; %% Equation for dv/dt
        dy(1) =  Mass*me*y(2);                          %% Equation for dx/dt
        %y = y + dz*dy ;                                %% integrate both equations with Euler
        y = y + dz*dy ;
        psi(i)= y(1);
    end

elseif method==2

    for i = 1:length(z)

        dy(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * y(1) ; %% Equation for dv/dt
        dy(1) =  Mass*me*y(2);                          %% Equation for dx/dt
      
        %K = y + 0.5*dz*dy;
        K = y + 0.5*dz*dy;
        dK(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * K(1) ; %% Equation for dv/dt
        dK(1) =  Mass*me*K(2);                          %% Equation for dx/dt
    
        %y =  y + dz*dK;
        y =  y + dz*dK;
        psi(i)= y(1);
    end
end

end