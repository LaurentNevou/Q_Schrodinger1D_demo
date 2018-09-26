%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vb=0.3;                     % potential barrier [eV]

%zt = [  10 10 10 9 10 8 10 7 10 6 10 ]*1E-9;  % thickness of each layer [nm]
%Vt = [  1 0 1 0 1 0 1 0 1 0 1   ]*Vb;        % barrier height of each layer [eV]

zt = [ 10 10 5 3 10 ]*1E-9;  % thickness of each layer [nm]
Vt = [  1  0 1 0  1  ]*Vb;        % barrier height of each layer [eV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z and the potential V0

for i=1:length(zt)
  if i==1
    zv{1} = 0:dz:zt(1);
    z=zv{1};
    V0=zv{1}*0+Vt(1);
  else
    zv{i} = (z(end)+dz):dz:(z(end) + zt(i));
    z  = [ z  zv{i} ];
    V0 = [ V0   (zv{i}*0+1) * Vt(i)  ];
  end
end

% here, I kind of patch the grid for the TMM

if TM_Method==1
  clear zv
  for i=1:length(zt)
    zzt((i-1)*Ns+1:i*Ns)=zt(i)/Ns;
    VVt((i-1)*Ns+1:i*Ns)=Vt(i);  
  end
  zzt=[1e-12  zzt 1e-12];
  VVt=[10 VVt 10];
  
  for i=1:length(zzt)
    if i==1
      ZZ(1) = zzt(1);
      zv{1} = 0:dz:zzt(1);
      zz    = zv{1};
    else
      ZZ(i) = ZZ(end)+zzt(i);
      zv{i} = (zz(end)+dz):dz:ZZ(end);
      zz    = [zz zv{i}];
    end
  end

end
