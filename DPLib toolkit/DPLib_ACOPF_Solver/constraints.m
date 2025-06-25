
% ----------------------------------------------------------------------
function c = constraints (x, auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);
  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;
  G = mpc.G;
  B = mpc.B;
global regnum loadedData
regID = ['R', num2str(regnum)];
regvar = loadedData.(['nit', regID]);
check         = loadedData.(['check', regID]);
region_buses  = loadedData.(['region_', regID]);
G_itl             = loadedData.(['G', regID]);
B_itl             = loadedData.(['B', regID]);

  c = zeros(2*nbuses + 3*nbranches+3*regvar,1);

  VM = 1:nbuses;
  c(1:2:2*nbuses) = -mpc.bus(:,5)/baseMVA.*x(VM).^2;
  c(2:2:2*nbuses) =  mpc.bus(:,6)/baseMVA.*x(VM).^2;

  k = 0;
  for i=1:nbuses
    if length(mpc.genids{i})
      Pg_count = 2*nbuses + mpc.genids{i};
      Qg_count = 2*nbuses + ngens + mpc.genids{i};
      c(k+2*(i-1)+1) = c(k+2*(i-1)+1) + sum(x(Pg_count));
      c(k+2*(i-1)+2) = c(k+2*(i-1)+2) + sum(x(Qg_count));
    end

    if length(mpc.frombranchids{i})
      j = mpc.frombranchids{i};
      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      Pf =  G(j,1).*x(VMf).^2 + G(j,2).*x(VMf).*x(VMt).*cos(x(VAf)-x(VAt)) + B(j,2).*x(VMf).*x(VMt).*sin(x(VAf)-x(VAt));
      Qf = -B(j,1).*x(VMf).^2 - B(j,2).*x(VMf).*x(VMt).*cos(x(VAf)-x(VAt)) + G(j,2).*x(VMf).*x(VMt).*sin(x(VAf)-x(VAt));

      c(k+2*(i-1)+1) = c(k+2*(i-1)+1) - sum(Pf);
      c(k+2*(i-1)+2) = c(k+2*(i-1)+2) - sum(Qf);
    end

    if length(mpc.tobranchids{i})
      j = mpc.tobranchids{i};
      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      Pt =  G(j,4).*x(VMt).^2 + G(j,3).*x(VMt).*x(VMf).*cos(x(VAt)-x(VAf)) + B(j,3).*x(VMt).*x(VMf).*sin(x(VAt)-x(VAf));
      Qt = -B(j,4).*x(VMt).^2 - B(j,3).*x(VMt).*x(VMf).*cos(x(VAt)-x(VAf)) + G(j,3).*x(VMt).*x(VMf).*sin(x(VAt)-x(VAf));

      c(k+2*(i-1)+1) = c(k+2*(i-1)+1) - sum(Pt);
      c(k+2*(i-1)+2) = c(k+2*(i-1)+2) - sum(Qt);
    end
  end

  k = 2*nbuses;
  j = find(mpc.branch(:,11));

  VMf = IDtoCountmap( mpc.branch(j,1) );
  VMt = IDtoCountmap( mpc.branch(j,2) );
  VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
  VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

  Pf(j) =  G(j,1).*x(VMf).^2 + G(j,2).*x(VMf).*x(VMt).*cos(x(VAf)-x(VAt)) + B(j,2).*x(VMf).*x(VMt).*sin(x(VAf)-x(VAt));
  Qf(j) = -B(j,1).*x(VMf).^2 - B(j,2).*x(VMf).*x(VMt).*cos(x(VAf)-x(VAt)) + G(j,2).*x(VMf).*x(VMt).*sin(x(VAf)-x(VAt));
  Pt(j) =  G(j,4).*x(VMt).^2 + G(j,3).*x(VMt).*x(VMf).*cos(x(VAt)-x(VAf)) + B(j,3).*x(VMt).*x(VMf).*sin(x(VAt)-x(VAf));
  Qt(j) = -B(j,4).*x(VMt).^2 - B(j,3).*x(VMt).*x(VMf).*cos(x(VAt)-x(VAf)) + G(j,3).*x(VMt).*x(VMf).*sin(x(VAt)-x(VAf));

  c(k+2*(j-1)+1) = Pf(j).^2 + Qf(j).^2;
  c(k+2*(j-1)+2) = Pt(j).^2 + Qt(j).^2;

  k = 2*nbuses + 2*nbranches;
  VAf_count = nbuses + mpc.IDtoCountmap(mpc.branch(:,1));
  VAt_count = nbuses + mpc.IDtoCountmap(mpc.branch(:,2));
  c(k+(1:nbranches)) = mpc.branch(:,11).*( x(VAf_count) - x(VAt_count) );





  for i=1:regvar


      if check(i)==0 %from

            VMf = region_buses{i};
            VMt = 2*nbuses+2*ngens+2*i-1;
            VAf = region_buses{i}+nbuses;
            VAt = 2*nbuses+2*ngens+2*i;

      elseif check(i)==1 %to

            VMt = region_buses{i};
            VMf = 2*nbuses+2*ngens+2*i-1;
            VAt = region_buses{i}+nbuses;
            VAf = 2*nbuses+2*ngens+2*i;

      end

  Pfv =  G_itl(i,1).*x(VMf).^2 + G_itl(i,2).*x(VMf).*x(VMt).*cos(x(VAf)-x(VAt)) + B_itl(i,2).*x(VMf).*x(VMt).*sin(x(VAf)-x(VAt));
  Qfv = -B_itl(i,1).*x(VMf).^2 - B_itl(i,2).*x(VMf).*x(VMt).*cos(x(VAf)-x(VAt)) + G_itl(i,2).*x(VMf).*x(VMt).*sin(x(VAf)-x(VAt));
  Ptv =  G_itl(i,4).*x(VMt).^2 + G_itl(i,3).*x(VMt).*x(VMf).*cos(x(VAt)-x(VAf)) + B_itl(i,3).*x(VMt).*x(VMf).*sin(x(VAt)-x(VAf));
  Qtv = -B_itl(i,4).*x(VMt).^2 - B_itl(i,3).*x(VMt).*x(VMf).*cos(x(VAt)-x(VAf)) + G_itl(i,3).*x(VMt).*x(VMf).*sin(x(VAt)-x(VAf));


      if check(i)==0 %from

  c(2*(VMf-1)+1) = c(2*(VMf-1)+1) - Pfv;
  c(2*(VMf-1)+2) = c(2*(VMf-1)+2) - Qfv;

      elseif check(i)==1 %to

  c(2*(VMt-1)+1) = c(2*(VMt-1)+1) - Ptv;
  c(2*(VMt-1)+2) = c(2*(VMt-1)+2) - Qtv;

      end

c(2*nbuses+3*nbranches+(i - 1) * 3 +1) =Pfv.^2 + Qfv.^2;
c(2*nbuses+3*nbranches+(i - 1) * 3 +2) =Ptv.^2 + Qtv.^2;
c(2*nbuses+3*nbranches+(i - 1) * 3 +3) =x(VAf) - x(VAt);



  end


