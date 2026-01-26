
% ----------------------------------------------------------------------
function J = jacobian (x, auxdata)
  mpc = auxdata{1};
global regnum loadedData
regID = ['R', num2str(regnum)];
regvar = loadedData.(['nit', regID]);
check         = loadedData.(['check', regID]);
region_buses  = loadedData.(['region_', regID]);
G_itl             = loadedData.(['G', regID]);
B_itl             = loadedData.(['B', regID]);
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);
  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;
  G = mpc.G;
  B = mpc.B;

  J = zeros(2*nbuses+3*nbranches+3*regvar, 2*nbuses+2*ngens+2*regvar);

  k = 0;
  J(k+2*((1:nbuses)-1)+1,1:nbuses) = J(k+2*((1:nbuses)-1)+1,1:nbuses) - diag(2*mpc.bus(:,5)/baseMVA.*x(1:nbuses));
  J(k+2*((1:nbuses)-1)+2,1:nbuses) = J(k+2*((1:nbuses)-1)+2,1:nbuses) + diag(2*mpc.bus(:,6)/baseMVA.*x(1:nbuses));

  for i=1:nbuses
    if length(mpc.genids{i})
      j = mpc.genids{i};
      Pg_count = 2*nbuses + j;
      Qg_count = 2*nbuses + ngens + j;
      J(k+2*(i-1)+1,Pg_count) = J(k+2*(i-1)+1,Pg_count) + 1;
      J(k+2*(i-1)+2,Qg_count) = J(k+2*(i-1)+2,Qg_count) + 1;
    end

    for j=mpc.frombranchids{i}
      if mpc.branch(j,11) == 0; continue; end

      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      J(k+2*(i-1)+1,VMf) = J(k+2*(i-1)+1,VMf) -( 2*G(j,1)*x(VMf) + G(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) + B(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(i-1)+1,VMt) = J(k+2*(i-1)+1,VMt) -( G(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) + B(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(i-1)+1,VAf) = J(k+2*(i-1)+1,VAf) -(-G(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + B(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';
      J(k+2*(i-1)+1,VAt) = J(k+2*(i-1)+1,VAt) -( G(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - B(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';

      J(k+2*(i-1)+2,VMf) = J(k+2*(i-1)+2,VMf) -(-2*B(j,1)*x(VMf) - B(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) + G(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(i-1)+2,VMt) = J(k+2*(i-1)+2,VMt) -(-B(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) + G(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(i-1)+2,VAf) = J(k+2*(i-1)+2,VAf) -( B(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + G(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';
      J(k+2*(i-1)+2,VAt) = J(k+2*(i-1)+2,VAt) -(-B(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - G(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';
    end

    for j=mpc.tobranchids{i}
      if mpc.branch(j,11) == 0; continue; end

      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      J(k+2*(i-1)+1,VMt) = J(k+2*(i-1)+1,VMt) -( 2*G(j,4)*x(VMt) + G(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) + B(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) );
      J(k+2*(i-1)+1,VMf) = J(k+2*(i-1)+1,VMf) -( G(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) + B(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) );
      J(k+2*(i-1)+1,VAt) = J(k+2*(i-1)+1,VAt) -(-G(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + B(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );
      J(k+2*(i-1)+1,VAf) = J(k+2*(i-1)+1,VAf) -( G(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - B(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );

      J(k+2*(i-1)+2,VMt) = J(k+2*(i-1)+2,VMt) -(-2*B(j,4)*x(VMt) - B(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) + G(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) );
      J(k+2*(i-1)+2,VMf) = J(k+2*(i-1)+2,VMf) -(-B(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) + G(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) );
      J(k+2*(i-1)+2,VAt) = J(k+2*(i-1)+2,VAt) -( B(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + G(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );
      J(k+2*(i-1)+2,VAf) = J(k+2*(i-1)+2,VAf) -(-B(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - G(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );
    end
  end


for i=1:regvar


      if check(i)==0 %from

            VMf = region_buses{i};
            VMt = 2*nbuses+2*ngens+2*i-1;
            VAf = region_buses{i}+nbuses;
            VAt = 2*nbuses+2*ngens+2*i;j=i;

      J(k+2*(VMf-1)+1,VMf) = J(k+2*(VMf-1)+1,VMf) -( 2*G_itl(j,1)*x(VMf) + G_itl(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) + B_itl(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(VMf-1)+1,VMt) = J(k+2*(VMf-1)+1,VMt) -( G_itl(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) + B_itl(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(VMf-1)+1,VAf) = J(k+2*(VMf-1)+1,VAf) -(-G_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + B_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';
      J(k+2*(VMf-1)+1,VAt) = J(k+2*(VMf-1)+1,VAt) -( G_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - B_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';

      J(k+2*(VMf-1)+2,VMf) = J(k+2*(VMf-1)+2,VMf) -(-2*B_itl(j,1)*x(VMf) - B_itl(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) + G_itl(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(VMf-1)+2,VMt) = J(k+2*(VMf-1)+2,VMt) -(-B_itl(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) + G_itl(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) )';
      J(k+2*(VMf-1)+2,VAf) = J(k+2*(VMf-1)+2,VAf) -( B_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + G_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';
      J(k+2*(VMf-1)+2,VAt) = J(k+2*(VMf-1)+2,VAt) -(-B_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - G_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) )';

      elseif check(i)==1 %to

            VMt = region_buses{i};
            VMf = 2*nbuses+2*ngens+2*i-1;
            VAt = region_buses{i}+nbuses;
            VAf = 2*nbuses+2*ngens+2*i;j=i;

      J(k+2*(VMt-1)+1,VMt) = J(k+2*(VMt-1)+1,VMt) -( 2*G_itl(j,4)*x(VMt) + G_itl(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) + B_itl(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) );
      J(k+2*(VMt-1)+1,VMf) = J(k+2*(VMt-1)+1,VMf) -( G_itl(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) + B_itl(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) );
      J(k+2*(VMt-1)+1,VAt) = J(k+2*(VMt-1)+1,VAt) -(-G_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + B_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );
      J(k+2*(VMt-1)+1,VAf) = J(k+2*(VMt-1)+1,VAf) -( G_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - B_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );

      J(k+2*(VMt-1)+2,VMt) = J(k+2*(VMt-1)+2,VMt) -(-2*B_itl(j,4)*x(VMt) - B_itl(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) + G_itl(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) );
      J(k+2*(VMt-1)+2,VMf) = J(k+2*(VMt-1)+2,VMf) -(-B_itl(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) + G_itl(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) );
      J(k+2*(VMt-1)+2,VAt) = J(k+2*(VMt-1)+2,VAt) -( B_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + G_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );
      J(k+2*(VMt-1)+2,VAf) = J(k+2*(VMt-1)+2,VAf) -(-B_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - G_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) );


      end


end




  k = 2*nbuses;
  for i=1:nbranches
    if mpc.branch(i,11) == 0; continue; end

    VMf = IDtoCountmap( mpc.branch(i,1) );
    VMt = IDtoCountmap( mpc.branch(i,2) );
    VAf = nbuses + IDtoCountmap( mpc.branch(i,1) );
    VAt = nbuses + IDtoCountmap( mpc.branch(i,2) );

    Pf      =  G(i,1)*x(VMf)^2 + G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    dPfdVMf =  2*G(i,1)*x(VMf) + G(i,2)*x(VMt)*cos(x(VAf)-x(VAt)) + B(i,2)*x(VMt)*sin(x(VAf)-x(VAt));
    dPfdVMt =  G(i,2)*x(VMf)*cos(x(VAf)-x(VAt)) + B(i,2)*x(VMf)*sin(x(VAf)-x(VAt));
    dPfdVAf = -G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    dPfdVAt =  G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));

    Qf      = -B(i,1)*x(VMf)^2 - B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMf = -2*B(i,1)*x(VMf) - B(i,2)*x(VMt)*cos(x(VAf)-x(VAt)) + G(i,2)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMt = -B(i,2)*x(VMf)*cos(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*sin(x(VAf)-x(VAt));
    dQfdVAf =  B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    dQfdVAt = -B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));

    J(k+2*(i-1)+1,VMf) = J(k+2*(i-1)+1,VMf) + 2*Pf*dPfdVMf + 2*Qf*dQfdVMf;
    J(k+2*(i-1)+1,VMt) = J(k+2*(i-1)+1,VMt) + 2*Pf*dPfdVMt + 2*Qf*dQfdVMt;
    J(k+2*(i-1)+1,VAf) = J(k+2*(i-1)+1,VAf) + 2*Pf*dPfdVAf + 2*Qf*dQfdVAf;
    J(k+2*(i-1)+1,VAt) = J(k+2*(i-1)+1,VAt) + 2*Pf*dPfdVAt + 2*Qf*dQfdVAt;

    Pt      =  G(i,4)*x(VMt)^2 + G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMt =  2*G(i,4)*x(VMt) + G(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + B(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMf =  G(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dPtdVAt = -G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dPtdVAf =  G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));

    Qt      = -B(i,4)*x(VMt)^2 - B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMt = -2*B(i,4)*x(VMt) - B(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + G(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMf = -B(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dQtdVAt =  B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dQtdVAf = -B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));

    J(k+2*(i-1)+2,VMt) = J(k+2*(i-1)+2,VMt) + 2*Pt*dPtdVMt + 2*Qt*dQtdVMt;
    J(k+2*(i-1)+2,VMf) = J(k+2*(i-1)+2,VMf) + 2*Pt*dPtdVMf + 2*Qt*dQtdVMf;
    J(k+2*(i-1)+2,VAt) = J(k+2*(i-1)+2,VAt) + 2*Pt*dPtdVAt + 2*Qt*dQtdVAt;
    J(k+2*(i-1)+2,VAf) = J(k+2*(i-1)+2,VAf) + 2*Pt*dPtdVAf + 2*Qt*dQtdVAf;
  end

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

    Pf      =  G_itl(i,1)*x(VMf)^2 + G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    dPfdVMf =  2*G_itl(i,1)*x(VMf) + G_itl(i,2)*x(VMt)*cos(x(VAf)-x(VAt)) + B_itl(i,2)*x(VMt)*sin(x(VAf)-x(VAt));
    dPfdVMt =  G_itl(i,2)*x(VMf)*cos(x(VAf)-x(VAt)) + B_itl(i,2)*x(VMf)*sin(x(VAf)-x(VAt));
    dPfdVAf = -G_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + B_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    dPfdVAt =  G_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - B_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));

    Qf      = -B_itl(i,1)*x(VMf)^2 - B_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMf = -2*B_itl(i,1)*x(VMf) - B_itl(i,2)*x(VMt)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMt = -B_itl(i,2)*x(VMf)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*sin(x(VAf)-x(VAt));
    dQfdVAf =  B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    dQfdVAt = -B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));

    J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VMf) = J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VMf) + 2*Pf*dPfdVMf + 2*Qf*dQfdVMf;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VMt) = J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VMt) + 2*Pf*dPfdVMt + 2*Qf*dQfdVMt;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VAf) = J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VAf) + 2*Pf*dPfdVAf + 2*Qf*dQfdVAf;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VAt) = J(2*nbuses+3*nbranches+(i - 1) * 3 +1,VAt) + 2*Pf*dPfdVAt + 2*Qf*dQfdVAt;

    Pt      =  G_itl(i,4)*x(VMt)^2 + G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMt =  2*G_itl(i,4)*x(VMt) + G_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMf =  G_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dPtdVAt = -G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dPtdVAf =  G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));

    Qt      = -B_itl(i,4)*x(VMt)^2 - B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMt = -2*B_itl(i,4)*x(VMt) - B_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMf = -B_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dQtdVAt =  B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dQtdVAf = -B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));

    J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VMt) = J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VMt) + 2*Pt*dPtdVMt + 2*Qt*dQtdVMt;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VMf) = J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VMf) + 2*Pt*dPtdVMf + 2*Qt*dQtdVMf;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VAt) = J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VAt) + 2*Pt*dPtdVAt + 2*Qt*dQtdVAt;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VAf) = J(2*nbuses+3*nbranches+(i - 1) * 3 +2,VAf) + 2*Pt*dPtdVAf + 2*Qt*dQtdVAf;

end



  k = 2*nbuses + 2*nbranches;
  for i=1:nbranches
    if mpc.branch(i,11) == 0; continue; end

    VAf_count = nbuses + mpc.IDtoCountmap(mpc.branch(i,1));
    VAt_count = nbuses + mpc.IDtoCountmap(mpc.branch(i,2));
    J(k+i,VAf_count) = J(k+i,VAf_count) + 1;
    J(k+i,VAt_count) = J(k+i,VAt_count) - 1;
  end


  for i=1:regvar


      if check(i)==0 %from

            VAf = region_buses{i}+nbuses;
            VAt = 2*nbuses+2*ngens+2*i;

      elseif check(i)==1 %to

            VAt = region_buses{i}+nbuses;
            VAf = 2*nbuses+2*ngens+2*i;

      end


    J(2*nbuses+3*nbranches+(i - 1) * 3 +3 ,VAf) = J(2*nbuses+3*nbranches+(i - 1) * 3 +3 ,VAf) + 1;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +3 ,VAt) = J(2*nbuses+3*nbranches+(i - 1) * 3 +3 ,VAt) - 1;

  end
  

  J = sparse(J);
