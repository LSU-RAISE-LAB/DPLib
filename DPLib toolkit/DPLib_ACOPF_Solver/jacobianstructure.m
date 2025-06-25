
% ----------------------------------------------------------------------
function J = jacobianstructure (auxdata)
  mpc = auxdata{1};
global regnum loadedData
regID = ['R', num2str(regnum)];
regvar = loadedData.(['nit', regID]);
check         = loadedData.(['check', regID]);
region_buses  = loadedData.(['region_', regID]);
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);
  IDtoCountmap = mpc.IDtoCountmap;

  J = zeros(2*nbuses+3*nbranches+3*regvar, 2*nbuses+2*ngens+2*regvar);

  k = 0;
  for i=1:nbuses
    J(k+2*(i-1)+1,i) = J(k+2*(i-1)+1,i) || (mpc.bus(i,5)~=0);
    J(k+2*(i-1)+2,i) = J(k+2*(i-1)+2,i) || (mpc.bus(i,6)~=0);

    for genid=mpc.genids{i}
      if mpc.gen(genid,8)==0; continue; end

      Pg_count = 2*nbuses + genid;
      Qg_count = 2*nbuses + ngens + genid;
      J(k+2*(i-1)+1,Pg_count) = 1;
      J(k+2*(i-1)+2,Qg_count) = 1;
    end

    for j=mpc.frombranchids{i}
      if mpc.branch(j,11) == 0; continue; end

      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      J(k+2*(i-1)+1,[VMf VMt VAf VAt]) = 1;
      J(k+2*(i-1)+2,[VMf VMt VAf VAt]) = 1;
    end

    for j=mpc.tobranchids{i}
      if mpc.branch(j,11) == 0; continue; end

      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      J(k+2*(i-1)+1,[VMf VMt VAf VAt]) = 1;
      J(k+2*(i-1)+2,[VMf VMt VAf VAt]) = 1;
    end
  end


for i=1:regvar


      if check(i)==0 %from

            VMf = region_buses{i};
            VMt = 2*nbuses+2*ngens+2*i-1;
            VAf = region_buses{i}+nbuses;
            VAt = 2*nbuses+2*ngens+2*i;

      J(k+2*(VMf-1)+1,[VMf VMt VAf VAt]) = 1;
      J(k+2*(VMf-1)+2,[VMf VMt VAf VAt]) = 1;

      elseif check(i)==1 %to

            VMt = region_buses{i};
            VMf = 2*nbuses+2*ngens+2*i-1;
            VAt = region_buses{i}+nbuses;
            VAf = 2*nbuses+2*ngens+2*i;

      J(k+2*(VMt-1)+1,[VMf VMt VAf VAt]) = 1;
      J(k+2*(VMt-1)+2,[VMf VMt VAf VAt]) = 1;


      end


end


  k = 2*nbuses;
  for i=1:nbranches
    if mpc.branch(i,11) == 0; continue; end

    VMf = IDtoCountmap( mpc.branch(i,1) );
    VMt = IDtoCountmap( mpc.branch(i,2) );
    VAf = nbuses + IDtoCountmap( mpc.branch(i,1) );
    VAt = nbuses + IDtoCountmap( mpc.branch(i,2) );

    J(k+2*(i-1)+1,[VMf VMt VAf VAt]) = 1;
    J(k+2*(i-1)+2,[VMf VMt VAf VAt]) = 1;
  end

  k = 2*nbuses + 2*nbranches;
  for i=1:nbranches
    if mpc.branch(i,11) == 0; continue; end

    VAf_count = nbuses + mpc.IDtoCountmap(mpc.branch(i,1));
    VAt_count = nbuses + mpc.IDtoCountmap(mpc.branch(i,2));
    J(k+i,[VAf_count VAt_count]) = 1;
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

    J(2*nbuses+3*nbranches+(i - 1) * 3 +1,[VMf VMt VAf VAt]) = 1;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +2,[VMf VMt VAf VAt]) = 1;
    J(2*nbuses+3*nbranches+(i - 1) * 3 +3,[VAf VAt]) = 1;
end

  J = sparse(J);
