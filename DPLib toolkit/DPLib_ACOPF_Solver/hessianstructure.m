
% ----------------------------------------------------------------------
function H = hessianstructure (auxdata)
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
  H = zeros(2*nbuses+2*ngens+2*regvar, 2*nbuses+2*ngens+2*regvar);
if ngens~=0
  for i=1:ngens
    if mpc.gen(i,8)==0; continue; end

    Pg_count = 2*nbuses + i;
    if mpc.gencost(1,4)==3
      H(Pg_count,Pg_count) = (mpc.gencost(i,5)~=0);
    end
  end
end
  for i=1:nbuses
    H(i,i) = H(i,i) || (mpc.bus(i,5)~=0) || (mpc.bus(i,6)~=0);
  end

  for i=1:nbranches
    if mpc.branch(i,11) == 0; continue; end

    VMf = IDtoCountmap( mpc.branch(i,1) );
    VMt = IDtoCountmap( mpc.branch(i,2) );
    VAf = nbuses + IDtoCountmap( mpc.branch(i,1) );
    VAt = nbuses + IDtoCountmap( mpc.branch(i,2) );

      H(VMf,VMf) = 1;
      H(VMt,VMt) = 1;
      H(VAf,VAf) = 1;
      H(VAt,VAt) = 1;
      H(max(VMf,VMt),min(VMf,VMt)) = 1;
      H(max(VMf,VAf),min(VMf,VAf)) = 1;
      H(max(VMf,VAt),min(VMf,VAt)) = 1;
      H(max(VMt,VAf),min(VMt,VAf)) = 1;
      H(max(VMt,VAt),min(VMt,VAt)) = 1;
      H(max(VAf,VAt),min(VAf,VAt)) = 1;
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

      H(VMf,VMf) = 1;
      H(VMt,VMt) = 1;
      H(VAf,VAf) = 1;
      H(VAt,VAt) = 1;
      H(max(VMf,VMt),min(VMf,VMt)) = 1;
      H(max(VMf,VAf),min(VMf,VAf)) = 1;
      H(max(VMf,VAt),min(VMf,VAt)) = 1;
      H(max(VMt,VAf),min(VMt,VAf)) = 1;
      H(max(VMt,VAt),min(VMt,VAt)) = 1;
      H(max(VAf,VAt),min(VAf,VAt)) = 1;

 

  end


  H = sparse(H);