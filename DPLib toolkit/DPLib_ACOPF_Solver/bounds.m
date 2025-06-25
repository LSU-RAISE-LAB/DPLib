
% ----------------------------------------------------------------------
function [lb,ub] = bounds (auxdata)
global regnum loadedData
regID = ['R', num2str(regnum)];
regvar = loadedData.(['nit', regID]);
  mpc = auxdata{1};bmpc = mpc.bus;
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;

  lb = zeros(2*nbuses+2*ngens+2*regvar,1);
  ub = zeros(2*nbuses+2*ngens+2*regvar,1);

  for i=1:regvar 
   lb(2*nbuses+2*ngens+(2*i)-1) =0.9;
   ub(2*nbuses+2*ngens+(2*i)-1) = 1.1;
   lb(2*nbuses+2*ngens+(2*i)) = -pi;
   ub(2*nbuses+2*ngens+(2*i)) = pi;
  end

  

  lb(1:nbuses) = mpc.bus(:,13);
  ub(1:nbuses) = mpc.bus(:,12);

  k = nbuses;
  lb(k+(1:nbuses)) = -pi;
  ub(k+(1:nbuses)) = pi;

      flag_enable_refbus = 0;
    row_refbus = 0;

    % Loop through each row of matrix A
    for i = 1:size(bmpc, 1)
        if bmpc(i, 2) == 3
            % If a '3' is found in the second column
            flag_enable_refbus = 1;
            row_refbus = i;
        end
    end

    if flag_enable_refbus==1
  lb(k+row_refbus) = 0;
  ub(k+row_refbus) = 0;
    end


  k = 2*nbuses;
  if ngens~=0
  lb(k+(1:ngens)) = mpc.gen(:,10).*mpc.gen(:,8)/baseMVA;
  ub(k+(1:ngens)) = mpc.gen(:,9).*mpc.gen(:,8)/baseMVA;

  k = 2*nbuses + ngens;
  lb(k+(1:ngens)) = mpc.gen(:,5).*mpc.gen(:,8)/baseMVA;
  ub(k+(1:ngens)) = mpc.gen(:,4).*mpc.gen(:,8)/baseMVA;

  end