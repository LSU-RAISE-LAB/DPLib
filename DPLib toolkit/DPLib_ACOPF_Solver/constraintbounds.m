% ----------------------------------------------------------------------
function [cl,cu] = constraintbounds (auxdata)
global regnum loadedData

  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);
  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;

region = ['R', num2str(regnum)];
regvar = loadedData.(['nit', region]);
tie_table_name = ['interregional_tielines', region];
dfjhdf = cell2table(loadedData.(tie_table_name)).Var3;


  cl = zeros(1,2*nbuses+3*nbranches+3*regvar);
  cu = zeros(1,2*nbuses+3*nbranches+3*regvar);




  k = 0;
  cl(k+2*((1:nbuses)-1)+1) = mpc.bus(:,3)/baseMVA;  %Pd
  cu(k+2*((1:nbuses)-1)+1) = mpc.bus(:,3)/baseMVA;  %Pd
  cl(k+2*((1:nbuses)-1)+2) = mpc.bus(:,4)/baseMVA;  %Qd
  cu(k+2*((1:nbuses)-1)+2) = mpc.bus(:,4)/baseMVA;  %Qd

  k = 2*nbuses;
  j = find(mpc.branch(:,11));
  cl(k+2*(j-1)+1) = -Inf;
  cl(k+2*(j-1)+2) = -Inf;
  cu(k+2*(j-1)+1) = Inf;
  cu(k+2*(j-1)+2) = Inf;
  j = find(mpc.branch(:,6)>0);
  cu(k+2*(j-1)+1) = (mpc.branch(j,6)/baseMVA).^2;
  cu(k+2*(j-1)+2) = (mpc.branch(j,6)/baseMVA).^2;

  k = 2*nbuses + 2*nbranches;
  j = find(mpc.branch(:,11)>0);
  cl(k+j) = mpc.anglelim_rad(j,1);
  cu(k+j) = mpc.anglelim_rad(j,2);
  j = find(mpc.branch(:,11)==0);
  cl(k+j) = -Inf;
  cu(k+j) = Inf;

  for i=1:regvar

      var = dfjhdf(i,6);
      if var == 0
        var = 100000;
      end

   cl(2*nbuses+3*nbranches+(i - 1) * 3 +1) =-Inf;
   cu(2*nbuses+3*nbranches+(i - 1) * 3 +1) =(var/baseMVA).^2;
   cl(2*nbuses+3*nbranches+(i - 1) * 3 +2) =-Inf;
   cu(2*nbuses+3*nbranches+(i - 1) * 3 +2) =(var/baseMVA).^2;

    minangle = dfjhdf(i,12);
    maxangle = dfjhdf(i,13);

   cl(2*nbuses+3*nbranches+(i - 1) * 3 +3) =pi*minangle/180;
   cu(2*nbuses+3*nbranches+(i - 1) * 3 +3) =pi*maxangle/180;


  end


  