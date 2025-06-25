function x0 = initialx0 (auxdata)
global regnum loadedData
regID = ['R', num2str(regnum)];
regvar = loadedData.(['nit', regID]);

  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  x0 = zeros(2*nbuses+2*ngens+2*regvar,1);
  x0(1:nbuses) = 1;

  for i=1:regvar
   x0(2*nbuses+2*ngens+(2*i)-1)=1;
  end

  