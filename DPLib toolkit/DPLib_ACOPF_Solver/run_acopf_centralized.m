function [x_opf, obj_value,stat] = run_acopf_centralized(mfilename)
%ACOPF  Runs an AC optimal power flow using IPOPT case files must have Matpower format
%
%  Inputs:
%     MFILENAME: Matpower case file name (with .m extension), it could be just the
%                file name or the file path together with the file name
%
%  Outputs:
%     No outputs, but create two files: One .out containing detailed the solution
%     of the optimal power flow, and other .mat containing the matlab variables of 
%     the solution
%
%  Calling syntax:
%     acopf(mfilename) 
%
% 
%  File created by Mauro Escobar (www.columbia.edu/~me2533/)
%  Based on MATPOWER formats (www.pserc.cornell.edu/matpower/)
%
  
  t0 = cputime;
  mpc = preprocess1(mfilename);
  options.auxdata = {mpc};

  % Set the IPOPT options.
  options.ipopt.print_level              = 0;
  options.ipopt.tol                      = 1e-8;
  options.ipopt.max_iter                 = 250;
  options.ipopt.dual_inf_tol             = 1e-1;
  options.ipopt.compl_inf_tol            = 1e-5;
  options.ipopt.acceptable_tol           = 1e-8;
  options.ipopt.acceptable_compl_inf_tol = 1e-3;
  options.ipopt.mu_strategy              = 'adaptive';
  %options.ipopt.hessian_approximation    = 'limited-memory';
  %options.ipopt.derivative_test          = 'second-order';

  x1 = initialx01 (options.auxdata);
  x0_info = [];
  nbuses = size(mpc.bus,1);
  ngens = size(mpc.gen,1);
  [options.lb, options.ub] = bounds1 (options.auxdata);
  [options.cl, options.cu] = constraintbounds1 (options.auxdata);

  % The callback functions.
  funcs.objective            = @objective1;
  funcs.constraints          = @constraints1;
  funcs.gradient             = @gradient1;
  funcs.jacobian             = @jacobian1;
  funcs.jacobianstructure    = @jacobianstructure1;
  funcs.hessian              = @hessian1;
  funcs.hessianstructure     = @hessianstructure1;

  % Run IPOPT.
  [x_opf, info_opf] = ipopt_auxdata(x1, funcs, options);obj_value = objective1(x_opf, options.auxdata);stat = info_opf.status;



 
% ----------------------------------------------------------------------
function mpc = preprocess1 (filename)
%PREPROCESS1  Load and preprocess a MATPOWER case.
%   Accepts:
%     - 'case9'
%     - 'case9.m'
%     - 'C:\path\case9'
%     - 'C:\path\case9.m'

  % Make sure filename is a char/string and trim spaces
  if ~ischar(filename) && ~isstring(filename)
      error('preprocess1: filename must be a char or string.');
  end
  filename = strtrim(char(filename));

  % Split path, name, extension
  [filepath,name,ext] = fileparts(filename);

  % Handle different input forms:
  %   - no extension      -> assume .m
  %   - .m extension      -> OK
  %   - other extension   -> error
  if isempty(ext)
      ext = '.m';
  elseif ~strcmpi(ext, '.m')
      error('Case filename must be a MATPOWER .m file (got "%s").', ext);
  end

  % Rebuild full filename (for existence check, if desired)
  fullfile_name = fullfile(filepath, [name ext]);

  % Optional: check that the file actually exists
  if ~isempty(filepath) && ~exist(fullfile_name, 'file')
      error('Case file "%s" not found.', fullfile_name);
  end

  % Change directory (if path is provided), call case function, then return
  currFld = pwd;
  if ~isempty(filepath)
      cd(filepath);
  end
  % feval uses only the function name (no .m and no path)
  mpc = feval(name);
  cd(currFld);

  % ---- rest of your original code below ----

  if max(mpc.gencost(:,1))>2 || min(mpc.gencost(:,1))<2
    error('Only generation cost functions of type 2 are supported')
  end

  if max(mpc.gencost(:,4))>min(mpc.gencost(:,4))
    error('Generation cost functions of type 2 are supportedmust be of the same degree (2 or 3) for all generators')
  end

  mpc.casename = name;
  mpc.bus = mpc.bus(mpc.bus(:,2)<4,:);
  mpc.branch = mpc.branch(mpc.branch(:,11)==1,:);
  baseMVA = mpc.baseMVA;
  nbuses = size(mpc.bus,1);
  ngens = size(mpc.gen,1);
  nbranches = size(mpc.branch,1);

  IDtoCountmap = zeros(1,nbuses);
  genids = cell(1,nbuses);
  frombranchids = cell(1,nbuses);
  tobranchids = cell(1,nbuses);
  for i=1:nbuses
    IDtoCountmap( mpc.bus(i,1) ) = i;
    if mpc.bus(i,2)==3
      mpc.refbuscount = i;
    end
  end

  for i=1:ngens
    if mpc.gen(i,8)==0; continue; end
    genids{IDtoCountmap(mpc.gen(i,1))} = [genids{IDtoCountmap(mpc.gen(i,1))} i];
  end

  G = zeros(nbranches,4);
  B = zeros(nbranches,4);
  Y = zeros(nbranches,4);
  Z = zeros(nbranches,1);
  anglelim_rad = zeros(nbranches,2);
  for i=1:nbranches
    if mpc.branch(i,11)==0; continue; end   % status==0
    r = mpc.branch(i,3);
    x = mpc.branch(i,4);
    bc = mpc.branch(i,5);
    ratio = mpc.branch(i,9);
    if ratio==0
      ratio = 1;
    end
    angle = mpc.branch(i,10);
    angle_rad = pi*angle/180;

    invratio2 = 1/ratio^2;
    multtf = 1/(ratio*exp(1j*angle_rad));
    multft = 1/(ratio*exp(-1j*angle_rad));
    z = r + 1j*x;  Z(i) = z;
    y = 1/z;

    Yff = (y + bc/2*1j)*invratio2;
    Yft = -y*multft;
    Ytf = -y*multtf;
    Ytt = y + bc/2*1j;

    Y(i,1) = Yff;  Y(i,2) = Yft;  Y(i,3) = Ytf;  Y(i,4) = Ytt;

    G(i,1) = real(Yff);  B(i,1) = imag(Yff);
    G(i,2) = real(Yft);  B(i,2) = imag(Yft);
    G(i,3) = real(Ytf);  B(i,3) = imag(Ytf);
    G(i,4) = real(Ytt);  B(i,4) = imag(Ytt);

    minangle = mpc.branch(i,12);
    maxangle = mpc.branch(i,13);
    if (minangle==0 && maxangle==0) || (minangle==-360 && maxangle==360)
      minangle = -180; maxangle = 180;
    end
    anglelim_rad(i,1) = pi*minangle/180;
    anglelim_rad(i,2) = pi*maxangle/180;

    frombranchids{IDtoCountmap(mpc.branch(i,1))} = [frombranchids{IDtoCountmap(mpc.branch(i,1))} i];
    tobranchids{IDtoCountmap(mpc.branch(i,2))} = [tobranchids{IDtoCountmap(mpc.branch(i,2))} i];
  end

  mpc.G = G;
  mpc.B = B;
  mpc.Y = Y;
  mpc.Z = Z;
  mpc.anglelim_rad = anglelim_rad;
  mpc.genids = genids;
  mpc.frombranchids = frombranchids;
  mpc.tobranchids = tobranchids;
  mpc.IDtoCountmap = IDtoCountmap;



% ----------------------------------------------------------------------
function x0 = initialx01 (auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  x0 = zeros(2*nbuses+2*ngens,1);
  x0(1:nbuses) = 1;

% ----------------------------------------------------------------------
function [lb,ub] = bounds1 (auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;

  lb = zeros(2*nbuses+2*ngens,1);
  ub = zeros(2*nbuses+2*ngens,1);

  lb(1:nbuses) = mpc.bus(:,13);
  ub(1:nbuses) = mpc.bus(:,12);

  k = nbuses;
  lb(k+(1:nbuses)) = -pi;
  ub(k+(1:nbuses)) = pi;
  lb(k+mpc.refbuscount) = pi*mpc.bus(mpc.refbuscount,9)/180;
  ub(k+mpc.refbuscount) = pi*mpc.bus(mpc.refbuscount,9)/180;

  k = 2*nbuses;
  lb(k+(1:ngens)) = mpc.gen(:,10).*mpc.gen(:,8)/baseMVA;
  ub(k+(1:ngens)) = mpc.gen(:,9).*mpc.gen(:,8)/baseMVA;

  k = 2*nbuses + ngens;
  lb(k+(1:ngens)) = mpc.gen(:,5).*mpc.gen(:,8)/baseMVA;
  ub(k+(1:ngens)) = mpc.gen(:,4).*mpc.gen(:,8)/baseMVA;

% ----------------------------------------------------------------------
function f = objective1 (x, auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  ngens = size(mpc.gen,1);
  baseMVA = mpc.baseMVA;

  Pg = 2*nbuses + (1:ngens);
  actgen = mpc.gen(:,8);  
  if mpc.gencost(1,4)==3
    f = sum(actgen.*( baseMVA^2*mpc.gencost(:,5).*x(Pg).^2 + baseMVA*mpc.gencost(:,6).*x(Pg) + mpc.gencost(:,7) ));
  elseif mpc.gencost(1,4)==2
    f = sum(actgen.*( baseMVA*mpc.gencost(:,5).*x(Pg) + mpc.gencost(:,6) ));
  end
  
% ----------------------------------------------------------------------
function g = gradient1 (x, auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  ngens = size(mpc.gen,1);
  baseMVA = mpc.baseMVA;

  Pg = 2*nbuses + (1:ngens);
  actgen = mpc.gen(:,8);
  g = zeros(2*nbuses+2*ngens,1);
  if mpc.gencost(1,4)==3
    g(Pg) = actgen.*( 2*baseMVA^2*mpc.gencost(:,5).*x(Pg) + baseMVA*mpc.gencost(:,6) );
  elseif mpc.gencost(1,4)==2
    g(Pg) = actgen.*( baseMVA*mpc.gencost(:,5) );
  end

% ----------------------------------------------------------------------
function c = constraints1 (x, auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;
  G = mpc.G;
  B = mpc.B;

  c = zeros(2*nbuses + 3*nbranches,1);

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

% ----------------------------------------------------------------------
function [cl,cu] = constraintbounds1 (auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;

  cl = zeros(1,2*nbuses+3*nbranches);
  cu = zeros(1,2*nbuses+3*nbranches);

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

% ----------------------------------------------------------------------
function J = jacobian1 (x, auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;
  G = mpc.G;
  B = mpc.B;

  J = zeros(2*nbuses+3*nbranches, 2*nbuses+2*ngens);

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

  k = 2*nbuses + 2*nbranches;
  for i=1:nbranches
    if mpc.branch(i,11) == 0; continue; end

    VAf_count = nbuses + mpc.IDtoCountmap(mpc.branch(i,1));
    VAt_count = nbuses + mpc.IDtoCountmap(mpc.branch(i,2));
    J(k+i,VAf_count) = J(k+i,VAf_count) + 1;
    J(k+i,VAt_count) = J(k+i,VAt_count) - 1;
  end

  J = sparse(J);

% ----------------------------------------------------------------------
function J = jacobianstructure1 (auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;

  J = zeros(2*nbuses+3*nbranches, 2*nbuses+2*ngens);

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

  J = sparse(J);

% ----------------------------------------------------------------------
function H = hessian1 (x, sigma, lambda, auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;
  G = mpc.G;
  B = mpc.B;

  H = zeros(2*nbuses+2*ngens, 2*nbuses+2*ngens);

  for i=1:ngens
    if mpc.gen(i,8)==0; continue; end

    Pg_count = 2*nbuses + i;    
    if mpc.gencost(1,4)==3
      H(Pg_count,Pg_count) = H(Pg_count,Pg_count) + sigma*2*mpc.gencost(i,5)*baseMVA^2;
    end
  end

  k = 0;
  for i=1:nbuses
    H(i,i) = H(i,i) - lambda(k+2*(i-1)+1)*2*mpc.bus(i,5)/baseMVA;
    H(i,i) = H(i,i) + lambda(k+2*(i-1)+2)*2*mpc.bus(i,6)/baseMVA;

    for j=mpc.frombranchids{i}
      if mpc.branch(j,11) == 0; continue; end

      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      H(VMf,VMf) = H(VMf,VMf) -lambda(k+2*(i-1)+1)*( 2*G(j,1) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(i-1)+1)*(-G(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(i-1)+1)*(-G(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VMt),min(VMf,VMt)) = H(max(VMf,VMt),min(VMf,VMt)) -lambda(k+2*(i-1)+1)*( G(j,2)*cos(x(VAf)-x(VAt)) + B(j,2)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(i-1)+1)*(-G(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) + B(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(i-1)+1)*( G(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) - B(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(i-1)+1)*(-G(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) + B(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(i-1)+1)*( G(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) - B(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VAf,VAt),min(VAf,VAt)) = H(max(VAf,VAt),min(VAf,VAt)) -lambda(k+2*(i-1)+1)*( G(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + B(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );

      H(VMf,VMf) = H(VMf,VMf) -lambda(k+2*(i-1)+2)*(-2*B(j,1) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(i-1)+2)*( B(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(i-1)+2)*( B(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VMt),min(VMf,VMt)) = H(max(VMf,VMt),min(VMf,VMt)) -lambda(k+2*(i-1)+2)*(-B(j,2)*cos(x(VAf)-x(VAt)) + G(j,2)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(i-1)+2)*( B(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) + G(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(i-1)+2)*(-B(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) - G(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(i-1)+2)*( B(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) + G(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(i-1)+2)*(-B(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) - G(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VAf,VAt),min(VAf,VAt)) = H(max(VAf,VAt),min(VAf,VAt)) -lambda(k+2*(i-1)+2)*(-B(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
    end

    for j=mpc.tobranchids{i}
      if mpc.branch(j,11) == 0; continue; end

      VMf = IDtoCountmap( mpc.branch(j,1) );
      VMt = IDtoCountmap( mpc.branch(j,2) );
      VAf = nbuses + IDtoCountmap( mpc.branch(j,1) );
      VAt = nbuses + IDtoCountmap( mpc.branch(j,2) );

      H(VMt,VMt) = H(VMt,VMt) -lambda(k+2*(i-1)+1)*( 2*G(j,4) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(i-1)+1)*(-G(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(i-1)+1)*(-G(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VMf),min(VMt,VMf)) = H(max(VMt,VMf),min(VMt,VMf)) -lambda(k+2*(i-1)+1)*( G(j,3)*cos(x(VAt)-x(VAf)) + B(j,3)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(i-1)+1)*(-G(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) + B(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(i-1)+1)*( G(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) - B(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(i-1)+1)*(-G(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) + B(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(i-1)+1)*( G(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) - B(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VAt,VAf),min(VAt,VAf)) = H(max(VAt,VAf),min(VAt,VAf)) -lambda(k+2*(i-1)+1)*( G(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );

      H(VMt,VMt) = H(VMt,VMt) -lambda(k+2*(i-1)+2)*(-2*B(j,4) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(i-1)+2)*( B(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(i-1)+2)*( B(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VMf),min(VMt,VMf)) = H(max(VMt,VMf),min(VMt,VMf)) -lambda(k+2*(i-1)+2)*(-B(j,3)*cos(x(VAt)-x(VAf)) + G(j,3)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(i-1)+2)*( B(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) + G(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(i-1)+2)*(-B(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) - G(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(i-1)+2)*( B(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) + G(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(i-1)+2)*(-B(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) - G(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VAt,VAf),min(VAt,VAf)) = H(max(VAt,VAf),min(VAt,VAf)) -lambda(k+2*(i-1)+2)*(-B(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
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
    d2Pfd2VMf =  2*G(i,1);
    d2Pfd2VMt =  0;
    d2Pfd2VAf = -G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2Pfd2VAt = -G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2PfdVMfdVMt =  G(i,2)*cos(x(VAf)-x(VAt)) + B(i,2)*sin(x(VAf)-x(VAt));
    d2PfdVMfdVAf = -G(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) + B(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2PfdVMfdVAt =  G(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) - B(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2PfdVMtdVAf = -G(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) + B(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2PfdVMtdVAt =  G(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) - B(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2PfdVAfdVAt =  G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));

    Qf      = -B(i,1)*x(VMf)^2 - B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMf = -2*B(i,1)*x(VMf) - B(i,2)*x(VMt)*cos(x(VAf)-x(VAt)) + G(i,2)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMt = -B(i,2)*x(VMf)*cos(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*sin(x(VAf)-x(VAt));
    dQfdVAf =  B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    dQfdVAt = -B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    d2Qfd2VMf = -2*B(i,1);
    d2Qfd2VMt = 0;
    d2Qfd2VAf =  B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2Qfd2VAt =  B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2QfdVMfdVMt = -B(i,2)*cos(x(VAf)-x(VAt)) + G(i,2)*sin(x(VAf)-x(VAt));
    d2QfdVMfdVAf =  B(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) + G(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2QfdVMfdVAt = -B(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) - G(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2QfdVMtdVAf =  B(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2QfdVMtdVAt = -B(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) - G(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2QfdVAfdVAt = -B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));

    H(VMf,VMf) = H(VMf,VMf) + lambda(k+2*(i-1)+1)*( 2*dPfdVMf^2 + 2*Pf*d2Pfd2VMf + 2*dQfdVMf^2 + 2*Qf*d2Qfd2VMf );
    H(VMt,VMt) = H(VMt,VMt) + lambda(k+2*(i-1)+1)*( 2*dPfdVMt^2 + 2*Pf*d2Pfd2VMt + 2*dQfdVMt^2 + 2*Qf*d2Qfd2VMt );
    H(VAf,VAf) = H(VAf,VAf) + lambda(k+2*(i-1)+1)*( 2*dPfdVAf^2 + 2*Pf*d2Pfd2VAf + 2*dQfdVAf^2 + 2*Qf*d2Qfd2VAf );
    H(VAt,VAt) = H(VAt,VAt) + lambda(k+2*(i-1)+1)*( 2*dPfdVAt^2 + 2*Pf*d2Pfd2VAt + 2*dQfdVAt^2 + 2*Qf*d2Qfd2VAt );
    H(max(VMf,VMt),min(VMf,VMt)) = H(max(VMf,VMt),min(VMf,VMt)) + lambda(k+2*(i-1)+1)*( 2*dPfdVMt*dPfdVMf + 2*Pf*d2PfdVMfdVMt + 2*dQfdVMt*dQfdVMf + 2*Qf*d2QfdVMfdVMt );
    H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) + lambda(k+2*(i-1)+1)*( 2*dPfdVAf*dPfdVMf + 2*Pf*d2PfdVMfdVAf + 2*dQfdVAf*dQfdVMf + 2*Qf*d2QfdVMfdVAf );
    H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) + lambda(k+2*(i-1)+1)*( 2*dPfdVAt*dPfdVMf + 2*Pf*d2PfdVMfdVAt + 2*dQfdVAt*dQfdVMf + 2*Qf*d2QfdVMfdVAt );
    H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) + lambda(k+2*(i-1)+1)*( 2*dPfdVAf*dPfdVMt + 2*Pf*d2PfdVMtdVAf + 2*dQfdVAf*dQfdVMt + 2*Qf*d2QfdVMtdVAf );
    H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) + lambda(k+2*(i-1)+1)*( 2*dPfdVAt*dPfdVMt + 2*Pf*d2PfdVMtdVAt + 2*dQfdVAt*dQfdVMt + 2*Qf*d2QfdVMtdVAt );
    H(max(VAf,VAt),min(VAf,VAt)) = H(max(VAf,VAt),min(VAf,VAt)) + lambda(k+2*(i-1)+1)*( 2*dPfdVAt*dPfdVAf + 2*Pf*d2PfdVAfdVAt + 2*dQfdVAt*dQfdVAf + 2*Qf*d2QfdVAfdVAt );

    Pt      =  G(i,4)*x(VMt)^2 + G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMt =  2*G(i,4)*x(VMt) + G(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + B(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMf =  G(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dPtdVAt = -G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dPtdVAf =  G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    d2Ptd2VMt =  2*G(i,4);
    d2Ptd2VMf =  0;
    d2Ptd2VAt = -G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2Ptd2VAf = -G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2PtdVMtdVMf =  G(i,3)*cos(x(VAt)-x(VAf)) + B(i,3)*sin(x(VAt)-x(VAf));
    d2PtdVMtdVAt = -G(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) + B(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2PtdVMtdVAf =  G(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) - B(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2PtdVMfdVAt = -G(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2PtdVMfdVAf =  G(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) - B(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2PtdVAtdVAf =  G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));

    Qt      = -B(i,4)*x(VMt)^2 - B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMt = -2*B(i,4)*x(VMt) - B(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + G(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMf = -B(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dQtdVAt =  B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dQtdVAf = -B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    d2Qtd2VMt = -2*B(i,4);
    d2Qtd2VMf =  0;
    d2Qtd2VAt =  B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2Qtd2VAf =  B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2QtdVMtdVMf = -B(i,3)*cos(x(VAt)-x(VAf)) + G(i,3)*sin(x(VAt)-x(VAf));
    d2QtdVMtdVAt =  B(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) + G(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2QtdVMtdVAf = -B(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) - G(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2QtdVMfdVAt =  B(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2QtdVMfdVAf = -B(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) - G(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2QtdVAtdVAf = -B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));

    H(VMt,VMt) = H(VMt,VMt) + lambda(k+2*(i-1)+2)*( 2*dPtdVMt^2 + 2*Pt*d2Ptd2VMt + 2*dQtdVMt^2 + 2*Qt*d2Qtd2VMt );
    H(VMf,VMf) = H(VMf,VMf) + lambda(k+2*(i-1)+2)*( 2*dPtdVMf^2 + 2*Pt*d2Ptd2VMf + 2*dQtdVMf^2 + 2*Qt*d2Qtd2VMf );
    H(VAt,VAt) = H(VAt,VAt) + lambda(k+2*(i-1)+2)*( 2*dPtdVAt^2 + 2*Pt*d2Ptd2VAt + 2*dQtdVAt^2 + 2*Qt*d2Qtd2VAt );
    H(VAf,VAf) = H(VAf,VAf) + lambda(k+2*(i-1)+2)*( 2*dPtdVAf^2 + 2*Pt*d2Ptd2VAf + 2*dQtdVAf^2 + 2*Qt*d2Qtd2VAf );
    H(max(VMt,VMf),min(VMt,VMf)) = H(max(VMt,VMf),min(VMt,VMf)) + lambda(k+2*(i-1)+2)*( 2*dPtdVMf*dPtdVMt + 2*Pt*d2PtdVMtdVMf + 2*dQtdVMf*dQtdVMt + 2*Qt*d2QtdVMtdVMf );
    H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) + lambda(k+2*(i-1)+2)*( 2*dPtdVAt*dPtdVMt + 2*Pt*d2PtdVMtdVAt + 2*dQtdVAt*dQtdVMt + 2*Qt*d2QtdVMtdVAt );
    H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) + lambda(k+2*(i-1)+2)*( 2*dPtdVAf*dPtdVMt + 2*Pt*d2PtdVMtdVAf + 2*dQtdVAf*dQtdVMt + 2*Qt*d2QtdVMtdVAf );
    H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) + lambda(k+2*(i-1)+2)*( 2*dPtdVAt*dPtdVMf + 2*Pt*d2PtdVMfdVAt + 2*dQtdVAt*dQtdVMf + 2*Qt*d2QtdVMfdVAt );
    H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) + lambda(k+2*(i-1)+2)*( 2*dPtdVAf*dPtdVMf + 2*Pt*d2PtdVMfdVAf + 2*dQtdVAf*dQtdVMf + 2*Qt*d2QtdVMfdVAf );
    H(max(VAt,VAf),min(VAt,VAf)) = H(max(VAt,VAf),min(VAt,VAf)) + lambda(k+2*(i-1)+2)*( 2*dPtdVAf*dPtdVAt + 2*Pt*d2PtdVAtdVAf + 2*dQtdVAf*dQtdVAt + 2*Qt*d2QtdVAtdVAf );
  end

  H = sparse(H);

% ----------------------------------------------------------------------
function H = hessianstructure1 (auxdata)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;

  H = zeros(2*nbuses+2*ngens, 2*nbuses+2*ngens);

  for i=1:ngens
    if mpc.gen(i,8)==0; continue; end

    Pg_count = 2*nbuses + i;
    if mpc.gencost(1,4)==3
      H(Pg_count,Pg_count) = (mpc.gencost(i,5)~=0);
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

  H = sparse(H);

% ----------------------------------------------------------------------
function doanalysis1 (x, options, info, auxdata, x0_info, t0)
  mpc = auxdata{1};
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);

  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;
  G = mpc.G;
  B = mpc.B;

  fid = fopen(['acopf_solution_' mpc.casename '.out'],'w');

  k = 2*nbuses;
  lossQbr = zeros(1,nbranches);
  lossPbr = zeros(1,nbranches);
  fchg = zeros(1,nbranches);
  tchg = zeros(1,nbranches);
  for i=1:nbranches
    if mpc.branch(i,11) == 0; continue; end

    VMf = IDtoCountmap( mpc.branch(i,1) );
    VMt = IDtoCountmap( mpc.branch(i,2) );
    VAf = nbuses + IDtoCountmap( mpc.branch(i,1) );
    VAt = nbuses + IDtoCountmap( mpc.branch(i,2) );

    Pf = 2*nbuses + 2*ngens + i;
    x(Pf) = G(i,1)*x(VMf)^2 +G(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) +B(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));

    Qf = 2*nbuses + 2*ngens + nbranches + i;
    x(Qf) =-B(i,1)*x(VMf)^2 -B(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) +G(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));

    Pt = 2*nbuses + 2*ngens + 2*nbranches + i;
    x(Pt) = G(i,4)*x(VMt)^2 +G(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) +B(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));

    Qt = 2*nbuses + 2*ngens + 3*nbranches + i;
    x(Qt) =-B(i,4)*x(VMt)^2 -B(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) +G(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    Ybr = [mpc.Y(i,1), mpc.Y(i,2); mpc.Y(i,3), mpc.Y(i,4)];
    Vbr = [ x(IDtoCountmap(mpc.branch(i,1)))*exp(1j*x(nbuses+IDtoCountmap(mpc.branch(i,1)))); ...
            x(IDtoCountmap(mpc.branch(i,2)))*exp(1j*x(nbuses+IDtoCountmap(mpc.branch(i,2)))) ];
    Ibr = Ybr*Vbr;
    lossbranch = (-Vbr(1)*Ybr(1,2)*mpc.Z(i)-Vbr(2))*conj(-Vbr(1)*Ybr(1,2)*mpc.Z(i)-Vbr(2))/conj(mpc.Z(i));
    lossQbr(i) = baseMVA*imag(lossbranch);
    lossPbr(i) = baseMVA*real(lossbranch);
    fchg(i) = baseMVA*(-Vbr(1)*Ybr(1,2)*mpc.Z(i))*conj(-Vbr(1)*Ybr(1,2)*mpc.Z(i))*mpc.branch(i,5)/2;
    tchg(i) = baseMVA*Vbr(2)*conj(Vbr(2))*mpc.branch(i,5)/2;
  end

  sumPd = 0; sumQd = 0;
  PG = zeros(1,nbuses); QG = zeros(1,nbuses);
  for i=1:nbuses
    sumPd = sumPd + mpc.bus(i,3);
    sumQd = sumQd + mpc.bus(i,4);

    for genid=mpc.genids{i}
      Pg_count = 2*nbuses + genid;
      Qg_count = 2*nbuses + ngens + genid;
      PG(i) = PG(i) + baseMVA*x(Pg_count);
      QG(i) = QG(i) + baseMVA*x(Qg_count);
    end
  end


  fprintf(fid,'\n');
  fprintf(fid,'Objective Function Value = %0.2f $/hr\n',objective1(x,auxdata));
  fprintf(fid,'================================================================================\n');
  fprintf(fid,'|     System Summary                                                           |\n');
  fprintf(fid,'================================================================================\n');
  fprintf(fid,'\n');
  fprintf(fid,'How many?                How much?              P (MW)            Q (MVAr)\n');
  fprintf(fid,'---------------------    -------------------  -------------  -----------------\n');
  fprintf(fid,'Buses          %5d',nbuses);
  fprintf(fid,'     Total Gen Capacity %9.1f',sum(mpc.gen(:,9)));
  fprintf(fid,'        %0.1f to %0.1f\n',sum(mpc.gen(:,5)),sum(mpc.gen(:,4)));
  fprintf(fid,'Generators     %5d',ngens);
  fprintf(fid,'     On-line Capacity   %9.1f',sum((mpc.gen(:,8)==1).*mpc.gen(:,9)));
  fprintf(fid,'        %0.1f to %0.1f\n',sum((mpc.gen(:,8)==1).*mpc.gen(:,5)),sum((mpc.gen(:,8)==1).*mpc.gen(:,4)));
  fprintf(fid,'Committed Gens %5d     Generation (actual)%9.1f         %9.1f\n',sum(mpc.gen(:,8)==1),sum(PG),sum(QG));
  fprintf(fid,'Loads          %5d     Load               %9.1f         %9.1f\n',sum(sum(abs(mpc.bus(:,[3 4])'))~=0),sumPd,sumQd);
  fprintf(fid,'  Fixed        %5d       Fixed            %9.1f         %9.1f\n',sum(sum(abs(mpc.bus(:,[3 4])'))~=0),sumPd,sumQd);
  fprintf(fid,'  Dispatchable             Dispatchable     \n');
  fprintf(fid,'Shunts         %5d', sum(sum(abs(mpc.bus(:,[5 6])'))~=0) );
  fprintf(fid,'     Shunt (inj)        %9.1f         %9.1f\n',-sum(mpc.bus(:,5).*x(1:nbuses).*x(1:nbuses)),sum(mpc.bus(:,6).*x(1:nbuses).*x(1:nbuses)));
  fprintf(fid,'Branches       %5d     Losses (I^2 * Z)    %9.2f         %9.2f\n',nbranches,sum(lossPbr),sum(lossQbr));
  fprintf(fid,'Transformers   %5d     Branch Charging (inj)     -          %9.1f\n',sum(mpc.branch(:,9)~=0),sum(fchg)+sum(tchg));
  fprintf(fid,'Inter-ties               Total Inter-tie Flow\n');
  fprintf(fid,'Areas          %5d\n',max(mpc.bus(:,7))+(max(mpc.bus(:,7)))==0);
  fprintf(fid,'\n');
  fprintf(fid,'                          Minimum                      Maximum\n');
  fprintf(fid,'                 -------------------------  --------------------------------\n');
  fprintf(fid,'Voltage Magnitude  %6.3f p.u. @ bus %-5d',min(x(1:nbuses)),mpc.bus(find(x(1:nbuses)==min(x(1:nbuses)),1),1));
  fprintf(fid,'     %6.3f p.u. @ bus %-5d\n',max(x(1:nbuses)),mpc.bus(find(x(1:nbuses)==max(x(1:nbuses)),1),1));
  fprintf(fid,'Voltage Angle    %7.2f deg   @ bus %-5d',180/pi*min(x((nbuses+1):2*nbuses)),mpc.bus(find(x((nbuses+1):2*nbuses)==min(x((nbuses+1):2*nbuses)),1),1));
  fprintf(fid,'   %7.2f deg   @ bus %-5d\n',180/pi*max(x((nbuses+1):2*nbuses)),mpc.bus(find(x((nbuses+1):2*nbuses)==max(x((nbuses+1):2*nbuses)),1),1));
  fprintf(fid,'P Losses (I^2*R)             -             %9.2f MW    @ line %d-%d\n',max(lossPbr),mpc.branch(find(lossPbr==max(lossPbr),1),1),mpc.branch(find(lossPbr==max(lossPbr),1),2));
  fprintf(fid,'Q Losses (I^2*X)             -             %9.2f MVar  @ line %d-%d\n',max(lossQbr),mpc.branch(find(lossQbr==max(lossQbr),1),1),mpc.branch(find(lossQbr==max(lossQbr),1),2));
  fprintf(fid,'Lambda P         %7.2f $/MWh @ bus %-5d',min(-info.lambda(1:2:2*nbuses)/baseMVA),mpc.bus(find(info.lambda(1:2:2*nbuses)==max(info.lambda(1:2:2*nbuses)),1),1));
  fprintf(fid,'   %7.2f $/MWh @ bus %-5d\n',max(-info.lambda(1:2:2*nbuses)/baseMVA),mpc.bus(find(info.lambda(1:2:2*nbuses)==min(info.lambda(1:2:2*nbuses)),1),1));
  fprintf(fid,'Lambda Q         %7.2f $/MWh @ bus %-5d',min(-info.lambda(2:2:2*nbuses)/baseMVA),mpc.bus(find(info.lambda(2:2:2*nbuses)==max(info.lambda(2:2:2*nbuses)),1),1));
  fprintf(fid,'   %7.2f $/MWh @ bus %-5d\n',max(-info.lambda(2:2:2*nbuses)/baseMVA),mpc.bus(find(info.lambda(2:2:2*nbuses)==min(info.lambda(2:2:2*nbuses)),1),1));
  fprintf(fid,'\n');
  fprintf(fid,'==========================================================================================================\n');
  fprintf(fid,'|     Bus Data    (*: reference bus, fixed voltage phase angle)                                          |\n');
  fprintf(fid,'==========================================================================================================\n');
  fprintf(fid,' Bus  ------- Voltage --------------------------      Generation             Load         Lambda($/MVA-hr)\n');
  fprintf(fid,'  #   Mag(pu) Ang(deg) Ang(rad)   Real    Imag     P (MW)   Q (MVAr)   P (MW)   Q (MVAr)     P        Q   \n');
  fprintf(fid,'----- ------- -------- --------  ------  -------  --------  --------  --------  --------  -------  -------\n');
  for i=1:nbuses
    if mpc.bus(i,2)==3
      fprintf(fid,'%4d*  %5.3f  %7.3f  %7.4f  %7.4f  %7.4f',mpc.bus(i,1),x(i),x(nbuses+i)*180/pi,x(nbuses+i),real(x(i)*exp(1j*x(nbuses+i))),imag(x(i)*exp(1j*x(nbuses+i))));
    else
      fprintf(fid,'%4d   %5.3f  %7.3f  %7.4f  %7.4f  %7.4f',mpc.bus(i,1),x(i),x(nbuses+i)*180/pi,x(nbuses+i),real(x(i)*exp(1j*x(nbuses+i))),imag(x(i)*exp(1j*x(nbuses+i))));
    end
    if length(mpc.genids{i})>0
      fprintf(fid,'  %8.2f  %8.2f',PG(i),QG(i));
    else
      fprintf(fid,'       -         -  ');
    end
    if mpc.bus(i,3)~=0 || mpc.bus(i,4)~=0
      fprintf(fid,'  %8.2f  %8.2f  ',mpc.bus(i,3),mpc.bus(i,4));
    else
      fprintf(fid,'       -         -    ');
    end
    if abs(info.lambda(2*i-1)/baseMVA)<0.0005
      fprintf(fid,'     -  ');
    else
      fprintf(fid,' %7.3f',-info.lambda(2*i-1)/baseMVA);
    end
    if abs(info.lambda(2*i)/baseMVA)<0.0005
      fprintf(fid,'     -  ');
    else
      fprintf(fid,' %7.3f',-info.lambda(2*i)/baseMVA);
    end
    fprintf(fid,'\n');
  end
  fprintf(fid,'                                                  --------  --------  --------  --------\n');
  fprintf(fid,'                                          Total:  %8.2f  %8.2f',sum(PG),sum(QG));
  fprintf(fid,'  %8.2f  %8.2f\n',sumPd,sumQd);
  fprintf(fid,'\n');
  fprintf(fid,'\n');
  fprintf(fid,'==========================================================================================================\n');
  fprintf(fid,'|     Branch Data                                                                                        |\n');
  fprintf(fid,'==========================================================================================================\n');
  fprintf(fid,'Brnch   From   To           -----From Bus Injection-----  ------To Bus Injection------    Loss (I^2 * Z)  \n');
  fprintf(fid,'  #     Bus    Bus    Lim    P (MW)   Q (MVAr)     S       P (MW)   Q (MVAr)     S       P (MW)   Q (MVAr)\n');
  fprintf(fid,'-----  -----  -----  -----  --------  --------  --------  --------  --------  --------  --------  --------\n');
  for i=1:nbranches
    Pfcount = 2*nbuses+2*ngens+i;
  Qfcount = 2*nbuses+2*ngens+nbranches+i;	    Ptcount = 2*nbuses+2*ngens+2*nbranches+i;
 Qtcount = 2*nbuses+2*ngens+3*nbranches+i;
    fprintf(fid,'%4d %6d %6d %7.0f',i,mpc.branch(i,1),mpc.branch(i,2),baseMVA*sqrt(options.cu(2*nbuses+2*i)));
    fprintf(fid,' %8.2f  %8.2f  %8.2f ',baseMVA*x(Pfcount),baseMVA*x(Qfcount),baseMVA*sqrt(x(Pfcount)^2+x(Qfcount)^2));
    fprintf(fid,' %8.2f  %8.2f  %8.2f ',baseMVA*x(Ptcount),baseMVA*x(Qtcount),baseMVA*sqrt(x(Ptcount)^2+x(Qtcount)^2));
    fprintf(fid,' %8.3f  %8.2f\n',lossPbr(i),lossQbr(i));

  end
  fprintf(fid,'                                                                                        --------  --------\n');
  fprintf(fid,'                                                                               Total:  %8.3f  %8.2f\n',sum(lossPbr),sum(lossQbr));
  if sum(info.zu(1:nbuses)>1e-3)+sum(info.zl(1:nbuses)>1e-3)>0
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'================================================================================\n');
    fprintf(fid,'|     Voltage Constraints                                                      |\n');
    fprintf(fid,'================================================================================\n');
    fprintf(fid,'Bus #  Vmin mu    Vmin    |V|   Vmax    Vmax mu\n');
    fprintf(fid,'-----  --------   -----  -----  -----   --------\n');
    for i = 1:nbuses
      if info.zu(i)>1e-3
        fprintf(fid,'%5d      -      %5.3f  %5.3f  %5.3f %10.3f\n',mpc.bus(i,1),options.lb(i),x(i),options.ub(i),info.zu(i));
      elseif info.zl(i)>1e-3
        fprintf(fid,'%5d %10.3f  %5.3f  %5.3f  %5.3f      -       \n',mpc.bus(i,1),info.zl(i),options.lb(i),x(i),options.ub(i));
      end
    end
  end

  if sum(info.zl(2*nbuses+1:2*nbuses+ngens)>1e-3)+sum(info.zu(2*nbuses+1:2*nbuses+ngens)>1e-3)+sum(info.zl(2*nbuses+ngens+1:2*nbuses+2*ngens)>1e-3)+sum(info.zu(2*nbuses+ngens+1:2*nbuses+2*ngens)>1e-3)>0
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'================================================================================\n');
    fprintf(fid,'|     Generation Constraints                                                   |\n');
    fprintf(fid,'================================================================================\n');
    if sum(info.zl(2*nbuses+1:2*nbuses+ngens)>1e-3)+sum(info.zu(2*nbuses+1:2*nbuses+ngens)>1e-3)>0
      fprintf(fid,' Gen   Bus                  Active Power Limits\n');
      fprintf(fid,'  #     #     Pmin mu     Pmin       Pg       Pmax    Pmax mu\n');
      fprintf(fid,'----  -----   -------   --------  --------  --------  -------\n');
      for i = 2*nbuses+1:2*nbuses+ngens
        if info.zl(i)>1e-3
          fprintf(fid,'%4d %5d  %9.3f %9.2f %9.2f %9.2f      -\n',i-2*nbuses,mpc.gen(i-2*nbuses,1),info.zl(i)/baseMVA,baseMVA*options.lb(i),baseMVA*x(i),baseMVA*options.ub(i));
        elseif info.zu(i)>1e-3
          fprintf(fid,'%4d %5d       -    %9.2f %9.2f %9.2f %9.3f\n',i-2*nbuses,mpc.gen(i-2*nbuses,1),baseMVA*options.lb(i),baseMVA*x(i),baseMVA*options.ub(i),info.zu(i)/baseMVA);
        end
      end
    fprintf(fid,'\n');
    end
    if sum(info.zl(2*nbuses+ngens+1:2*nbuses+2*ngens)>1e-3)+sum(info.zu(2*nbuses+ngens+1:2*nbuses+2*ngens)>1e-3)>0
      fprintf(fid,' Gen   Bus                 Reactive Power Limits\n');
      fprintf(fid,'  #     #     Qmin mu     Qmin       Qg       Qmax    Qmax mu\n');
      fprintf(fid,'----  -----   -------   --------  --------  --------  -------\n');
      for i = 2*nbuses+ngens+1:2*nbuses+2*ngens
        if info.zl(i)>1e-3
          fprintf(fid,'%4d %5d  %9.3f %9.2f %9.2f %9.2f      -\n',i-(2*nbuses+ngens),mpc.gen(i-(2*nbuses+ngens),1),info.zl(i)/baseMVA,baseMVA*options.lb(i),baseMVA*x(i),baseMVA*options.ub(i));
        elseif info.zu(i)>1e-3
          fprintf(fid,'%4d %5d       -    %9.2f %9.2f %9.2f %9.3f\n',i-(2*nbuses+ngens),mpc.gen(i-(2*nbuses+ngens),1),baseMVA*options.lb(i),baseMVA*x(i),baseMVA*options.ub(i),info.zu(i)/baseMVA);
        end
      end
      fprintf(fid,'\n');
    end
  end

  if sum(info.lambda((2*nbuses+1):(2*nbuses+2*nbranches))>1e-3)>0
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'================================================================================\n');
    fprintf(fid,'|     Branch Flow Constraints                                                  |\n');
    fprintf(fid,'================================================================================\n');
    fprintf(fid,'Brnch   From     "From" End        Limit       "To" End        To\n');
    fprintf(fid,'  #     Bus   |Sf| mu    |Sf|     |Smax|     |St|    |St| mu   Bus\n');
    fprintf(fid,'-----  -----  -------  --------  --------  --------  -------  -----\n');
    for i = (2*nbuses+1):2:(2*nbuses+2*nbranches)
      if info.lambda(i)>1e-3 || info.lambda(i+1)>1e-3
        brcount = (i+1-(2*nbuses))/2;
        factor = (1/baseMVA)^2 + 2*sqrt(options.cu(i))/baseMVA;
        fprintf(fid,'%4d  %5d  ',brcount,mpc.branch(brcount,1));
        if info.lambda(i)>1e-3
          fprintf(fid,'%7.3f ',info.lambda(i)*factor);
        else
          fprintf(fid,'   -    ');
        end
        fprintf(fid,'%9.2f ',baseMVA*sqrt(x(2*nbuses+2*ngens+brcount)^2+x(2*nbuses+2*ngens+nbranches+brcount)^2));
        fprintf(fid,'%9.2f ',baseMVA*sqrt(options.cu(i)));
        fprintf(fid,'%9.2f ',baseMVA*sqrt(x(2*nbuses+2*ngens+2*nbranches+brcount)^2+x(2*nbuses+2*ngens+3*nbranches+brcount)^2));
        if info.lambda(i+1)>1e-3
          fprintf(fid,' %7.3f ',info.lambda(i+1)*factor);
        else
          fprintf(fid,'    -    ');
        end
        fprintf(fid,' %5d  ',mpc.branch(brcount,2));
        fprintf(fid,'\n');
      end
    end
  end

  fprintf(fid,'\n');
  fprintf(fid,'\n');
  fprintf(fid,'\n');
  fprintf(fid,'\n');
  fprintf(fid,'Objective Function Value = %0.2f $/hr\n',objective(x,auxdata));
  if isempty(x0_info)
    fprintf(fid,'Converged in %0.4f seconds with status %d (%s)\n',info.cpu,info.status,whichstatus1(info.status));
  else
    fprintf(fid,'First fase converged in %0.4f seconds with status %d (%s)\n',x0_info.cpu,x0_info.status,whichstatus1(x0_info.status));
    fprintf(fid,'Second fase converged in %0.4f seconds with status %d (%s)\n',info.cpu,info.status,whichstatus1(info.status));
  end
  fprintf(fid,'Total time %0.4f seconds \n\n',cputime-t0);
  fprintf(fid,'END\n');
  fclose(fid);

% ----------------------------------------------------------------------
function r = whichstatus1 (n)
  r = '';
  if n==0
    r = 'solved';
  elseif n==1
    r = 'solved to acceptable level';
  elseif n==2
    r = 'infeasible problem detected';
  elseif n==3
    r = 'search direction becomes too small';
  elseif n==4
    r = 'diverging iterates';
  elseif n==5
    r = 'user requested stop';
  elseif n==-1
    r =  'maximum number of iterations exceeded';
  elseif n==-2
    r = 'restoration phase failed';
  elseif n==-3
    r = 'error in step computation';
  elseif n==-10
    r = 'not enough degrees of freedom';
  elseif n==-11
    r = 'invalid problem definition';
  elseif n==-12
    r = 'invalid option';
  elseif n==-13
    r = 'invalid number detected';
  elseif n==-100
    r = 'unrecoverable exception';
  elseif n==-101
    r = 'non-IPOPT exception thrown';
  elseif n==-102
    r = 'insufficient memory';
  elseif n==-199
    r = 'internal error';
  end

