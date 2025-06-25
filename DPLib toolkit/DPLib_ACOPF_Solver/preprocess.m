% ----------------------------------------------------------------------
function mpc = preprocess (filename)
global regnum
reg = ['mpcR' num2str(regnum)];
 % Check if the input is a structure and has required fields
if ~isstruct(filename)
    error('Input must be a MATPOWER case struct in the workspace.');
end

% Assign the input structure to mpc
mpc = filename;

% Optional: Validate required MATPOWER fields
requiredFields = {'bus', 'gen', 'branch', 'gencost', 'baseMVA'};
for k = 1:length(requiredFields)
    if ~isfield(mpc, requiredFields{k})
        error('Field "%s" is missing from the MATPOWER case structure.', requiredFields{k});
    end
end

% Validate generator cost types
if size(mpc.gencost, 1) ~= 0
    if max(mpc.gencost(:, 1)) > 2 || min(mpc.gencost(:, 1)) < 2
        error('Only generation cost functions of type 2 are supported.');
    end

    if max(mpc.gencost(:, 4)) > min(mpc.gencost(:, 4))
        error(['Generation cost functions of type 2 must be of the same degree (2 or 3) ' ...
               'for all generators.']);
    end
end


  mpc.casename = reg;
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

    G(i,1) = real(Yff);	 B(i,1) = imag(Yff);
    G(i,2) = real(Yft);	 B(i,2) = imag(Yft);
    G(i,3) = real(Ytf);	 B(i,3) = imag(Ytf);
    G(i,4) = real(Ytt);	 B(i,4) = imag(Ytt);

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



  