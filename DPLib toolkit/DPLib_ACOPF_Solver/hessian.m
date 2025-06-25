
% ----------------------------------------------------------------------
function H = hessian (x, sigma, lambda, auxdata)
  mpc = auxdata{1};
global regnum loadedData
regID = ['R', num2str(regnum)];
regvar = loadedData.(['nit', regID]);
check         = loadedData.(['check', regID]);
region_buses  = loadedData.(['region_', regID]);
G_itl             = loadedData.(['G', regID]);
B_itl             = loadedData.(['B', regID]);
global rho PuScale
  nbuses = size(mpc.bus,1);
  nbranches = size(mpc.branch,1);
  ngens = size(mpc.gen,1);
  IDtoCountmap = mpc.IDtoCountmap;
  baseMVA = mpc.baseMVA;
  G = mpc.G;
  B = mpc.B;

  H = zeros(2*nbuses+2*ngens+2*regvar, 2*nbuses+2*ngens+2*regvar);
if ngens~=0
  for i=1:ngens
    if mpc.gen(i,8)==0; continue; end

    Pg_count = 2*nbuses + i;    
    if mpc.gencost(1,4)==3
      H(Pg_count,Pg_count) = H(Pg_count,Pg_count) + sigma*2*mpc.gencost(i,5)*baseMVA^2;
    end
  end
end

if ngens~=0
VarVal=PuScale*rho*sigma;
else
VarVal=rho*sigma;
end


for i=1:regvar
  H(region_buses{i}+nbuses,region_buses{i}+nbuses) =H(region_buses{i}+nbuses,region_buses{i}+nbuses)+VarVal;
  H(2*nbuses+2*ngens+(2*i),2*nbuses+2*ngens+(2*i)) =H(2*nbuses+2*ngens+(2*i),2*nbuses+2*ngens+(2*i))+VarVal;
  H(region_buses{i} ,region_buses{i}) =H(region_buses{i} ,region_buses{i})+VarVal;
  H(2*nbuses+2*ngens+(2*i-1),2*nbuses+2*ngens+(2*i-1)) =H(2*nbuses+2*ngens+(2*i-1),2*nbuses+2*ngens+(2*i-1))+VarVal;
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
      H(max(VAt,VAf),min(VAt,VAf)) = H(max(VAt,VAf),min(VAt,VAf)) -lambda(k+2*(i-1)+2)*(-B(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) +G(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
    end
  end

  k = 0;
 
for i=1:regvar

      if check(i)==0 %from

            VMf = region_buses{i};
            VMt = 2*nbuses+2*ngens+2*i-1;
            VAf = region_buses{i}+nbuses;
            VAt = 2*nbuses+2*ngens+2*i;j=i;ii=VMf;


      H(VMf,VMf) = H(VMf,VMf) -lambda(k+2*(ii-1)+1)*( 2*G_itl(j,1) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(ii-1)+1)*(-G_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(ii-1)+1)*(-G_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VMt),min(VMf,VMt)) = H(max(VMf,VMt),min(VMf,VMt)) -lambda(k+2*(ii-1)+1)*( G_itl(j,2)*cos(x(VAf)-x(VAt)) + B_itl(j,2)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(ii-1)+1)*(-G_itl(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) + B_itl(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(ii-1)+1)*( G_itl(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) - B_itl(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(ii-1)+1)*(-G_itl(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) + B_itl(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(ii-1)+1)*( G_itl(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) - B_itl(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VAf,VAt),min(VAf,VAt)) = H(max(VAf,VAt),min(VAf,VAt)) -lambda(k+2*(ii-1)+1)*( G_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + B_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );

      H(VMf,VMf) = H(VMf,VMf) -lambda(k+2*(ii-1)+2)*(-2*B_itl(j,1) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(ii-1)+2)*( B_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(ii-1)+2)*( B_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VMt),min(VMf,VMt)) = H(max(VMf,VMt),min(VMf,VMt)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,2)*cos(x(VAf)-x(VAt)) + G_itl(j,2)*sin(x(VAf)-x(VAt)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(ii-1)+2)*( B_itl(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) + G_itl(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,2)*x(VMt)*sin(x(VAf)-x(VAt)) - G_itl(j,2)*x(VMt)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(ii-1)+2)*( B_itl(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) + G_itl(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,2)*x(VMf)*sin(x(VAf)-x(VAt)) - G_itl(j,2)*x(VMf)*cos(x(VAf)-x(VAt)) );
      H(max(VAf,VAt),min(VAf,VAt)) = H(max(VAf,VAt),min(VAf,VAt)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G_itl(j,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) );

      elseif check(i)==1 %to

            VMt = region_buses{i};
            VMf = 2*nbuses+2*ngens+2*i-1;
            VAt = region_buses{i}+nbuses;
            VAf = 2*nbuses+2*ngens+2*i;j=i;ii=VMt;

      H(VMt,VMt) = H(VMt,VMt) -lambda(k+2*(ii-1)+1)*( 2*G_itl(j,4) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(ii-1)+1)*(-G_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(ii-1)+1)*(-G_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VMf),min(VMt,VMf)) = H(max(VMt,VMf),min(VMt,VMf)) -lambda(k+2*(ii-1)+1)*( G_itl(j,3)*cos(x(VAt)-x(VAf)) + B_itl(j,3)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(ii-1)+1)*(-G_itl(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) + B_itl(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(ii-1)+1)*( G_itl(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) - B_itl(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(ii-1)+1)*(-G_itl(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) + B_itl(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(ii-1)+1)*( G_itl(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) - B_itl(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VAt,VAf),min(VAt,VAf)) = H(max(VAt,VAf),min(VAt,VAf)) -lambda(k+2*(ii-1)+1)*( G_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );

      H(VMt,VMt) = H(VMt,VMt) -lambda(k+2*(ii-1)+2)*(-2*B_itl(j,4) );
      H(VAt,VAt) = H(VAt,VAt) -lambda(k+2*(ii-1)+2)*( B_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(VAf,VAf) = H(VAf,VAf) -lambda(k+2*(ii-1)+2)*( B_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VMf),min(VMt,VMf)) = H(max(VMt,VMf),min(VMt,VMf)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,3)*cos(x(VAt)-x(VAf)) + G_itl(j,3)*sin(x(VAt)-x(VAf)) );
      H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) -lambda(k+2*(ii-1)+2)*( B_itl(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) + G_itl(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,3)*x(VMf)*sin(x(VAt)-x(VAf)) - G_itl(j,3)*x(VMf)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) -lambda(k+2*(ii-1)+2)*( B_itl(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) + G_itl(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,3)*x(VMt)*sin(x(VAt)-x(VAf)) - G_itl(j,3)*x(VMt)*cos(x(VAt)-x(VAf)) );
      H(max(VAt,VAf),min(VAt,VAf)) = H(max(VAt,VAf),min(VAt,VAf)) -lambda(k+2*(ii-1)+2)*(-B_itl(j,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G_itl(j,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) );

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
    d2Pfd2VMf =  2*G_itl(i,1);
    d2Pfd2VMt =  0;
    d2Pfd2VAf = -G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2Pfd2VAt = -G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2PfdVMfdVMt =  G_itl(i,2)*cos(x(VAf)-x(VAt)) + B_itl(i,2)*sin(x(VAf)-x(VAt));
    d2PfdVMfdVAf = -G_itl(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) + B_itl(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2PfdVMfdVAt =  G_itl(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) - B_itl(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2PfdVMtdVAf = -G_itl(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) + B_itl(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2PfdVMtdVAt =  G_itl(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) - B_itl(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2PfdVAfdVAt =  G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));

    Qf      = -B_itl(i,1)*x(VMf)^2 - B_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMf = -2*B_itl(i,1)*x(VMf) - B_itl(i,2)*x(VMt)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMt)*sin(x(VAf)-x(VAt));
    dQfdVMt = -B_itl(i,2)*x(VMf)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*sin(x(VAf)-x(VAt));
    dQfdVAf =  B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    dQfdVAt = -B_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt)) - G_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt));
    d2Qfd2VMf = -2*B_itl(i,1);
    d2Qfd2VMt = 0;
    d2Qfd2VAf =  B_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2Qfd2VAt =  B_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) - G_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));
    d2QfdVMfdVMt = -B_itl(i,2)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*sin(x(VAf)-x(VAt));
    d2QfdVMfdVAf =  B_itl(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2QfdVMfdVAt = -B_itl(i,2)*x(VMt)*sin(x(VAf)-x(VAt)) - G_itl(i,2)*x(VMt)*cos(x(VAf)-x(VAt));
    d2QfdVMtdVAf =  B_itl(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2QfdVMtdVAt = -B_itl(i,2)*x(VMf)*sin(x(VAf)-x(VAt)) - G_itl(i,2)*x(VMf)*cos(x(VAf)-x(VAt));
    d2QfdVAfdVAt = -B_itl(i,2)*x(VMf)*x(VMt)*cos(x(VAf)-x(VAt)) + G_itl(i,2)*x(VMf)*x(VMt)*sin(x(VAf)-x(VAt));

    H(VMf,VMf) = H(VMf,VMf) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVMf^2 + 2*Pf*d2Pfd2VMf + 2*dQfdVMf^2 + 2*Qf*d2Qfd2VMf );
    H(VMt,VMt) = H(VMt,VMt) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVMt^2 + 2*Pf*d2Pfd2VMt + 2*dQfdVMt^2 + 2*Qf*d2Qfd2VMt );
    H(VAf,VAf) = H(VAf,VAf) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVAf^2 + 2*Pf*d2Pfd2VAf + 2*dQfdVAf^2 + 2*Qf*d2Qfd2VAf );
    H(VAt,VAt) = H(VAt,VAt) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVAt^2 + 2*Pf*d2Pfd2VAt + 2*dQfdVAt^2 + 2*Qf*d2Qfd2VAt );
    H(max(VMf,VMt),min(VMf,VMt)) = H(max(VMf,VMt),min(VMf,VMt)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVMt*dPfdVMf + 2*Pf*d2PfdVMfdVMt + 2*dQfdVMt*dQfdVMf + 2*Qf*d2QfdVMfdVMt );
    H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVAf*dPfdVMf + 2*Pf*d2PfdVMfdVAf + 2*dQfdVAf*dQfdVMf + 2*Qf*d2QfdVMfdVAf );
    H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVAt*dPfdVMf + 2*Pf*d2PfdVMfdVAt + 2*dQfdVAt*dQfdVMf + 2*Qf*d2QfdVMfdVAt );
    H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVAf*dPfdVMt + 2*Pf*d2PfdVMtdVAf + 2*dQfdVAf*dQfdVMt + 2*Qf*d2QfdVMtdVAf );
    H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVAt*dPfdVMt + 2*Pf*d2PfdVMtdVAt + 2*dQfdVAt*dQfdVMt + 2*Qf*d2QfdVMtdVAt );
    H(max(VAf,VAt),min(VAf,VAt)) = H(max(VAf,VAt),min(VAf,VAt)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +1)*( 2*dPfdVAt*dPfdVAf + 2*Pf*d2PfdVAfdVAt + 2*dQfdVAt*dQfdVAf + 2*Qf*d2QfdVAfdVAt );

    Pt      =  G_itl(i,4)*x(VMt)^2 + G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMt =  2*G_itl(i,4)*x(VMt) + G_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dPtdVMf =  G_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dPtdVAt = -G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dPtdVAf =  G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    d2Ptd2VMt =  2*G_itl(i,4);
    d2Ptd2VMf =  0;
    d2Ptd2VAt = -G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2Ptd2VAf = -G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2PtdVMtdVMf =  G_itl(i,3)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*sin(x(VAt)-x(VAf));
    d2PtdVMtdVAt = -G_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2PtdVMtdVAf =  G_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) - B_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2PtdVMfdVAt = -G_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2PtdVMfdVAf =  G_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) - B_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2PtdVAtdVAf =  G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));

    Qt      = -B_itl(i,4)*x(VMt)^2 - B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMt = -2*B_itl(i,4)*x(VMt) - B_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf));
    dQtdVMf = -B_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf));
    dQtdVAt =  B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    dQtdVAf = -B_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf)) - G_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf));
    d2Qtd2VMt = -2*B_itl(i,4);
    d2Qtd2VMf =  0;
    d2Qtd2VAt =  B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2Qtd2VAf =  B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) - G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));
    d2QtdVMtdVMf = -B_itl(i,3)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*sin(x(VAt)-x(VAf));
    d2QtdVMtdVAt =  B_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2QtdVMtdVAf = -B_itl(i,3)*x(VMf)*sin(x(VAt)-x(VAf)) - G_itl(i,3)*x(VMf)*cos(x(VAt)-x(VAf));
    d2QtdVMfdVAt =  B_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2QtdVMfdVAf = -B_itl(i,3)*x(VMt)*sin(x(VAt)-x(VAf)) - G_itl(i,3)*x(VMt)*cos(x(VAt)-x(VAf));
    d2QtdVAtdVAf = -B_itl(i,3)*x(VMt)*x(VMf)*cos(x(VAt)-x(VAf)) + G_itl(i,3)*x(VMt)*x(VMf)*sin(x(VAt)-x(VAf));

    H(VMt,VMt) = H(VMt,VMt) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVMt^2 + 2*Pt*d2Ptd2VMt + 2*dQtdVMt^2 + 2*Qt*d2Qtd2VMt );
    H(VMf,VMf) = H(VMf,VMf) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVMf^2 + 2*Pt*d2Ptd2VMf + 2*dQtdVMf^2 + 2*Qt*d2Qtd2VMf );
    H(VAt,VAt) = H(VAt,VAt) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVAt^2 + 2*Pt*d2Ptd2VAt + 2*dQtdVAt^2 + 2*Qt*d2Qtd2VAt );
    H(VAf,VAf) = H(VAf,VAf) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVAf^2 + 2*Pt*d2Ptd2VAf + 2*dQtdVAf^2 + 2*Qt*d2Qtd2VAf );
    H(max(VMt,VMf),min(VMt,VMf)) = H(max(VMt,VMf),min(VMt,VMf)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVMf*dPtdVMt + 2*Pt*d2PtdVMtdVMf + 2*dQtdVMf*dQtdVMt + 2*Qt*d2QtdVMtdVMf );
    H(max(VMt,VAt),min(VMt,VAt)) = H(max(VMt,VAt),min(VMt,VAt)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVAt*dPtdVMt + 2*Pt*d2PtdVMtdVAt + 2*dQtdVAt*dQtdVMt + 2*Qt*d2QtdVMtdVAt );
    H(max(VMt,VAf),min(VMt,VAf)) = H(max(VMt,VAf),min(VMt,VAf)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVAf*dPtdVMt + 2*Pt*d2PtdVMtdVAf + 2*dQtdVAf*dQtdVMt + 2*Qt*d2QtdVMtdVAf );
    H(max(VMf,VAt),min(VMf,VAt)) = H(max(VMf,VAt),min(VMf,VAt)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVAt*dPtdVMf + 2*Pt*d2PtdVMfdVAt + 2*dQtdVAt*dQtdVMf + 2*Qt*d2QtdVMfdVAt );
    H(max(VMf,VAf),min(VMf,VAf)) = H(max(VMf,VAf),min(VMf,VAf)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVAf*dPtdVMf + 2*Pt*d2PtdVMfdVAf + 2*dQtdVAf*dQtdVMf + 2*Qt*d2QtdVMfdVAf );
    H(max(VAt,VAf),min(VAt,VAf)) = H(max(VAt,VAf),min(VAt,VAf)) + lambda(2*nbuses+3*nbranches+(i - 1) * 3 +2)*( 2*dPtdVAf*dPtdVAt + 2*Pt*d2PtdVAtdVAf + 2*dQtdVAf*dQtdVAt + 2*Qt*d2QtdVAtdVAf );


  end



  H = sparse(H);
