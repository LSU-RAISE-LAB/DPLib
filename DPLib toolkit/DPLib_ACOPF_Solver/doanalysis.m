

% ----------------------------------------------------------------------
function doanalysis (x, options, info, auxdata, x0_info, t0)
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
  fprintf(fid,'Objective Function Value = %0.2f $/hr\n',objective(x,auxdata));
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
    fprintf(fid,'Converged in %0.4f seconds with status %d (%s)\n',info.cpu,info.status,whichstatus(info.status));
  else
    fprintf(fid,'First fase converged in %0.4f seconds with status %d (%s)\n',x0_info.cpu,x0_info.status,whichstatus(x0_info.status));
    fprintf(fid,'Second fase converged in %0.4f seconds with status %d (%s)\n',info.cpu,info.status,whichstatus(info.status));
  end
  fprintf(fid,'Total time %0.4f seconds \n\n',cputime-t0);
  fprintf(fid,'END\n');
  fclose(fid);
