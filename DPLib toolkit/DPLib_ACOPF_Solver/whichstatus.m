
% ----------------------------------------------------------------------
function r = whichstatus (n)
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
