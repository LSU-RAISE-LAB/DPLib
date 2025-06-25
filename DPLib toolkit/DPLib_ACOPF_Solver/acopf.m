% clear all
% clc
function [out,co,co2] = acopf(reg)


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
%save_as_mpc('mpc_regionA.mat', 'mpcA');
filenamesp = ['mpc', reg];           % 'mpcR1', 'mpcR2', etc.
mpc_struct = evalin('base', filenamesp);  % gets the struct from workspace
mpc = preprocess(mpc_struct);             % pass actual structure to the function


  t0 = cputime;co2 = 1;
  ngens = size(mpc.gen,1);
  options.auxdata = {mpc};
  if ngens~=0
  co = mpc.gencost(:, [5 6 7]);
  else
      co =[];
co2 = 0;
  end
  % Set the IPOPT options.
  options.ipopt.print_level              = 0;
  options.ipopt.tol                      = 1e-6;
  options.ipopt.max_iter                 = 500;
  options.ipopt.dual_inf_tol             = 1e-1;
  options.ipopt.compl_inf_tol            = 1e-5;
  options.ipopt.acceptable_tol           = 1e-8;
  options.ipopt.acceptable_compl_inf_tol = 1e-3;
  options.ipopt.mu_strategy              = 'adaptive';
    options.ipopt.nlp_scaling_method = 'none'; % No scaling
  %options.ipopt.hessian_approximation    = 'limited-memory';
  %options.ipopt.derivative_test          = 'second-order';

  x1 = initialx0 (options.auxdata);
  x0_info = [];
  nbuses = size(mpc.bus,1);
  ngens = size(mpc.gen,1);
  [options.lb, options.ub] = bounds (options.auxdata);
  [options.cl, options.cu] = constraintbounds (options.auxdata);

  % The callback functions.
  funcs.objective            = @objective;
  funcs.constraints          = @constraints;
  funcs.gradient             = @gradient;
  funcs.jacobian             = @jacobian;
  funcs.jacobianstructure    = @jacobianstructure;
  funcs.hessian              = @hessian;
  funcs.hessianstructure     = @hessianstructure;

  % Run IPOPT.
  [out, info_opf] = ipopt_auxdata(x1, funcs, options);


 


