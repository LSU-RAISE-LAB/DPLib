% -------------------------------------------------------------------------
% Partitioned MATPOWER Case File
% Copyright (c) 2025 Milad Hasanzadeh
% Generated on: 2026-06-14
% This file contains partitioned data for 2 regions


% -------------------- mpc_regionR1 --------------------

mpc_regionR1.version = '2';
mpc_regionR1.baseMVA = 100.0000;
mpc_regionR1.bus = [
%	bus_i       type        Pd          Qd          Gs          Bs          area        Vm          Va          baseKV      zone        Vmax        Vmin        
	 1	 1	 806.16	 98.61	 0	 0	 1	 1	 0	 230	 1	 1.1	 0.9;
	 2	 2	 806.16	 98.61	 0	 0	 1	 1	 0	 230	 1	 1.1	 0.9;
];

mpc_regionR1.gen = [
%	bus         Pg          Qg          Qmax        Qmin        Vg          mBase       status      Pmax        Pmin        Pc1         Pc2         Qc1min      Qc1max      Qc2min      Qc2max      ramp_agc    ramp_10     ramp_30     ramp_q      apf         
	 2	 575	 0	 575	-575	 1	 100	 1	 1150	 0;
];

mpc_regionR1.branch = [
%	fbus        tbus        r           x           b           rateA       rateB       rateC       ratio       angle       status      angmin      angmax      
	 1	 2	 0.00108	 0.0108	 0.01852	 426	 426	 426	 0	 0	 1	-30	 30;
];

mpc_regionR1.gencost = [
%	2           startup     shutdown    n           c(n-1)...c0 
	 2	 0	 0	 3	 0	 30	 0;
];


% -------------------- mpc_regionR2 --------------------

mpc_regionR2.version = '2';
mpc_regionR2.baseMVA = 100.0000;
mpc_regionR2.bus = [
%	bus_i       type        Pd          Qd          Gs          Bs          area        Vm          Va          baseKV      zone        Vmax        Vmin        
	 1	 2	 0	 0	 0	 0	 1	 1	 0	 230	 1	 1.1	 0.9;
	 2	 3	 1074.88	 131.47	 0	 0	 1	 1	 0	 230	 1	 1.1	 0.9;
	 3	 2	 0	 0	 0	 0	 1	 1	 0	 230	 1	 1.1	 0.9;
];

mpc_regionR2.gen = [
%	bus         Pg          Qg          Qmax        Qmin        Vg          mBase       status      Pmax        Pmin        Pc1         Pc2         Qc1min      Qc1max      Qc2min      Qc2max      ramp_agc    ramp_10     ramp_30     ramp_q      apf         
	 1	 111	 0	 177.6	-177.6	 1	 100	 1	 222	 0;
	 1	 170.5	 0	 177.6	-177.6	 1	 100	 1	 341	 0;
	 2	 586	 0	 586	-586	 1	 100	 1	 1172	 0;
	 3	 146	 0	 450	-450	 1	 100	 1	 292	 0;
];

mpc_regionR2.branch = [
%	fbus        tbus        r           x           b           rateA       rateB       rateC       ratio       angle       status      angmin      angmax      
	 1	 2	 0.00304	 0.0304	 0.00658	 426	 426	 426	 0	 0	 1	-30	 30;
	 1	 3	 0.00064	 0.0064	 0.03126	 426	 426	 426	 0	 0	 1	-30	 30;
	 2	 3	 0.00297	 0.0297	 0.00674	 240	 240	 240	 0	 0	 1	-30	 30;
];

mpc_regionR2.gencost = [
%	2           startup     shutdown    n           c(n-1)...c0 
	 2	 0	 0	 3	 0	 14	 0;
	 2	 0	 0	 3	 0	 15	 0;
	 2	 0	 0	 3	 0	 40	 0;
	 2	 0	 0	 3	 0	 10	 0;
];


% ---------- Interregional Tie-Line Information ----------
% Each row: {From_Region (string), To_Region (string), Branch_Data (1x13 array)}
% Branch_Data: {fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax}

interregional_tielines_total = {
	{'mpc_regionR2', 'mpc_regionR1', [ 1	 1	 0.00281	 0.0281	 0.00712	 400	 400	 400	 0	 0	 1	-30	 30]},
	{'mpc_regionR1', 'mpc_regionR2', [ 2	 2	 0.00297	 0.0297	 0.00674	 426	 426	 426	 0	 0	 1	-30	 30]},
};

