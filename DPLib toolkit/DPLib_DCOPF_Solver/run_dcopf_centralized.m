function [res_dc, fval] = run_dcopf_centralized(case_name)
%RUN_DCOPF_CENTRALIZED  Centralized DC-OPF using YALMIP + MATPOWER (Pg in p.u.)
%
%   [res_dc, fval] = run_dcopf_centralized(case_name)
%
%   DC model (p.u.): P_ij = (theta_i - theta_j) / x_ij
yalmip('clear');
%% Load MATPOWER case
mpc = loadcase(case_name);

baseMVA = mpc.baseMVA;
nb = size(mpc.bus, 1);
ng = size(mpc.gen, 1);
nl = size(mpc.branch, 1);

% Bus IDs (may not be 1:nb, so we map IDs -> indices)
bus_ids     = mpc.bus(:, 1);           % e.g., [101; 102; ...]
max_bus_id  = max(bus_ids);
bus_map     = zeros(max_bus_id, 1);    % bus_map(bus_id) = row index in bus matrix

for idx = 1:nb
    bus_map(bus_ids(idx)) = idx;
end

%% Decision variables (Pg in p.u., Theta in rad)
Pg    = sdpvar(ng, 1);     % generator active powers (p.u.)
Theta = sdpvar(nb, 1);     % bus voltage angles (rad), indexed by bus ROW

%% Objective function (quadratic cost in $ using Pg in MW)
gencost   = mpc.gencost;
gen_coeff = gencost(:, [4 5 6 7]);  % [NCOST, c2, c1, c0]

Objective = 0;
for g = 1:ng
    n_poly = gen_coeff(g, 1);
    a      = gen_coeff(g, 2);   % c2
    b      = gen_coeff(g, 3);   % c1
    c      = gen_coeff(g, 4);   % c0

    Pg_MW = Pg(g) * baseMVA;

    if n_poly == 3
        Objective = Objective + a*Pg_MW^2 + b*Pg_MW + c;
    elseif n_poly == 2
        Objective = Objective + a*Pg_MW + b;
    else
        % if NCOST > 3, you can ignore extra terms or extend logic
        Objective = Objective + a*Pg_MW^2 + b*Pg_MW + c;
    end
end


%% Constraints
Constraints = [];

% Reference bus (slack)
ref_idx = find(mpc.bus(:, 2) == 3, 1);   % row index of REF bus
if isempty(ref_idx)
    ref_idx = 1;
end
Constraints = [Constraints, Theta(ref_idx) == 0];

% Generator limits (Pg in p.u.)
gen_status = mpc.gen(:, 8);     % 1 = online, 0 = offline

Pmax_MW = mpc.gen(:, 9);
Pmin_MW = mpc.gen(:, 10);

Pmax_MW(gen_status == 0) = 0;
Pmin_MW(gen_status == 0) = 0;

Pmax_pu = Pmax_MW / baseMVA;
Pmin_pu = Pmin_MW / baseMVA;

Constraints = [Constraints, Pmin_pu <= Pg <= Pmax_pu];
Constraints = [Constraints, -pi <= Theta <= pi];

%% Bus loads in p.u.
Pd_MW = mpc.bus(:, 3);
Pd_pu = Pd_MW / baseMVA;

%% Power balance at each bus in p.u.
for i = 1:nb
    bus_id_i = bus_ids(i);

    % Generators at this bus (p.u.)
    gen_at_bus = find(mpc.gen(:, 1) == bus_id_i);
    Pg_i_pu = sum(Pg(gen_at_bus));

    Pd_i_pu = Pd_pu(i);

    flow_sum_pu = 0;

    % From-bus lines
    from_lines = find(mpc.branch(:, 1) == bus_id_i);
    for k = from_lines'
        if mpc.branch(k, 11) == 0, continue; end

        x = mpc.branch(k, 4);      % p.u.

        j_bus_id = mpc.branch(k, 2);
        j = bus_map(j_bus_id);

        flow_sum_pu = flow_sum_pu + (Theta(i) - Theta(j)) / x;
    end

    % To-bus lines
    to_lines = find(mpc.branch(:, 2) == bus_id_i);
    for k = to_lines'
        if mpc.branch(k, 11) == 0, continue; end

        x = mpc.branch(k, 4);

        j_bus_id = mpc.branch(k, 1);
        j = bus_map(j_bus_id);

        flow_sum_pu = flow_sum_pu + (Theta(i) - Theta(j)) / x;
    end

    Constraints = [Constraints, Pg_i_pu - Pd_i_pu == flow_sum_pu];
end

%% Line flow limits in p.u.
Pflows = sdpvar(nl, 1);   % p.u.

for k = 1:nl
    if mpc.branch(k, 11) == 0
        continue;
    end

    from_bus_id = mpc.branch(k, 1);
    to_bus_id   = mpc.branch(k, 2);

    x = mpc.branch(k, 4);     % p.u.
    rateA_MW = mpc.branch(k, 6);

    i = bus_map(from_bus_id);
    j = bus_map(to_bus_id);

    Pij_pu = (Theta(i) - Theta(j)) / x;
    Pflows(k) = Pij_pu;

    if rateA_MW > 0
        rateA_pu = rateA_MW / baseMVA;
        Constraints = [Constraints, -rateA_pu <= Pij_pu <= rateA_pu];
    end
end

%% Solve the optimization problem
options = sdpsettings('solver', 'quadprog', 'verbose', 0);
sol = optimize(Constraints, Objective, options);

%% Build result structure (output MW for Pg, Pflows)
res_dc = struct();
res_dc.status = sol.problem;

if sol.problem == 0
    Pg_pu     = value(Pg);
    Theta_val = value(Theta);
    Pflows_pu = value(Pflows);

    res_dc.Pg     = Pg_pu * baseMVA;      % MW
    res_dc.Theta  = Theta_val;            % rad
    res_dc.Pflows = Pflows_pu * baseMVA;  % MW
    res_dc.cost   = value(Objective);
    fval          = res_dc.cost;
else
    res_dc.Pg     = [];
    res_dc.Theta  = [];
    res_dc.Pflows = [];
    res_dc.cost   = Inf;
    fval          = Inf;
end

end
