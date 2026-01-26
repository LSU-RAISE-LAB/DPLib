% ========================================================================
% Local DCOPF subproblem for region "reg" (YALMIP version)
% Compatible with DistDCOPF.m
%
% Assumes the following are available in the caller workspace (DistDCOPF):
%   - loadedData       : struct with fields:
%       .(region_var)  : mpc_regionR#   (MATPOWER struct for region)
%       .filename      : original full-case filename (string)
%       .(tie_var)     : interregional_tielinesR# (cell, per-region tielines)
%   - reg              : e.g. 'R1', 'R2', ...
%   - lambda           : [2*nit x 1] dual variables (global)
%   - rho              : ADMM penalty parameter
%   - nit              : total number of global interfaces
%   - nbR#, nitR#, conR#, YonR#, region_R# in base workspace
%
% On exit, leaves:
%   - Pg               : sdpvar(#gen,1)
%   - Delta            : sdpvar(nb_local + nit_local,1)
%   - gen_coeff        : cost coeffs used by DistDCOPF for reporting
% ========================================================================
yalmip('clear');
%% --- Identify region and load region-specific data ----------------------
region_var = ['mpc_region' reg];
tie_var    = ['interregional_tielines' reg];

mpc  = loadedData.(region_var);  % regional MATPOWER case
bmpc = mpc.bus;
br_status = mpc.branch(:, 11);   % 1 = in-service, 0 = out-of-service

x = mpc.branch(:, 4);
invX = 1 ./ x;


% for RTE cases, slightly adjust line limits later
case_name = loadedData.filename;
hasRTE    = contains(case_name, 'rte', 'IgnoreCase', true);

% region-specific tie-line descriptor cell (as built in DistDCOPF)
dataMatrix = loadedData.(tie_var);

% Sort tie-line rows so that region1_num < region2_num (consistent with DistDCOPF)
for i = 1:size(dataMatrix, 1)
    region1_num = sscanf(char(dataMatrix{i, 1}), 'mpc_regionR%d');
    region2_num = sscanf(char(dataMatrix{i, 2}), 'mpc_regionR%d');
    if region1_num > region2_num
        [dataMatrix{i, 1}, dataMatrix{i, 2}] = deal(dataMatrix{i, 2}, dataMatrix{i, 1});
    end
end
Tablek = cell2table(dataMatrix);

%% --- Load region-specific ADMM data from base workspace -----------------
region_busesq = evalin('base', ['region_' reg]);   % cell: local bus index for each tie-line
nbq           = evalin('base', ['nb' reg]);        % # local buses
nitq          = evalin('base', ['nit' reg]);       % # local interfaces (tie-lines in this region)
conq          = evalin('base', ['con' reg]);       % mapping to neighbor local index
Yonq          = evalin('base', ['Yon' reg]);       % global interface index for each local interface

%% --- Define decision variables -----------------------------------------
Pg    = sdpvar(size(mpc.gen, 1), 1);          % generator active power [p.u.]
Delta = sdpvar(size(mpc.bus, 1) + nitq, 1);   % bus angles [rad], local + remote copies

%% --- Local generation cost objective -----------------------------------
% MATPOWER poly cost: [ model startup shutdown n c_{n-1} c_{n-2} ... c0 ]
% Here we assume model=2, n up to 3, so we use columns [4 5 6 7] => [n, c2, c1, c0]
gen_coeff = mpc.gencost(:, [4 5 6 7]);
baseMW    = mpc.baseMVA;PScale = 1;

PgMW = Pg * baseMW;   % convert all generators to MW in one shot

n_poly = gen_coeff(:, 1);
a      = gen_coeff(:, 2);
b      = gen_coeff(:, 3);
c      = gen_coeff(:, 4);

quad_idx = (n_poly == 3);
lin_idx  = (n_poly == 2);

Objective = 0;

% Quadratic costs
if any(quad_idx)
    Objective = Objective + sum( ...
        a(quad_idx) .* (PgMW(quad_idx).^2) + ...
        b(quad_idx) .* PgMW(quad_idx)      + ...
        c(quad_idx) );
end

% Linear costs
if any(lin_idx)
    Objective = Objective + sum( ...
        a(lin_idx) .* PgMW(lin_idx) + ...
        b(lin_idx) );
end


%% --- ADMM consistency residuals for tie-line angles --------------------
val_num = sscanf(reg, 'R%d');   % numeric region index for this region

Rgrad      = [];  % primal consistency residuals
Rgrad_dual = [];  % "dual" residuals (second angle pair per interface)

for i = 1:nitq
    region1 = sscanf(char(Tablek{i, 1}), 'mpc_regionR%d');
    region2 = sscanf(char(Tablek{i, 2}), 'mpc_regionR%d');

    % Identify neighbor region (the "other" one)
    if val_num == region1
        neighbor_reg = ['R' num2str(region2)];
    else
        neighbor_reg = ['R' num2str(region1)];
    end

    % Fetch neighbor phase angles and nit_neighbor from base
    phase_angles_neighbor = evalin('base', ['phase_angles_' neighbor_reg]);
    nit_neighbor          = evalin('base', ['nit' neighbor_reg]);

    % conq(i) = local interface index in neighbor region (as built in DistDCOPF)
    neigh_local_idx = conq(i);

    % angle_2 = neighbor local bus angle, angle_1 = neighbor "remote copy" angle
    angle_2 = phase_angles_neighbor(neigh_local_idx, 1);
    angle_1 = phase_angles_neighbor(nit_neighbor + neigh_local_idx, 1);

    % local bus and remote-bus indices in this region
    local_bus  = region_busesq{i};    % bus index in [1..nbq]
    remote_bus = nbq + i;             % dummy bus index for remote side

    % Define residuals so that global ADMM updates stay consistent
    if val_num == region1
        % This region is "region1" in the (region1, region2) ordering
        %   r(i)        = Delta(local)     - angle_remote(neighbor)
        %   r_dual(i)   = Delta(remote)    - angle_local(neighbor)
        Rgrad      = [Rgrad;      Delta(local_bus)  - angle_1];
        Rgrad_dual = [Rgrad_dual; Delta(remote_bus) - angle_2];
    else
        % This region is "region2"
        % Choose complementary signs so that global sum of residuals is well-defined
        %   r(i)        = angle_local(neighbor)  - Delta(remote)
        %   r_dual(i)   = angle_remote(neighbor) - Delta(local)
        Rgrad      = [Rgrad;      angle_2 - Delta(remote_bus)];
        Rgrad_dual = [Rgrad_dual; angle_1 - Delta(local_bus)];
    end
end

%% --- Add ADMM penalty terms to objective (augmented Lagrangian) --------
% Standard form per interface k:
%   L_aug += lambda_k * r_k + (rho/2)*r_k^2
%          + lambda_{k+nit} * r_dual_k + (rho/2)*r_dual_k^2

if ~isempty(Yonq)
    for k = 1:numel(Yonq)
        y = Yonq(k);   % global interface index for this local interface

        % ---- Normalize local residuals by the same theta_scale ----
        r_loc      = Rgrad(k)      / theta_scale;
        r_dual_loc = Rgrad_dual(k) / theta_scale;

        Objective = Objective + PScale*(lambda(y)       * r_loc      + (rho / 2) * r_loc^2);
        Objective = Objective + PScale*(lambda(y + nit) * r_dual_loc + (rho / 2) * r_dual_loc^2);
    end
end


%% --- Constraints -------------------------------------------------------
Constraints = [];

% Angle bounds for all local + remote buses
Constraints = [Constraints, -pi <= Delta <= pi];

% Reference bus (if exists in this region)
flag_enable_refbus = 0;
row_refbus         = 0;
for i = 1:size(bmpc, 1)
    if bmpc(i, 2) == 3    % bus type 3 = reference
        flag_enable_refbus = 1;
        row_refbus         = i;
        break;
    end
end
if flag_enable_refbus == 1
    Constraints = [Constraints, Delta(row_refbus) == 0];
end

% Generator limits (Pg in p.u.)
gen_status = mpc.gen(:, 8);     % 1 = online, 0 = offline

Pg_min = mpc.gen(:, 10)/baseMW;
Pg_max = mpc.gen(:,  9)/baseMW;

Pg_max(gen_status == 0) = 0;
Pg_min(gen_status == 0) = 0;

Constraints = [Constraints, Pg_min <= Pg <= Pg_max];

% Nodal power balance for local buses (DC power flow)
for i = 1:size(mpc.bus, 1)
    % Generation at this bus (Note: here bus index == bus column 1 in this region)
    power_injections = sum(Pg(mpc.gen(:, 1) == i));  % p.u.

    % Demand at this bus (MW -> p.u.)
    power_demands = mpc.bus(i, 3) / baseMW;

    % Branches incident to this bus, only in-service ones
    branch_out = find(mpc.branch(:, 1) == i & br_status == 1);
    branch_in  = find(mpc.branch(:, 2) == i & br_status == 1);

    % DC flows: (θ_from - θ_to)/x, index-based
    power_flows_out = sum( (Delta(i) - Delta(mpc.branch(branch_out, 2))) ./ mpc.branch(branch_out, 4) );
    power_flows_in  = sum( (Delta(mpc.branch(branch_in, 1)) - Delta(i)) ./ mpc.branch(branch_in, 4) );

    Constraints = [Constraints, power_injections + power_flows_in == power_demands + power_flows_out];
end

for k = 1:size(mpc.branch, 1)
    if ~br_status(k)
        continue;
    end

    from  = mpc.branch(k, 1);
    to    = mpc.branch(k, 2);
    rateA = mpc.branch(k, 6);

    if rateA > 0
        wq = rateA / baseMW;
        flow_k = (Delta(from) - Delta(to)) * invX(k);
        Constraints = [Constraints, -wq <= flow_k <= wq];
    end
end


%% --- Solve local optimization ------------------------------------------
options = sdpsettings('solver', 'quadprog', 'verbose', 0);
sol     = optimize(Constraints, Objective, options);
