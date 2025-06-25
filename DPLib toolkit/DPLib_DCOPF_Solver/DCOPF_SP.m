
% === Identify region and load region-specific data ===
region_var = ['mpc_region' reg];
tie_var = ['interregional_tielines' reg];


mpc = loadedData.(region_var);
bmpc = mpc.bus;numbus=loadedData.filename;hasRTE = contains(numbus, 'rte', 'IgnoreCase', true);
dataMatrix = loadedData.(tie_var);

% === Sort tie-line rows based on region indices ===
for i = 1:size(dataMatrix, 1)
    region1_num = sscanf(char(dataMatrix{i, 1}), 'mpc_regionR%d');
    region2_num = sscanf(char(dataMatrix{i, 2}), 'mpc_regionR%d');
    if region1_num > region2_num
        [dataMatrix{i, 1}, dataMatrix{i, 2}] = deal(dataMatrix{i, 2}, dataMatrix{i, 1});
    end
end
Tablek = cell2table(dataMatrix);

% === Load region-specific parameters from base workspace ===
region_busesq = evalin('base', ['region_' reg]);
nbq = evalin('base', ['nb' reg]);
nitq = evalin('base', ['nit' reg]);
conq = evalin('base', ['con' reg]);
Yonq = evalin('base', ['Yon' reg]);

% === Define symbolic decision variables ===
Pg = sdpvar(size(mpc.gen, 1), 1);
Delta = sdpvar(size(mpc.bus, 1) + nitq, 1);

% === Local objective function for generation cost ===
gen_coeff = mpc.gencost(:, [5 6 7]);
Objective = sum(gen_coeff(:,1) .* Pg.^2 + gen_coeff(:,2) .* Pg + gen_coeff(:,3));


val_num = sscanf(reg, 'R%d');

% === Symbolic gradient computation before solving ===
for i = 1:nitq
    region1 = sscanf(char(Tablek{i, 1}), 'mpc_regionR%d');
    region2 = sscanf(char(Tablek{i, 2}), 'mpc_regionR%d');

    % Identify neighbor region
    if val_num == region1
        neighbor_reg = ['R' num2str(region2)];
    elseif val_num == region2
        neighbor_reg = ['R' num2str(region1)];
    end

    % Fetch neighbor phase angles and tie-bus index offset
    phase_angles_neighbor = evalin('base', ['phase_angles_' neighbor_reg]);
    nit_neighbor = evalin('base', ['nit' neighbor_reg]);

    % Safely extract constant values from neighbor phase angles
    angle_1 = phase_angles_neighbor(nit_neighbor + conq(i), 1);
    angle_2 = phase_angles_neighbor(conq(i), 1);

    if any(isnan([angle_1, angle_2]))
        error('NaN in phase angle data at tie-line %d', i);
    end

    % Compute symbolic gradients
    if val_num == region1
        Rgrad(i) = Delta(region_busesq{i}) - angle_1;
        Rgrad_dual(i) = Delta(nbq + i) - angle_2;
    else
        Rgrad(i) = angle_2 - Delta(nbq + i);
        Rgrad_dual(i) = angle_1 - Delta(region_busesq{i});
    end
end

% === Add consistency penalties to the objective ===
for k = 1:numel(Yonq)
    y = Yonq(k);
    Objective = Objective + PuScale * (lambda(abs(y),1) * Rgrad(k) + (rho / 2) * Rgrad(k)^2);
    Objective = Objective + PuScale * (lambda(abs(y) + nit, 1) * Rgrad_dual(k) + (rho / 2) * Rgrad_dual(k)^2);
end
                        

%Constraints for areas
Constraints = [];

% constraints on delta
Constraints = [Constraints,-pi <= Delta <= pi];

    flag_enable_refbus = 0;
    row_refbus = 0;

    % Loop through each row of matrix A
    for i = 1:size(bmpc, 1)
        if bmpc(i, 2) == 3
            % If a '3' is found in the second column
            flag_enable_refbus = 1;
            row_refbus = i;
        end
    end

    if flag_enable_refbus==1
Constraints = [Constraints,Delta(row_refbus) == 0];
    end
% Correcting the generator constraints for Area 1
Constraints = [Constraints, mpc.gen(:,10)/mpc.baseMVA <= Pg <= mpc.gen(:,9)/mpc.baseMVA];


%Power balance constraints for area 1

    for i = 1:size(mpc.bus, 1)
        power_injections = sum(Pg(mpc.gen(:, 1) == i));
        power_demands = mpc.bus(i, 3)/mpc.baseMVA;
        branch_indices_out = find(mpc.branch(:, 1) == i);
        branch_indices_in = find(mpc.branch(:, 2) == i);
        power_flows_out = sum((Delta(i) - Delta(mpc.branch(branch_indices_out, 2))) ./ mpc.branch(branch_indices_out, 4));
        power_flows_in = sum((Delta(mpc.branch(branch_indices_in, 1)) - Delta(i)) ./ mpc.branch(branch_indices_in, 4));
        Constraints = [Constraints, power_injections + power_flows_in == power_demands + power_flows_out];
    end



%Line flow limit of area 1
for k = 1:size(mpc.branch, 1)
    from = mpc.branch(k, 1);
    to = mpc.branch(k, 2);
    x = mpc.branch(k, 4);
    line_capacity = mpc.branch(k, 6);
    wq=line_capacity/mpc.baseMVA;
    if wq==0
        wq=Inf;
    end

   if hasRTE==1
        if wq<10
        wq=20;
        end
   end

    Constraints = [Constraints,-wq <= (Delta(from) - Delta(to)) / x <= wq];
end
options = sdpsettings('solver', 'fmincon', 'verbose', 0);
% Solve the local optimization problems for each region
sol = optimize(Constraints, Objective, options);