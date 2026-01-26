function out = DistDCOPF(partitionedDataFile, rho, residualThreshold, maxIterations, centralizedCost, holds)


% =========================================================================
% Title       : Distributed DC Optimal Power Flow (DCOPF) Solver using ADMM
% Author      : Milad Hasanzadeh
% Email       : e.mhasanzadeh1377@yahoo.com
% Affiliation : Department of Electrical and Computer Engineering,
%               Louisiana State University, Baton Rouge, LA, USA
% Date        : May 29, 2025
%
% Description :
% This function implements a distributed DC Optimal Power Flow (DCOPF) 
% solver based on the Alternating Direction Method of Multipliers (ADMM). 
% The method supports large-scale and privacy-preserving OPF computation 
% by solving regional subproblems in parallel and coordinating through 
% tie-line consistency enforcement.
%
% FUNCTION SIGNATURE:
%   out = DistDCOPF(partitionedDataFile, rho, residualThreshold, maxIterations, centralizedCost)
%
% INPUTS:
%   partitionedDataFile : string, name of the partitioned .mat file (with or without .mat)
%   rho                 : ADMM penalty parameter
%   residualThreshold   : stopping threshold on worst primal residual
%   maxIterations       : maximum number of ADMM iterations
%   centralizedCost     : centralized optimal total cost (dollars); 
%                         set to 0 to skip optimality gap computation
%
% OUTPUT:
%   out : struct containing key results:
%       .converged          - logical, true if converged before maxIterations
%       .iterations         - number of iterations actually run
%       .errorLog           - vector of max error per iteration
%       .optimalityGap      - vector of optimality gap per iteration (NaN if centralizedCost==0)
%       .lambda             - final dual variables (2*nit x 1)
%       .currentObjective   - objective value at final iteration
%       .tieLines           - tie-line region pairing cell array
%       .total_Pd           - total demand (MW)
%       .total_Pg           - total generation (MW, as computed in script)
%       .nit                - total number of interfaces
%       .numRegions         - number of regions
%
% NOTE:
%   This function still uses evalin/assignin to interact with DCOPF_SP,
%   as in the original script.
% =========================================================================


%% ----------------------------- Load Case -------------------------------
% Append .mat if not already included
if ~endsWith(partitionedDataFile, '.mat')
    partitionedDataFile = [partitionedDataFile, '.mat'];
end

% Check if the file exists
if ~isfile(partitionedDataFile)
    error('Required file "%s" not found. Please ensure it is in the current directory.', ...
        partitionedDataFile);
end

% Load partitioned data
loadedData   = load(partitionedDataFile);
numRegions   = loadedData.num_regions;
totalInterfaces = 0;

% Handle optional centralizedCost (if someone calls with 4 inputs)
if nargin < 5 || isempty(centralizedCost)
    centralizedCost = 0;
end

% Handle optional holds (if someone calls without it)
if nargin < 6 || isempty(holds)
    holds = 1*10^(8);   % default value (you can pick whatever you like)
end




%% ------------- Augment Each Region with Interregional Tie-Lines --------
% Add tie-lines to each region’s branch matrix for distributed computation

for regionIdx = 1:numRegions
    regionID      = ['R' num2str(regionIdx)];
    regionName    = ['mpc_region', regionID];
    tielinesName  = ['interregional_tielines', regionID];
    baseBranch    = loadedData.(regionName).branch;

    % Convert global interregional table to cell
    dataMatrix = table2cell(loadedData.interregional_tielines_total);

    % Rows where this region participates in the tie-lines
    row_indices = strcmp(dataMatrix(:, 1), regionName) | ...
                  strcmp(dataMatrix(:, 2), regionName);

    filteredData = dataMatrix(row_indices, :);
    nbus        = size(loadedData.(regionName).bus, 1);
    nTieLocal   = size(filteredData, 1);

    % Safety: handle regions with no tie-lines
    if nTieLocal == 0
        loadedData.(tielinesName) = cell(0, size(filteredData, 2));
        continue;
    end

    % Initialize check matrix per region (direction indicator)
    checkmattrix = zeros(nTieLocal, 1);

    for i = 1:nTieLocal
        region1_str = filteredData{i, 1};
        region2_str = filteredData{i, 2};

        region1_num = sscanf(region1_str, 'mpc_regionR%d');
        region2_num = sscanf(region2_str, 'mpc_regionR%d');

        if region1_num == regionIdx
            checkmattrix(i) = 1;
        else
            checkmattrix(i) = -1;
        end
    end

    % Numeric tie-line rows (branch columns)
    output_rows = cell2mat(filteredData(:, 3:end));

    % Adjust from/to bus indices for out-of-region buses
    for i = 1:nTieLocal
        if checkmattrix(i) == 1
            % This region is "region1": to-bus is external
            output_rows(i, 2) = nbus + i;
        else
            % This region is "region2": from-bus is external
            output_rows(i, 1) = nbus + i;
        end
    end

    % Append tie-lines to region branch matrix
    tieLines = output_rows;
    loadedData.(regionName).branch = [baseBranch; tieLines];

    % Store the numeric tie-line row in column 3 for this region’s tie-line cell
    for i = 1:nTieLocal
        filteredData{i, 3} = output_rows(i, :);  % 1×13 double
    end

    % Save region-specific tie-line data
    loadedData.(tielinesName) = filteredData;
end

%% ---------------------- Initialize Region Variables --------------------
% Prepare regional tie-line identifiers and preallocate memory for ADMM

for regionIdx = 1:numRegions
    regionID      = ['R' num2str(regionIdx)];
    regionName    = ['mpc_region', regionID];
    tielinesName  = ['interregional_tielines', regionID];
    regionData    = loadedData.(regionName);
    tielinesData  = loadedData.(tielinesName);

    % Create variable names
    nbName        = ['nb', regionID];
    nitName       = ['nit', regionID];
    conName       = ['con', regionID];
    yonName       = ['Yon', regionID];
    phaseName     = ['phase_angles_', regionID];
    regionVar     = ['region_', regionID];
    RealpowerName = ['Realpowerflow_', regionID];

    % Extract region data
    nb  = size(regionData.bus, 1);
    nitLocal = size(tielinesData, 1);
    branchMatrix = regionData.branch;

    % Initialize values in base workspace (for DCOPF_SP script)
    assignin('base', nbName, nb);
    assignin('base', nitName, nitLocal);
    assignin('base', phaseName, zeros(2 * nitLocal, 1));
    assignin('base', RealpowerName, zeros(nitLocal, 1));
    assignin('base', conName, []);
    assignin('base', yonName, []);

    % Extract tie-line rows (buses connected to this region)
    outRegionRows = branchMatrix(:, 1) > nb | branchMatrix(:, 2) > nb;
    tieLineData   = branchMatrix(outRegionRows, :);

    regionBusIDs  = cell(size(tieLineData, 1), 1);
    for j = 1:size(tieLineData, 1)
        if tieLineData(j, 1) <= nb
            regionBusIDs{j} = tieLineData(j, 1);
        elseif tieLineData(j, 2) <= nb
            regionBusIDs{j} = tieLineData(j, 2);
        end
    end
    regionBusIDs = regionBusIDs(~cellfun('isempty', regionBusIDs));
    assignin('base', regionVar, regionBusIDs);

    totalInterfaces = totalInterfaces + nitLocal;
end

% Total number of interfaces (each interface appears twice across regions)
nit = totalInterfaces / 2;
assignin('base', 'nit', nit);
% Load the original MATPOWER case only to get baseMVA
% Load original MATPOWER case (for scaling and reporting)
mpc0   = loadcase(loadedData.filename);
nb_full = size(mpc0.bus, 1);
baseMW  = mpc0.baseMVA;
%rho = rho *baseMW*baseMW; 
% ---- Angle scaling for normalized residuals ----
% Approximate max admissible angle difference per line:
% |theta_from - theta_to| <= |x| * rateA / baseMVA
x      = mpc0.branch(:, 4);
rateA  = mpc0.branch(:, 6);

theta_max_approx = abs(x) .* (rateA / baseMW);

valid = theta_max_approx > 0 & ~isnan(theta_max_approx) & ~isinf(theta_max_approx);
if any(valid)
    theta_scale = median(theta_max_approx(valid));   % e.g., ~0.1–0.3 rad typically
else
    theta_scale = 0.1;   % safe fallback if no valid limits
end


%% ------------------------ Get ADMM Parameters -------------------------
% (Inputs are passed as arguments now, no interactive input.)
lambda        = zeros(2 * nit, 1);
errorValues   = zeros(2 * nit, 1);
errorLog      = zeros(maxIterations, 1);
optimalityGap = NaN(maxIterations, 1);   % NaN if centralizedCost == 0
grad          = zeros(nit, 1);
gradDual      = zeros(nit, 1);
currentObjective = 0;

% Residual vectors for adaptive rho
r_theta      = zeros(2 * nit, 1);   % current normalized residuals
r_theta_old  = zeros(2 * nit, 1);   % previous iteration residuals

% Heuristic parameters for adaptive rho
mu        = max(1,15/((numRegions^2)*0.2));   % residual balance factor
tau_incr  = 2;    % how much to increase rho
tau_decr  = .5;    % how much to decrease rho
eta = 1;
max_count=20;
if nb_full>2000
max_count=max_count*2;
if nb_full<3000
residualThreshold=residualThreshold*100;
end
end


rhoLog  = zeros(maxIterations, 1);   % to store rho per iteration
dualLog = zeros(maxIterations, 1);   % optional: log dual residual too
iter_count=0;
Flag_rho=0;
%% ---------------- Prepare Interregional Interface Buses ----------------
% Identify which buses participate in tie-lines for ADMM updates

tieLineTable = loadedData.interregional_tielines_total;
tieLines     = table2cell(tieLineTable);

% Reorder region strings so that region1_num < region2_num
for i = 1:nit
    region1_str = tieLines{i, 1};
    region2_str = tieLines{i, 2};

    region1_num = sscanf(region1_str, 'mpc_regionR%d');
    region2_num = sscanf(region2_str, 'mpc_regionR%d');

    if region1_num > region2_num
        [tieLines{i, 1}, tieLines{i, 2}] = deal(region2_str, region1_str);
    end
end

% Replace with sorted ordering
tieLineTable = cell2table(tieLines);

% Precompute interconnection counters (local interface index per region)
interfaceCounter = zeros(numRegions, 1);

for i = 1:nit
    region1_str = tieLines{i, 1};  % e.g., 'mpc_regionR3'
    region2_str = tieLines{i, 2};  % e.g., 'mpc_regionR5'

    region1_num = sscanf(region1_str, 'mpc_regionR%d');
    region2_num = sscanf(region2_str, 'mpc_regionR%d');

    for j = 1:numRegions
        if j == region1_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            a = interfaceCounter(j);
        elseif j == region2_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            b = interfaceCounter(j);
        end
    end

    evalin('base', sprintf('conR%d = [conR%d; %d];', region1_num, region1_num, b));
    evalin('base', sprintf('conR%d = [conR%d; %d];', region2_num, region2_num, a));
    evalin('base', sprintf('YonR%d = [YonR%d; %d];', region1_num, region1_num, i));
    evalin('base', sprintf('YonR%d = [YonR%d; %d];', region2_num, region2_num, i));
end

%% -------------------------- Begin ADMM Loop ----------------------------
fprintf('Starting ADMM iterations...\n');

for iter = 1:maxIterations
    currentObjective = 0;
    Pg_total = 0;

    % ----------- Local DCOPF Solve per Region ---------------------
    for r = 1:numRegions
        reg = ['R' num2str(r)];

        % DCOPF_SP is assumed to:
        %  - read reg, lambda, rho, conR*, YonR*, etc. from workspace
        %  - define Delta, Pg, gen_coeff, etc.

         
        run('DCOPF_SP.m');


        nitVal = evalin('base', ['nit', reg]);
        nbVal  = evalin('base', ['nb', reg]);
        regionBus = evalin('base', ['region_', reg]);
        phaseName = ['phase_angles_', reg];
        RealpowerName = ['Realpowerflow_', reg];

        % Reinitialize per region to avoid leftover values
        angle  = zeros(2 * nitVal, 1);
        Real_p = zeros(nitVal, 1);

        % Extract updated angles and real power "flow" proxy
        for k = 1:nitVal
            angle(k, 1)           = value(Delta(regionBus{k}));
            angle(k + nitVal, 1)  = value(Delta(k + nbVal));
            Real_p(k, 1)          = value(Delta(regionBus{k})) - ...
                                    value(Delta(k + nbVal));
        end
        assignin('base', phaseName, angle);
        assignin('base', RealpowerName, Real_p);

        % % Local objective cost
        % obj = 0;
        % for g = 1:size(gen_coeff, 1)
        %     coeff = gen_coeff(g, :);PgMW = Pg(g)*baseMW ;
        % 
        %         if coeff(1) == 3
        %             obj = obj + coeff(2)*PgMW^2 + coeff(3)*PgMW + coeff(4);
        %         elseif coeff(1) == 2
        %             obj = obj + coeff(2)*PgMW + coeff(3);
        %         end
        % 
        % end
        % currentObjective = currentObjective + value(obj);

            % --- Local objective cost using MATPOWER's totcost ---
    % Pg is in p.u., so convert to MW first
    Pg_MW = value(Pg)*baseMW;      % column vector [ngen_r x 1]

  
        % Compute regional generation cost (same definition as MATPOWER)
    cost_vec = totcost(mpc.gencost, Pg_MW);

    % Ensure scalar (totcost may return a vector)
    cost_reg = sum(cost_vec(:));

    % Accumulate global objective and total generation (still in p.u. here)
    currentObjective = currentObjective + cost_reg;
    Pg_total         = Pg_total + sum(value(Pg));

    end
       
    % --------------------- Dual Update Step ------------------------
    interfaceCounter = zeros(numRegions, 1);  % reset per iteration

    for i = 1:nit
        str1 = tieLines{i, 1};  % 'mpc_regionR#'
        str2 = tieLines{i, 2};

        region1_num = sscanf(str1, 'mpc_regionR%d');
        region2_num = sscanf(str2, 'mpc_regionR%d');

        for j = 1:numRegions
            if j == region1_num
                interfaceCounter(j) = interfaceCounter(j) + 1;
                a = interfaceCounter(j);
            elseif j == region2_num
                interfaceCounter(j) = interfaceCounter(j) + 1;
                b = interfaceCounter(j);
            end
        end

        reg1 = ['R' num2str(region1_num)];
        reg2 = ['R' num2str(region2_num)];

        angle1 = evalin('base', ['phase_angles_' reg1]);
        angle2 = evalin('base', ['phase_angles_' reg2]);
        nit1   = evalin('base', ['nit' reg1]);
        nit2   = evalin('base', ['nit' reg2]);

% Gradients (raw angle residuals in radians)
grad(i)     = angle1(a)        - angle2(nit2 + b);
gradDual(i) = angle1(nit1 + a) - angle2(b);

% ---- Normalize residuals by theta_scale ----
r_norm      = grad(i)     / theta_scale;
rdual_norm  = gradDual(i) / theta_scale;
% Store normalized residuals in r_theta (length 2*nit)
r_theta(i)      = r_norm;
r_theta(i+nit)  = rdual_norm;
% Dual updates with normalized residuals
lambda(i)      = lambda(i)      + eta *rho * r_norm;
lambda(i+nit)  = lambda(i+nit)  + eta *rho * rdual_norm;

% Track absolute *normalized* errors
errorValues(i)      = abs(r_norm);
errorValues(i+nit)  = abs(rdual_norm);

    end

    % ----------------------- Log and Check Errors -----------------
% ----------------------- Log and Check Errors -----------------
% Primal residual = max normalized consensus mismatch
primalRes       = max(abs(r_theta));
errorLog(iter)  = primalRes;

% Dual residual = change in residuals between iterations
dualRes         = norm(r_theta - r_theta_old, inf);

dualLog(iter)   = dualRes;
rhoLog(iter)    = rho;
    if (Flag_rho==0 && rho<holds)
        
% --- Adaptive rho update (Boyd-style heuristic on residuals) ---
if primalRes > mu * dualRes
    rho     = rho * tau_incr;Flag_rho=1;
elseif dualRes > mu * primalRes
    rho     = rho / tau_decr;Flag_rho=1;
end
    else
        iter_count=iter_count+1;
    end

    if iter_count==max_count
        iter_count=1;Flag_rho=0;
    end

% Store residual for next iteration
r_theta_old = r_theta;

if centralizedCost ~= 0
    optimalityGap(iter) = abs(((currentObjective) - centralizedCost) ...
                               / centralizedCost) * 100;
    fprintf('Iter %3d | Error: %.6e | Obj: %.4f | Gap: %.3f %% | rho: %.4f |\n', ...
        iter, primalRes, currentObjective, optimalityGap(iter), rho);
else
    fprintf('Iter %3d | Error: %.6e | Obj: %.4f | rho: %.4f |\n', ...
        iter, primalRes, currentObjective, rho);
end


    % Convergence check
    if primalRes < residualThreshold
        fprintf('Convergence achieved at iteration %d\n', iter);
        break;
    end
end

converged = (iter < maxIterations);

if ~converged
    disp('Maximum iterations reached without full convergence.');
end

%% ---------------------- Plot Error and Gap Results --------------------
figure;
semilogy(1:iter, errorLog(1:iter), 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Primal residual (max error)');
title('ADMM Error Progression');
grid on;
set(gca, 'FontSize', 22);

if centralizedCost ~= 0
    figure;
    semilogy(1:iter, optimalityGap(1:iter), 'k-', 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Optimality Gap [%]');
    title('Optimality Gap Progression');
    grid on;
    set(gca, 'FontSize', 22);
end
figure;
semilogy(1:iter, rhoLog(1:iter), 'LineWidth', 2);
xlabel('Iteration', 'Interpreter', 'latex');
ylabel('$\rho$', 'Interpreter', 'latex');   % <- LaTeX rho
title('ADMM Penalty Parameter $\rho$', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 22);


%% ---------------------- Post-Processing: Demand/Generation ------------
% Load MATPOWER case
mpc = loadcase(loadedData.filename);

% Extract the demand data
Pd = mpc.bus(:, 3);   % Real power demand (MW)
total_Pd = sum(Pd);   % Total real power demand (MW)

% If Pg_total is in MW already, you might NOT need to multiply by baseMVA.
% Kept as in your original logic.
total_Pg = (mpc.baseMVA) .* Pg_total;

fprintf('Total Real Power Demand (Pd): %.2f MW\n', total_Pd);
fprintf('Total Real Power Generation (Pg): %.2f MW\n', total_Pg);

%% ---------------------------- Build Output Struct ----------------------
out = struct();
out.converged        = converged;
out.iterations       = iter;
out.errorLog         = errorLog(1:iter);
out.optimalityGap    = optimalityGap(1:iter);
out.lambda           = lambda;
out.currentObjective = currentObjective;
out.tieLines         = tieLines;
out.total_Pd         = total_Pd;
out.total_Pg         = total_Pg;
out.nit              = nit;
out.numRegions       = numRegions;

end


