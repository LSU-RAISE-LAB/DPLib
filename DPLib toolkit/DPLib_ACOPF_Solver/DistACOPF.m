function out = DistACOPF(partitionedDataFile, rho, residualThreshold, maxIterations, centralizedCost, holds)

% Distributed ACOPF Solver - Clean Version
clear global

% =========================================================================
% Title       : Distributed AC Optimal Power Flow (ACOPF) Solver using ADMM
% Author      : Milad Hasanzadeh
% Email       : e.mhasanzadeh1377@yahoo.com
% Affiliation : Department of Electrical and Computer Engineering,
%               Louisiana State University, Baton Rouge, LA, USA
% Date        : May 29, 2025
%
% Description :
% This MATLAB script implements a distributed AC Optimal Power Flow (ACOPF) 
% solver using the Alternating Direction Method of Multipliers (ADMM).
% ...
% (header truncated for brevity – unchanged)
% =========================================================================


%% ----------------------------- Load Case -------------------------------

% Append .mat if not already included
if ~endsWith(partitionedDataFile, '.mat')
    partitionedDataFile = [partitionedDataFile, '.mat'];
end

% Handle optional centralizedCost
if nargin < 5 || isempty(centralizedCost)
    centralizedCost = 0;
end

% Handle optional holds
if nargin < 6 || isempty(holds)
    holds = 1*10^(8);   % default value (same as your old internal value)
end


global loadedData
global regnum

% Load partitioned data
loadedData = load(partitionedDataFile);
numRegions = loadedData.num_regions;

% --- ADMM globals (now with separate rho for angle and magnitude) ---
global rho_theta rho_V
global PuScale
global landa      % duals for angle consensus
global lambdav    % duals for magnitude consensus


%% ------------- Augment Each Region with Interregional Tie-Lines -------------

for regionIdx = 1:numRegions
    regionID   = ['R' num2str(regionIdx)];
    regionName = ['mpc_region', regionID];

    tielinesName = ['interregional_tielines', regionID];
    baseBranch   = loadedData.(regionName);

    dataMatrix  = table2cell(loadedData.interregional_tielines_total);
    row_indices = strcmp(dataMatrix(:, 1), regionName) | strcmp(dataMatrix(:, 2), regionName);
    filteredData = dataMatrix(row_indices, :);

    nbus = size(loadedData.(regionName).bus,1);
    checkmattrix = zeros(size(filteredData,1),1);

    for i = 1:size(filteredData,1)
        region1_str = filteredData{i, 1};
        region2_str = filteredData{i, 2};

        region1_num = sscanf(region1_str, 'mpc_regionR%d');
        region2_num = sscanf(region2_str, 'mpc_regionR%d');

        if region1_num == regionIdx
            checkmattrix(i)=1;
        else
            checkmattrix(i)=-1; 
        end
    end

    output_rows = cell2mat(filteredData(:, 3:end));

    for i = 1:size(filteredData,1)
        if checkmattrix(i)==1
            output_rows(i,2)=nbus+i;
        else
            output_rows(i,1)=nbus+i;
        end
    end

    for i = 1:size(output_rows, 1)
        filteredData{i, 3} = output_rows(i, :);  % Replace row i in column 3 with 1×13 double row
    end
    loadedData.(tielinesName) = filteredData;

    % Construct struct dynamically
    mpc_struct = struct;
    mpc_struct.version = '2';

    % Base MVA
    if isfield(baseBranch, 'baseMVA')
        mpc_struct.baseMVA = baseBranch.baseMVA;
    else
        mpc_struct.baseMVA = 100;
    end

    % Bus data
    if isfield(baseBranch, 'bus')
        mpc_struct.bus = baseBranch.bus;
    end

    % Generator data
    if isfield(baseBranch, 'gen')
        mpc_struct.gen = baseBranch.gen;
    end

    % Branch data
    if isfield(baseBranch, 'branch')
        mpc_struct.branch = baseBranch.branch;
    end

    % Generator cost data
    if isfield(baseBranch, 'gencost')
        mpc_struct.gencost = baseBranch.gencost;
    end

    % Assign to workspace with dynamic name
    structName = ['mpc', regionID];
    assignin('base', structName, mpc_struct);

    % Dynamic variable names
    Gname      = ['G', regionID];
    Bname      = ['B', regionID];
    checkname  = ['check', regionID];
    tielinesName = ['interregional_tielines', regionID];

    % Load tie lines baseBranch
    tieLines = output_rows;
    LL       = tieLines;

    nbranches = size(LL, 1);
    G = zeros(nbranches, 4);
    B = zeros(nbranches, 4);

    for i = 1:nbranches
        r     = LL(i, 3);
        x     = LL(i, 4);
        bc    = LL(i, 5);
        ratio = LL(i, 9);
        if ratio == 0
            ratio = 1;
        end
        angle     = LL(i, 10);
        angle_rad = pi * angle / 180;

        invratio2 = 1 / ratio^2;
        multtf    = 1 / (ratio * exp(1j * angle_rad));
        multft    = 1 / (ratio * exp(-1j * angle_rad));
        z         = r + 1j * x;
        y         = 1 / z;

        Yff = (y + bc / 2 * 1j) * invratio2;
        Yft = -y * multft;
        Ytf = -y * multtf;
        Ytt = y + bc / 2 * 1j;

        G(i, 1) = real(Yff);  B(i, 1) = imag(Yff);
        G(i, 2) = real(Yft);  B(i, 2) = imag(Yft);
        G(i, 3) = real(Ytf);  B(i, 3) = imag(Ytf);
        G(i, 4) = real(Ytt);  B(i, 4) = imag(Ytt);
    end

    % Define the direction check: 1 if from bus > to bus, else 0
    check = double(LL(:, 1) > LL(:, 2));

    % Save G, B, and check matrices in loadedData
    loadedData.(Gname)     = G;
    loadedData.(Bname)     = B;
    loadedData.(checkname) = check;
end


%% ---------------------- Initialize Region Variables ----------------------
totalInterfaces = 0;
Flag_rho_theta=0;Flag_rho_v=0;
for regionIdx = 1:numRegions
    regionID   = ['R' num2str(regionIdx)];
    regionName = ['mpc_region', regionID];
    regionData = loadedData.(regionName);
    tielinesName = ['interregional_tielines', regionID];
    tielinesData = cell2table(loadedData.(tielinesName));
    tieLines     = tielinesData.Var3;

    % Basic attributes
    nb  = size(regionData.bus, 1);
    nit = size(tielinesData, 1);
    phase_angles = zeros(2 * nit, 1);
    voltage_mag  = zeros(2 * nit, 1);

    % Assign Yon and con as empty arrays initially
    Yon = [];
    con = [];

    % Get tie line bus IDs
    regionBusIDs = cell(size(tieLines, 1), 1);
    for j = 1:size(tieLines, 1)
        if tieLines(j, 1) <= nb
            regionBusIDs{j} = tieLines(j, 1);
        elseif tieLines(j, 2) <= nb
            regionBusIDs{j} = tieLines(j, 2);
        end
    end
    regionBusIDs = regionBusIDs(~cellfun('isempty', regionBusIDs));

    % Construct field names
    nbName       = ['nb', regionID];
    nitName      = ['nit', regionID];
    conName      = ['con', regionID];
    yonName      = ['Yon', regionID];
    phaseName    = ['phase_angles_', regionID];
    magnitudeName = ['voltage_mag_', regionID];
    regionVar    = ['region_', regionID];

    % Assign to loadedData structure
    loadedData.(nbName)       = nb;
    loadedData.(nitName)      = nit;
    loadedData.(conName)      = con;
    loadedData.(yonName)      = Yon;
    loadedData.(phaseName)    = phase_angles;
    loadedData.(magnitudeName)= voltage_mag;
    loadedData.(regionVar)    = regionBusIDs;

    % Update total interfaces
    totalInterfaces = totalInterfaces + nit;
end

% Total inter-regional interfaces
nit = totalInterfaces / 2;
loadedData.nit = nit;


%% ------------------------ Get ADMM Parameters -------------------------

baseMV = loadedData.mpc_regionR1.baseMVA;
mpc0   = loadcase(loadedData.filename);
nb_full = size(mpc0.bus, 1); %#ok<NASGU>

% Here 'rho' (function input) is the *base* penalty;
% we create separate scaled ones for angle and magnitude:
rho_theta = rho;%* baseMV* baseMV;%  * baseMV;
rho_V     = rho;%* baseMV* baseMV;% * baseMV * baseMV;

% ---- Angle & voltage scaling for normalized residuals ----
x     = mpc0.branch(:, 4);
rateA = mpc0.branch(:, 6);
Vmin  = min(mpc0.bus(:,13));
Vmax  = max(mpc0.bus(:,12));

global V_scale
global theta_scale
V_scale = (Vmax - Vmin)/2;   % usually ≈0.05–0.1 p.u.
iter_count_theta=0;iter_count_v=0;
theta_max_approx = abs(x) .* (rateA / baseMV);
valid = theta_max_approx > 0 & ~isnan(theta_max_approx) & ~isinf(theta_max_approx);
if any(valid)
    theta_scale = median(theta_max_approx(valid));   % e.g., ~0.1–0.3 rad typically
else
    theta_scale = 0.1;
end


PuScale = 1;

gen_coeff0      = mpc0.gencost(:, [5 6 7]);
allQuadraticZero = all(gen_coeff0(:,1) == 0); %#ok<NASGU>

% Initialization
% Initialization
landa      = zeros(2 * nit, 1);   % angle duals
lambdav    = zeros(2 * nit, 1);   % voltage duals
errorValues  = zeros(2 * nit, 1);
errorValuesv = zeros(2 * nit, 1);
errorLog     = zeros(maxIterations, 1);
optimalityGap = zeros(maxIterations, 1);
grad       = zeros(nit, 1);       %#ok<NASGU>
gradDual   = zeros(nit, 1);       %#ok<NASGU>
gradv      = zeros(nit, 1);       %#ok<NASGU>
gradvDual  = zeros(nit, 1);       %#ok<NASGU>
currentObjective = 0;

% --- Residual vectors for adaptive rho (normalized) ---
r_theta_vec      = zeros(2 * nit, 1);   % current angle residuals (normalized)
r_theta_vec_old  = zeros(2 * nit, 1);   % previous angle residuals
r_V_vec          = zeros(2 * nit, 1);   % current magnitude residuals (normalized)
r_V_vec_old      = zeros(2 * nit, 1);   % previous magnitude residuals

% For adaptive rho logging
rho_theta_log    = zeros(maxIterations,1);
rho_V_log        = zeros(maxIterations,1);


% Heuristic parameters
mu_theta  = max(1,15/((numRegions^2)*0.2));   % balance factor for angle residuals
mu_V      = max(1,10/((numRegions^2)*0.2));   % balance factor for voltage residuals
tau_incr  = 2;    % rho increase factor
tau_decr  = .5;     % rho decrease factor
eta = 1;max_count=20;
if nb_full>2000
max_count=max_count*2;
end
%% ---------------- Prepare Interregional Interface Buses ----------------

tieLineTable = loadedData.interregional_tielines_total;
tieLines     = table2cell(tieLineTable);

% Ensure region order is ascending in each tie-line
for i = 1:nit
    region1_str = tieLines{i, 1};
    region2_str = tieLines{i, 2};

    region1_num = sscanf(region1_str, 'mpc_regionR%d');
    region2_num = sscanf(region2_str, 'mpc_regionR%d');

    if region1_num > region2_num
        [tieLines{i, 1}, tieLines{i, 2}] = deal(region2_str, region1_str);
    end
end
tieLineTable = cell2table(tieLines);

% Initialize storage
interfaceCounter = zeros(numRegions, 1);
for regionIdx = 1:numRegions
    regionID = ['R', num2str(regionIdx)];
    loadedData.(['con', regionID]) = [];
    loadedData.(['Yon', regionID]) = [];
end

% Compute con and Yon for each region
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

    reg1_id = ['R', num2str(region1_num)];
    reg2_id = ['R', num2str(region2_num)];

    loadedData.(['con', reg1_id]) = [loadedData.(['con', reg1_id]); b];
    loadedData.(['con', reg2_id]) = [loadedData.(['con', reg2_id]); a];
    loadedData.(['Yon', reg1_id]) = [loadedData.(['Yon', reg1_id]); i];
    loadedData.(['Yon', reg2_id]) = [loadedData.(['Yon', reg2_id]); i];
end


%% -------------------------- Begin ADMM Loop ---------------------------

for iter = 1:maxIterations
    currentObjective = 0;
    Pg_total = 0;

    % ---------- Local ACOPF solve per region ----------
    for r = 1:numRegions
        reg    = ['R' num2str(r)];
        regnum = r;

        nb      = loadedData.(['nb' reg]);
        nitall  = loadedData.(['nit' reg]);
        regionName = loadedData.(['region_' reg]);

        % Run ACOPF for this region (uses globals rho_theta, rho_V, landa, lambdav)
        [Optimized_Solution, gen_coeff, flag, obj, stat] = acopf(reg); %#ok<ASGLU>

        nvars = numel(Optimized_Solution) - 2 * nitall;

        phase_angles = zeros(2 * nitall, 1);
        voltage_mag  = zeros(2 * nitall, 1);

        for i = 1:nitall
            phase_angles(i, 1)        = Optimized_Solution(regionName{i} + nb);
            phase_angles(i+nitall, 1) = Optimized_Solution((2 * i) + nvars);

            voltage_mag(i, 1)         = Optimized_Solution(regionName{i});
            voltage_mag(i+nitall, 1)  = Optimized_Solution(((2 * i) - 1) + nvars);
        end

        loadedData.(['phase_angles_', reg]) = phase_angles;
        loadedData.(['voltage_mag_',  reg]) = voltage_mag;

        if flag ~= 0
            Pg = zeros(size(gen_coeff, 1), 1);
            for i = 1:numel(Pg)
                Pg(i) = Optimized_Solution(2 * nb + i);
            end
        else
            Pg = zeros(size(gen_coeff,1),1);
        end

        currentObjective = currentObjective + value(obj);
        Pg_total         = Pg_total + sum(value(Pg));
    end

    % -------------------------- Dual Update Step ------------------------
    interfaceCounter = zeros(numRegions, 1);

    angleGrad     = zeros(nit,1);
    angleGradDual = zeros(nit,1);
    magGrad       = zeros(nit,1);
    magGradDual   = zeros(nit,1);

    for i = 1:nit
        regionName1 = char(tieLines{i, 1});  % e.g., 'mpc_regionR3'
        regionName2 = char(tieLines{i, 2});  % e.g., 'mpc_regionR5'

        regionIdx1 = sscanf(regionName1, 'mpc_regionR%d');
        regionIdx2 = sscanf(regionName2, 'mpc_regionR%d');

        for regionIdx = 1:numRegions
            if regionIdx == regionIdx1
                interfaceCounter(regionIdx) = interfaceCounter(regionIdx) + 1;
                interfacePos1 = interfaceCounter(regionIdx);
            elseif regionIdx == regionIdx2
                interfaceCounter(regionIdx) = interfaceCounter(regionIdx) + 1;
                interfacePos2 = interfaceCounter(regionIdx);
            end
        end

        reg1 = ['R' num2str(regionIdx1)];
        reg2 = ['R' num2str(regionIdx2)];

        anglesReg1 = loadedData.(['phase_angles_' reg1]);
        anglesReg2 = loadedData.(['phase_angles_' reg2]);
        magsReg1   = loadedData.(['voltage_mag_' reg1]);
        magsReg2   = loadedData.(['voltage_mag_' reg2]);
        nitReg1    = loadedData.(['nit' reg1]);
        nitReg2    = loadedData.(['nit' reg2]);

        angleGrad(i)     = anglesReg1(interfacePos1)       - anglesReg2(nitReg2 + interfacePos2);
        angleGradDual(i) = anglesReg1(nitReg1 + interfacePos1) - anglesReg2(interfacePos2);

        magGrad(i)       = magsReg1(interfacePos1)         - magsReg2(nitReg2 + interfacePos2);
        magGradDual(i)   = magsReg1(nitReg1 + interfacePos1) - magsReg2(interfacePos2);

        if any(isnan([angleGrad(i), angleGradDual(i), magGrad(i), magGradDual(i)]))
            fprintf(' NaN detected at interface %d\n', i);
        end


        % Normalized residuals
        r_theta      = angleGrad(i)     / theta_scale;
        r_theta_dual = angleGradDual(i) / theta_scale;
        r_V          = magGrad(i)       / V_scale;
        r_V_dual     = magGradDual(i)   / V_scale;

        % --- Store normalized residuals in full vectors (2*nit each) ---
        r_theta_vec(i)      = r_theta;
        r_theta_vec(i+nit)  = r_theta_dual;
        r_V_vec(i)          = r_V;
        r_V_vec(i+nit)      = r_V_dual;

        % Dual updates with separate rho
        landa(i)        = landa(i)        + eta *rho_theta * r_theta;
        landa(i + nit)  = landa(i + nit)  + eta *rho_theta * r_theta_dual;
        lambdav(i)      = lambdav(i)      + eta *rho_V     * r_V;
        lambdav(i + nit)= lambdav(i + nit)+ eta *rho_V     * r_V_dual;

        % Unnormalized residuals (for backward-compatible errorLog)
        errorValues(i)        = abs(r_theta);
        errorValues(i + nit)  = abs(r_theta_dual);
        errorValuesv(i)       = abs(r_V);
        errorValuesv(i + nit) = abs(r_V_dual);

    end

    % ----------------------- Log and Check Errors -----------------------

       % ----------------------- Log and Check Errors -----------------------

    % Unnormalized error for convergence test (backward compatible):
    maxError1 = max(errorValues);
    maxError2 = max(errorValuesv);
    maxError  = max(maxError1, maxError2);
    errorLog(iter) = maxError;

    if centralizedCost ~= 0
        optimalityGap(iter) = abs(((currentObjective) - centralizedCost) / centralizedCost) * 100;
    else
        optimalityGap(iter) = NaN;
    end

    % --- Normalized primal residuals (theta & V) ---
    primal_theta = max(abs(r_theta_vec));
    primal_V     = max(abs(r_V_vec));

    % --- Dual residuals = change in residuals between iterations ---
    dual_theta   = norm(r_theta_vec - r_theta_vec_old, inf);
    dual_V       = norm(r_V_vec     - r_V_vec_old,     inf);


    if (Flag_rho_theta==0 && rho_theta<holds)
        
    % --- Adaptive rho for angle ---
    if primal_theta > mu_theta * dual_theta
        rho_theta     = rho_theta * tau_incr;Flag_rho_theta=1;
    elseif dual_theta > mu_theta * primal_theta
        rho_theta     = rho_theta / tau_decr;Flag_rho_theta=1;
    end

    else
        iter_count_theta=iter_count_theta+1;
    end


    if iter_count_theta==max_count
        iter_count_theta=1;Flag_rho_theta=0;
    end


     if (Flag_rho_v==0 && rho_V<holds)
        

    if primal_V > mu_V * dual_V
        rho_V     = rho_V * tau_incr;Flag_rho_v=1;
    elseif dual_V > mu_V * primal_V
        rho_V     = rho_V / tau_decr;Flag_rho_v=1;
    end
    else
        iter_count_v=iter_count_v+1;
    end


    if iter_count_v==max_count
        iter_count_v=1;Flag_rho_v=0;
    end


    % Save residuals for next iteration + log rhos
    r_theta_vec_old   = r_theta_vec;
    r_V_vec_old       = r_V_vec;
    rho_theta_log(iter)  = rho_theta;
    rho_V_log(iter)      = rho_V;


    % Convergence check (still using unnormalized maxError for consistency)
    if maxError < residualThreshold
        fprintf('Convergence achieved at iteration %d\n', iter);
        if centralizedCost ~= 0
            fprintf('Iter %3d | Error: %.6e | Obj: %.4f | Gap: %.3f %% | rho_theta=%.3g | rho_V=%.3g\n', ...
                iter, maxError, currentObjective, optimalityGap(iter), rho_theta, rho_V);
        else
            fprintf('Iter %3d | Error: %.6e | Obj: %.4f | rho_theta=%.3g | rho_V=%.3g\n', ...
                iter, maxError, currentObjective, rho_theta, rho_V);
        end
        break;
    end

    % Print objective at this iteration
    if centralizedCost ~= 0
        fprintf('Iter %3d | Error: %.6e | Obj: %.4f | Gap: %.3f %% | rho_theta=%.3g | rho_V=%.3g\n', ...
            iter, maxError, currentObjective, optimalityGap(iter), rho_theta, rho_V);
    else
        fprintf('Iter %3d | Obj = %.6f | MaxError = %.6f | rho_theta=%.3g | rho_V=%.3g\n', ...
            iter, currentObjective, maxError, rho_theta, rho_V);
    end
end

% Check for non-convergence
if iter >= maxIterations
    disp('Maximum iterations reached without full convergence');
end
fprintf('Final Objective after %d iterations = %.6f\n', iter, currentObjective);

%% ---------------------- Plot Error and Gap Results ---------------------

figure;
semilogy(1:iter, errorLog(1:iter), 'b-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('Error'); title('Error Progression'); grid on; set(gca,'FontSize', 22);

if centralizedCost~=0
    figure;
    semilogy(1:iter, optimalityGap(1:iter), 'k-', 'LineWidth', 2);
    xlabel('Iteration'); ylabel('Optimality Gap'); title('Optimality Gap Progression'); grid on; set(gca,'FontSize', 22);
end

% Optional: plot rho_theta and rho_V evolution (comment out if not needed)
figure;
semilogy(1:iter, rho_theta_log(1:iter), 'LineWidth', 2); hold on;
semilogy(1:iter, rho_V_log(1:iter), '--', 'LineWidth', 2);
xlabel('Iteration', 'Interpreter','latex');
ylabel('$\rho$', 'Interpreter','latex');
legend({'$\rho_\theta$','$\rho_V$'}, 'Interpreter','latex','Location','best');
title('Penalty Parameter Evolution', 'Interpreter','latex');
grid on; set(gca,'FontSize',22);

%% ---------------------- Demand/Generation Summary ----------------------

mpc = loadcase(loadedData.filename);
Pd  = mpc.bus(:, 3);  % Real power demand (MW)
total_Pd = sum(Pd);
total_Pg = (mpc.baseMVA).*Pg_total;

fprintf('Total Real Power Demand (Pd): %.2f MW\n', total_Pd);
fprintf('Total Real Power Generation (Pg): %.2f MW\n', total_Pg);

%% ----------------------- Pack outputs in a struct ----------------------

out = struct();
out.partitionedDataFile = partitionedDataFile;
out.rho_input          = rho;        % the base input rho
out.rho_theta_final    = rho_theta;  % final angle rho
out.rho_V_final        = rho_V;      % final voltage rho
out.residualThreshold  = residualThreshold;
out.maxIterations      = maxIterations;
out.iterationsRun      = iter;
out.errorLog           = errorLog(1:iter);
out.optimalityGap      = optimalityGap(1:iter);
out.rho_theta_log      = rho_theta_log(1:iter);
out.rho_V_log          = rho_V_log(1:iter);
out.landa              = landa;
out.lambdav            = lambdav;
out.total_Pd           = total_Pd;
out.total_Pg           = total_Pg;
out.currentObjective   = currentObjective;
out.loadedData         = loadedData;

end  % end of DistACOPF function

%% Helper function (unchanged) ------------------------------------------------

function fprintf_matrix(fid, matrix)
    [rows, cols] = size(matrix);
    for i = 1:rows
        fprintf(fid, '    ');
        for j = 1:cols
            fprintf(fid, '%g\t', matrix(i,j));
        end
        fprintf(fid, '\n');
    end
end
