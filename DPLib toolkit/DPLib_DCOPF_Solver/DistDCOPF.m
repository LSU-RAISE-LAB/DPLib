% =========================================================================
% Title       : Distributed DC Optimal Power Flow (DCOPF) Solver using ADMM
% Author      : Milad Hasanzadeh
% Email       : e.mhasanzadeh1377@yahoo.com
% Affiliation : Department of Electrical and Computer Engineering,
%               Louisiana State University, Baton Rouge, LA, USA
% Date        : May 29, 2025
%
% Description :
% This MATLAB script implements a distributed DC Optimal Power Flow (DCOPF) 
% solver based on the Alternating Direction Method of Multipliers (ADMM). 
% The method supports large-scale and privacy-preserving OPF computation 
% by solving regional subproblems in parallel and coordinating through tie-line
% consistency enforcement.
%
% Functional Overview :
% - Loads partitioned MATPOWER-compatible case data.
% - Initializes regional variables and solves local DCOPF problems iteratively.
% - Performs dual variable (lambda) updates to maintain inter-region consistency.
% - Monitors convergence via primal residuals and an optional optimality gap.
% - Accepts user input for ADMM tuning parameters (penalty ρ, residual threshold, etc.).
% - Visualizes the convergence process using residual and optimality gap plots.
%
% Key Features :
% - Fully supports systems partitioned into multiple regions.
% - Seamlessly integrates with graph-based partitioning outputs.
% - Modular and script-based interface designed for academic research use.
% - Easily customizable stopping criteria and convergence logic.
%
% Requirements :
% - MATLAB R2020a or newer.
% - YALMIP toolbox installed and accessible in the MATLAB path (Download from: https://yalmip.github.io/download/).
% - MATPOWER toolbox properly configured and on the MATLAB path.
%
% Usage :
% - Run this script or call `mainDCOPF()` directly.
% - The user is prompted to input:
%     * The partitioned system data file
%     * The ADMM penalty parameter (ρ)
%     * The residual error threshold
%     * The maximum number of ADMM iterations
%     * Optionally, the centralized optimal cost (to compute the optimality gap)
%
% License :
% This script is intended strictly for academic and research purposes.
% Redistribution or commercial use without explicit permission is prohibited.
% If used in any publication, proper acknowledgment of the author is required.
%
% Copyright (c) 2025, Milad Hasanzadeh
% =========================================================================



%% Distributed DCOPF Solver - Clean Version

clear;close all;clc;


%% ----------------------------- Load Case -------------------------------
% Prompt user to input the partitioned .mat file
% Append .mat if not included, and verify file existence


partitionedDataFile = input('Enter partitioned case data file name: ', 's');

% Append .mat if not already included
if ~endsWith(partitionedDataFile, '.mat')
    partitionedDataFile = [partitionedDataFile, '.mat'];
end

% Check if the file exists
% if ~isfile(partitionedDataFile)
%     error('Required file "%s" not found. Please ensure it is in the current directory.', partitionedDataFile);
% end
% Load partitioned data
loadedData = load(partitionedDataFile);
numRegions = loadedData.num_regions;
totalInterfaces = 0;

%% ------------- Augment Each Region with Interregional Tie-Lines -------------
% Add tie-lines to each region’s branch matrix for distributed computation


    for regionIdx = 1:numRegions
        regionID =  ['R' num2str(regionIdx)];
        regionName = ['mpc_region', regionID];
        tielinesName = ['interregional_tielines', regionID];
        baseBranch = loadedData.(regionName).branch;

        dataMatrix = table2cell(loadedData.interregional_tielines_total);
        row_indices = strcmp(dataMatrix(:, 1), regionName) | strcmp(dataMatrix(:, 2), regionName);
        filteredData = dataMatrix(row_indices, :);
        nbus = size(loadedData.(regionName).bus,1);

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
checkmattrix=checkmattrix';
        output_rows = cell2mat(filteredData(:, 3:end));

        for i = 1:size(filteredData,1)

            if checkmattrix(i)==1
               output_rows(i,2)=nbus+i;
            else
               output_rows(i,1)=nbus+i;
            end

        end

        tieLines = output_rows;
        loadedData.(regionName).branch = [baseBranch; tieLines];
        for i = 1:size(output_rows, 1)
    filteredData{i, 3} = output_rows(i, :);  % Replace row i in column 3 with 1×13 double row
        end

loadedData.(tielinesName) =filteredData;
        
    end

%% ---------------------- Initialize Region Variables ----------------------
% Prepare regional tie-line identifiers and preallocate memory for ADMM

for regionIdx = 1:numRegions
    regionID = ['R' num2str(regionIdx)];
    regionName = ['mpc_region', regionID];
    tielinesName = ['interregional_tielines', regionID];
    regionData = loadedData.(regionName);
    tielinesData = loadedData.(tielinesName);
    % Create variable names
    nbName = ['nb', regionID];
    nitName = ['nit', regionID];
    conName = ['con', regionID];
    yonName = ['Yon', regionID];
    phaseName = ['phase_angles_', regionID];
    regionVar = ['region_', regionID];
    RealpowerflowName = ['Realpowerflow_', regionID];
    % Extract region data
    nb = size(regionData.bus, 1);
    nit = size(tielinesData, 1);
    branchMatrix = regionData.branch;

    % Initialize values
    assignin('base', nbName, nb);
    assignin('base', nitName, nit);
    assignin('base', phaseName, zeros(2 * nit, 1));assignin('base', RealpowerflowName, zeros(nit, 1));
    assignin('base', conName, []);
    assignin('base', yonName, []);

    % Extract tie-line rows
    outRegionRows = branchMatrix(:,1) > nb | branchMatrix(:,2) > nb;
    tieLineData = branchMatrix(outRegionRows, :);

    regionBusIDs = cell(size(tieLineData, 1), 1);
    for j = 1:size(tieLineData, 1)
        if tieLineData(j, 1) <= nb
            regionBusIDs{j} = tieLineData(j, 1);
        elseif tieLineData(j, 2) <= nb
            regionBusIDs{j} = tieLineData(j, 2);
        end
    end
    regionBusIDs = regionBusIDs(~cellfun('isempty', regionBusIDs));
    assignin('base', regionVar, regionBusIDs);

    totalInterfaces = totalInterfaces + nit;
end

% Total number of interfaces
nit = totalInterfaces / 2;
assignin('base', 'nit', nit);
clc;


%% ------------------------ Get ADMM Parameters -------------------------
% Ask user to enter ADMM penalty (rho), stopping threshold, max iterations, and central cost

rho = input('Enter penalty parameter (rho): ');PuScale = 10;
residualThreshold = input('Enter worst primal residual threshold to be used as stopping criterion: ');
maxIterations = input('Enter maximum iterations as stopping criterion:');
centralizedCost = input('Enter the centralized total cost in dollars to be used for optimality gap calculation, or enter 0 to skip the calculation:');


lambda = ones(2 * nit, 1);
errorValues = zeros(2 * nit, 1);
errorLog = zeros(maxIterations, 1);
optimalityGap = zeros(maxIterations, 1);
grad = zeros(nit, 1);
gradDual = zeros(nit, 1);
currentObjective = 0;

%% ---------------- Prepare Interregional Interface Buses ----------------
% Identify which buses participate in tie-lines for ADMM updates


% Process and reorder interregional tie-lines
tieLineTable = loadedData.interregional_tielines_total;
tieLines = table2cell(tieLineTable);

for i = 1:nit
    % Extract numeric region index from strings like 'mpc_regionR1'
    region1_str = tieLines{i, 1};
    region2_str = tieLines{i, 2};

    region1_num = sscanf(region1_str, 'mpc_regionR%d');
    region2_num = sscanf(region2_str, 'mpc_regionR%d');

    if region1_num > region2_num
        [tieLines{i, 1}, tieLines{i, 2}] = deal(region2_str, region1_str);
    end
end

tieLineTable = cell2table(tieLines);

% Precompute interconnection counters
interfaceCounter = zeros(numRegions, 1);
for i = 1:nit
    region1_str = tieLines{i, 1};  % e.g., 'mpc_regionR3'
    region2_str = tieLines{i, 2};  % e.g., 'mpc_regionR5'

    region1_num = sscanf(region1_str, 'mpc_regionR%d');  % returns 3
    region2_num = sscanf(region2_str, 'mpc_regionR%d');  % returns 5

    for j = 1:numRegions
        if j == region1_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            a = interfaceCounter(j);
        elseif j == region2_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            b = interfaceCounter(j);
        end
    end

    evalin('base', ['conR' num2str(region1_num) ' = [conR' num2str(region1_num) '; b];']);
    evalin('base', ['conR' num2str(region2_num) ' = [conR' num2str(region2_num) '; a];']);
    evalin('base', ['YonR' num2str(region1_num) ' = [YonR' num2str(region1_num) '; i];']);
    evalin('base', ['YonR' num2str(region2_num) ' = [YonR' num2str(region2_num) '; i];']);
end
clc;

%% -------------------------- Begin ADMM Loop ---------------------------
% Solve local DCOPF per region, update dual variables, and check for convergence

for iter = 1:maxIterations
    currentObjective = 0;Pg_total=0;
    for r = 1:numRegions
        reg = ['R' num2str(r)];
        run('DCOPF_SP');  % pass the region string


        nitVal = evalin('base', ['nit', reg]);
        nbVal = evalin('base', ['nb', reg]);
        regionBus = evalin('base', ['region_', reg]);
        phaseName = ['phase_angles_', reg];
        RealpowerflowName = ['Realpowerflow_', reg];
        % Extract updated angles
        for k = 1:nitVal
            angle(k, 1) = value(Delta(regionBus{k}));
            angle(k + nitVal, 1) = value(Delta(k + nbVal));
            Real_p(k,1) = value(Delta(regionBus{k})) - value(Delta(k + nbVal));
        end
        assignin('base', phaseName, angle); assignin('base', RealpowerflowName, Real_p);

        % Objective cost
        obj = 0;
        for g = 1:size(gen_coeff, 1)
            coeff = gen_coeff(g, :);
            obj = obj + coeff(1)*Pg(g)^2 + coeff(2)*Pg(g) + coeff(3);
        end
        currentObjective = currentObjective + value(obj);Pg_total = Pg_total + sum(value(Pg));
    end

% -------------------------- Dual Update Step ---------------------------
% Compute consistency gradients and perform Lagrangian updates (lambda)

interfaceCounter = zeros(numRegions, 1);

for i = 1:nit
    % Get full region names from tieLines
    str1 = tieLines{i, 1};  % e.g., 'mpc_regionR3'
    str2 = tieLines{i, 2};  % e.g., 'mpc_regionR7'

    % Extract region indices
    region1_num = sscanf(str1, 'mpc_regionR%d');
    region2_num = sscanf(str2, 'mpc_regionR%d');

    % Update interface counters for each region
    for j = 1:numRegions
        if j == region1_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            a = interfaceCounter(j);
        elseif j == region2_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            b = interfaceCounter(j);
        end
    end

    % Construct region variable names
    reg1 = ['R' num2str(region1_num)];
    reg2 = ['R' num2str(region2_num)];

    angle1 = evalin('base', ['phase_angles_' reg1]);
    angle2 = evalin('base', ['phase_angles_' reg2]);
    nit1 = evalin('base', ['nit' reg1]);
    nit2 = evalin('base', ['nit' reg2]);

    % Compute gradients
    grad(i) = angle1(a) - angle2(nit2 + b);
    gradDual(i) = angle1(nit1 + a) - angle2(b);

    % Dual updates
    lambda(i) = lambda(i) + rho * grad(i);
    lambda(i + nit) = lambda(i + nit) + rho * gradDual(i);

    % Track errors
    errorValues(i) = abs(grad(i));
    errorValues(i + nit) = abs(gradDual(i));
end


% ----------------------- Log and Check Errors -------------------------
% Track primal residuals and compute optimality gap if centralized cost is known

    maxError = max(errorValues);
    errorLog(iter) = maxError;
    optimalityGap(iter) = abs((currentObjective - centralizedCost) / centralizedCost) * 100;

    % Convergence check
    if maxError < residualThreshold
        fprintf('Convergence achieved at iteration %d\n', iter);
        break;
    end
    fprintf('Iteration %d | Error: %.6f \n', iter, maxError);
end

% Check for non-convergence
if iter >= maxIterations
    disp('Maximum iterations reached without full convergence');
end

%% ---------------------- Plot Error and Gap Results ---------------------
% Display convergence performance using error and optimality plots

figure;
semilogy(1:iter, errorLog(1:iter), 'b-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('Error'); title('Error Progression'); grid on; set(gca,'FontSize', 22);
if centralizedCost~=0
% Plot Optimality Gap
figure;
semilogy(1:iter, optimalityGap(1:iter), 'k-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('Optimality Gap'); title('Optimality Gap Progression'); grid on; set(gca,'FontSize', 22);
end


% Load MATPOWER case
mpc = loadcase(loadedData.filename);

% Extract the demand data
Pd = mpc.bus(:, 3);  % Real power demand (MW)
% Sum the total power demands
total_Pd = sum(Pd);  % Total real power demand
total_Pg = (mpc.baseMVA).*Pg_total;
% Display results
fprintf('Total Real Power Demand (Pd): %.2f MW\n', total_Pd);
fprintf('Total Real Power Generation (Pg): %.2f MW\n', total_Pg);

