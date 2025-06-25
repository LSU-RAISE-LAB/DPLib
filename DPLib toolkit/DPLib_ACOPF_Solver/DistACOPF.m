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
% The method enables scalable and privacy-preserving power flow computation 
% by decomposing a large system into multiple electrically coherent regions 
% and solving subproblems locally, with coordination across tie-lines.
%
% Functional Overview :
% - Loads MATPOWER-compatible partitioned case data.
% - Initializes each region and solves ACOPF subproblems using IPOPT.
% - Updates dual variables (Lagrange multipliers) to ensure tie-line consistency.
% - Tracks convergence based on primal residual and optimality gap metrics.
% - Supports user-defined ADMM parameters such as penalty value (ρ), 
%   residual threshold, and maximum iteration count.
% - Visualizes convergence trends through error and optimality gap plots.
%
% Key Features :
% - Full support for systems partitioned into multiple regions.
% - Compatible with graph-based system partitioning frameworks.
% - Modular and script-driven structure, well-suited for research use.
% - Flexible convergence criteria and diagnostic outputs.
%
% Requirements :
% - MATLAB R2020a or newer.
% - MATPOWER toolbox (for data compatibility).
% - IPOPT optimizer installed and configured with MATLAB interface (Download IPOPT: https://coin-or.github.io/Ipopt/).
%  
%
% Usage :
% - Run the script or invoke `mainACOPF()` in MATLAB.
% - The user is prompted for:
%     * Penalty parameter (ρ)
%     * Convergence threshold (primal residual)
%     * Maximum number of iterations
%     * Centralized ACOPF benchmark cost (optional, for computing optimality gap)
%
% License :
% This script is provided for academic and research purposes only.
% Redistribution or commercial use without prior written permission is prohibited.
% Proper acknowledgment must be given if used in any published work.
%
% Copyright (c) 2025, Milad Hasanzadeh
% =========================================================================



%% Distributed ACOPF Solver - Clean Version

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


global loadedData
global regnum
% Load partitioned data
loadedData = load(partitionedDataFile);
numRegions = loadedData.num_regions;
global rho
global PuScale
global landa
global lambdav

%% ------------- Augment Each Region with Interregional Tie-Lines -------------

    for regionIdx = 1:numRegions
        regionID = ['R' num2str(regionIdx)];
        regionName = ['mpc_region', regionID];

        tielinesName = ['interregional_tielines', regionID];
        baseBranch = loadedData.(regionName);

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
                     for i = 1:size(output_rows, 1)
    filteredData{i, 3} = output_rows(i, :);  % Replace row i in column 3 with 1×13 double row
        end
loadedData.(tielinesName) =filteredData;
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
        Gname = ['G', regionID];
        Bname = ['B', regionID];
        checkname = ['check', regionID];
        tielinesName = ['interregional_tielines', regionID];

        % Load tie lines baseBranch
        tieLines = output_rows;
        LL = tieLines;

        nbranches = size(LL, 1);
        G = zeros(nbranches, 4);
        B = zeros(nbranches, 4);

        for i = 1:nbranches
            r = LL(i, 3);
            x = LL(i, 4);
            bc = LL(i, 5);
            ratio = LL(i, 9);
            if ratio == 0
                ratio = 1;
            end
            angle = LL(i, 10);
            angle_rad = pi * angle / 180;

            invratio2 = 1 / ratio^2;
            multtf = 1 / (ratio * exp(1j * angle_rad));
            multft = 1 / (ratio * exp(-1j * angle_rad));
            z = r + 1j * x;
            y = 1 / z;

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
        loadedData.(Gname) = G;
        loadedData.(Bname) = B;
        loadedData.(checkname) = check;
    end



%% ---------------------- Initialize Region Variables ----------------------
% Prepare regional tie-line identifiers and preallocate memory for ADMM
totalInterfaces = 0;

for regionIdx = 1:numRegions
    regionID = ['R' num2str(regionIdx)];
    regionName = ['mpc_region', regionID];
    regionData = loadedData.(regionName);
    tielinesName = ['interregional_tielines', regionID];
    tielinesData = cell2table(loadedData.(tielinesName));
    tieLines = tielinesData.Var3;

    % Basic attributes
    nb  = size(regionData.bus, 1);
    nit = size(tielinesData, 1);
    phase_angles = zeros(2 * nit, 1);
    voltage_mag = zeros(2 * nit, 1);

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
    nbName    = ['nb', regionID];
    nitName   = ['nit', regionID];
    conName   = ['con', regionID];
    yonName   = ['Yon', regionID];
    phaseName = ['phase_angles_', regionID];
    magnitudeName = ['voltage_mag_', regionID];
    regionVar = ['region_', regionID];

    % Assign to loadedData structure
    loadedData.(nbName) = nb;
    loadedData.(nitName) = nit;
    loadedData.(conName) = con;
    loadedData.(yonName) = Yon;
    loadedData.(phaseName) = phase_angles;
    loadedData.(magnitudeName) = voltage_mag;
    loadedData.(regionVar) = regionBusIDs;

    % Update total interfaces
    totalInterfaces = totalInterfaces + nit;
end

% Total inter-regional interfaces
nit = totalInterfaces / 2;
loadedData.nit = nit;
clc;


%% ------------------------ Get ADMM Parameters -------------------------
% Ask user to enter ADMM penalty (rho), stopping threshold, max iterations, and central cost


rho = input('Enter penalty parameter (rho): ');PuScale = 10*(loadedData.mpc_regionR1.baseMVA);
residualThreshold = input('Enter worst primal residual threshold to be used as stopping criterion: ');
maxIterations = input('Enter maximum iterations as stopping criterion:');
centralizedCost = input('Enter the centralized total cost in dollars to be used for optimality gap calculation, or enter 0 to skip the calculation:');

clc;


% Initialization
landa = ones(2 * nit, 1);lambdav = ones(2 * nit, 1);
errorValues = zeros(2 * nit, 1);errorValuesv = zeros(2 * nit, 1);
errorLog = zeros(maxIterations, 1);
optimalityGap = zeros(maxIterations, 1);
grad = zeros(nit, 1);
gradDual = zeros(nit, 1);
gradv = zeros(nit, 1);
gradvDual = zeros(nit, 1);
currentObjective = 0;


%% ---------------- Prepare Interregional Interface Buses ----------------
% Identify which buses participate in tie-lines for ADMM updates

tieLineTable = loadedData.interregional_tielines_total;
tieLines = table2cell(tieLineTable);

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

    % Append to loadedData
    reg1_id = ['R', num2str(region1_num)];
    reg2_id = ['R', num2str(region2_num)];

    loadedData.(['con', reg1_id]) = [loadedData.(['con', reg1_id]); b];
    loadedData.(['con', reg2_id]) = [loadedData.(['con', reg2_id]); a];
    loadedData.(['Yon', reg1_id]) = [loadedData.(['Yon', reg1_id]); i];
    loadedData.(['Yon', reg2_id]) = [loadedData.(['Yon', reg2_id]); i];
end


%% -------------------------- Begin ADMM Loop ---------------------------
% Solve local DCOPF per region, update dual variables, and check for convergence


for iter = 1:maxIterations
    currentObjective = 0;Pg_total=0;
    for r = 1:numRegions
        reg = ['R' num2str(r)];regnum = r;

        nb = loadedData.(['nb' reg]);
        nitall = loadedData.(['nit' reg]);
        regionName = loadedData.(['region_' reg]);

            % --- Run ACOPF for this region ---
    [Optimized_Solution, gen_coeff,flag] = acopf(reg);
    nvars = numel(Optimized_Solution) - 2 * nitall;
    
    phase_angles = zeros(2 * nitall, 1);
    voltage_mag = zeros(2 * nitall, 1);
    % Extract optimized values
    for i = 1:nitall
        phase_angles(i, 1)        = Optimized_Solution(regionName{i} + nb);
        phase_angles(i+nitall, 1)    = Optimized_Solution((2 * i) + nvars);
        
        voltage_mag(i, 1)         = Optimized_Solution(regionName{i});
        voltage_mag(i+nitall, 1)     = Optimized_Solution(((2 * i) - 1) + nvars);
    end
    % Store in loadedData instead of base workspace
    loadedData.(['phase_angles_', reg]) = phase_angles;
    loadedData.(['voltage_mag_',  reg]) = voltage_mag;
if flag~=0
    Pg = zeros(size(gen_coeff, 1), 1);
    for i = 1:numel(Pg)
        Pg(i) = Optimized_Solution(2 * nb + i);
    end
end
    % Compute subproblem objective
    SP_objective = 0;
    if flag~=0
    for i = 1:size(gen_coeff, 1)
        coeffs = gen_coeff(i, :);
        SP_objective = SP_objective + coeffs(1) * Pg(i)^2 + coeffs(2) * Pg(i) + coeffs(3);
    end
    end
    currentObjective = currentObjective + value(SP_objective);Pg_total = Pg_total + sum(value(Pg));


    end
% -------------------------- Dual Update Step ---------------------------
% Compute consistency gradients and perform Lagrangian updates (lambda)


    interfaceCounter = zeros(numRegions, 1);
    for i = 1:nit
    regionName1 = char(tieLines{i, 1});  % e.g., 'mpc_regionR3'
    regionName2 = char(tieLines{i, 2});  % e.g., 'mpc_regionR5'

    % Extract numeric region IDs
    regionIdx1 = sscanf(regionName1, 'mpc_regionR%d');
    regionIdx2 = sscanf(regionName2, 'mpc_regionR%d');

    % Update interface counters and determine interface positions
    for regionIdx = 1:numRegions
        if regionIdx == regionIdx1
            interfaceCounter(regionIdx) = interfaceCounter(regionIdx) + 1;
            interfacePos1 = interfaceCounter(regionIdx);
        elseif regionIdx == regionIdx2
            interfaceCounter(regionIdx) = interfaceCounter(regionIdx) + 1;
            interfacePos2 = interfaceCounter(regionIdx);
        end
    end

    % Build region key strings
    reg1 = ['R' num2str(regionIdx1)];
    reg2 = ['R' num2str(regionIdx2)];
    % Load region-specific data from loadedData
    anglesReg1 = loadedData.(['phase_angles_' reg1]);
    anglesReg2 = loadedData.(['phase_angles_' reg2]);
    magsReg1   = loadedData.(['voltage_mag_' reg1]);
    magsReg2   = loadedData.(['voltage_mag_' reg2]);
    nitReg1    = loadedData.(['nit' reg1]);
    nitReg2    = loadedData.(['nit' reg2]);

    % Compute interface gradient terms
    angleGrad(i)     = anglesReg1(interfacePos1) - anglesReg2(nitReg2 + interfacePos2);
    angleGradDual(i) = anglesReg1(nitReg1 + interfacePos1) - anglesReg2(interfacePos2);

    magGrad(i)       = magsReg1(interfacePos1) - magsReg2(nitReg2 + interfacePos2);
    magGradDual(i)   = magsReg1(nitReg1 + interfacePos1) - magsReg2(interfacePos2);
  if any(isnan([angleGrad(i), angleGradDual(i), magGrad(i), magGradDual(i)]))
        fprintf(' NaN detected at interface %d\n', i);
  end


        landa(i) = landa(i) + rho * angleGrad(i);
        landa(i + nit) = landa(i + nit) + rho * angleGradDual(i);
        lambdav(i) = lambdav(i) + rho * magGrad(i);
        lambdav(i + nit) = lambdav(i + nit) + rho * magGradDual(i);
        errorValues(i) = abs(angleGrad(i));
        errorValues(i + nit) = abs(angleGradDual(i));
        errorValuesv(i) = abs(magGrad(i));
        errorValuesv(i + nit) = abs(magGradDual(i));
    end


% ----------------------- Log and Check Errors -------------------------
% Track primal residuals and compute optimality gap if centralized cost is known


    maxError1 = max(errorValues);maxError2 = max(errorValuesv);maxError=max(maxError1,maxError2);
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



%% Helper function------------------------------------------------------------------------------------

function fprintf_matrix(fid, matrix)
    % Helper function to print a matrix in MATLAB .m file format
    [rows, cols] = size(matrix);
    for i = 1:rows
        fprintf(fid, '    ');
        for j = 1:cols
            fprintf(fid, '%g\t', matrix(i,j));
        end
        fprintf(fid, '\n');
    end
end
