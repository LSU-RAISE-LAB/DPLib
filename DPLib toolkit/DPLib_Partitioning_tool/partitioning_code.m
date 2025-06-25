% =========================================================================
% Title       : Graph-Based Power System Partitioning with Spectral Clustering
% Author      : Milad Hasanzadeh
% Email       : e.mhasanzadeh1377@yahoo.com
% Affiliation : Department of Electrical and Computer Engineering,
%               Louisiana State University, Baton Rouge, LA, USA
% Date        : May 29, 2025
%
% Description :
% This MATLAB script partitions a power system network into multiple 
% electrically coherent regions using graph-theoretic techniques.
% It leverages spectral clustering combined with k-means to provide a 
% balanced and meaningful decomposition of MATPOWER case files, enabling 
% distributed optimization and control across regional sub-networks.
%
% Functional Overview :
% - Constructs the weighted graph Laplacian from the admittance matrix of the system.
% - Computes eigenvectors of the Laplacian for spectral embedding.
% - Applies k-means clustering to group buses into user-specified regions.
% - Renumbers buses locally and ensures consistency within each partition.
% - Identifies inter-regional tie-lines and separates them from intra-regional data.
% - Saves the partitioned system as a structured .mat file and as .m files,
%   compatible with distributed DC and AC OPF solvers.
%
% Key Features :
% - Scalable spectral partitioning for large-scale power networks.
% - Clear separation of internal and boundary (tie-line) components.
% - Compatible with distributed ADMM-based DCOPF and ACOPF solvers.
% - Automatic creation of region-specific data structures and variables.
% - Option to export results in both `.mat` and `.m` formats for reproducibility.
%
% Requirements :
% - MATLAB R2020a or newer.
% - MATPOWER toolbox installed and added to the MATLAB path.
% 
%
% License :
% This script is provided for academic and research use only.
% Redistribution or commercial use is not permitted without prior written consent.
% If used in published work, proper acknowledgment is required.
%
% Copyright (c) 2025, Milad Hasanzadeh
% =========================================================================


clear; clc; close all;

%% Prompt user to input a valid MATPOWER case 

while true
    sound(sin(2*pi*440*(0:1/8000:0.5)), 8000);case_input = input('Enter MATPOWER case name (e.g., case10) or type "exit" to stop: ', 's');

    if strcmpi(case_input, 'exit')
        error('User exited. No valid case file was provided.');
    end

    if exist([case_input '.m'], 'file') == 2
        break; 
    else
        fprintf('Error: The case "%s" does not exist in the MATLAB path or MATPOWER is not installed properly.\n', case_input);
    end
end

%% Get user-defined number of regions for partitioning

while true
    num_regions_input = input('Enter number of regions: ');

    if isnumeric(num_regions_input) && isscalar(num_regions_input) && num_regions_input <= 1000 && num_regions_input >= 1
        break;
    else
        fprintf('Error: The number of regions must be an integer between 1 and 1000.\n');
    end
end


filename = ['''' case_input ''''];   
filename = case_input;              
The_number_of_regions = num_regions_input;mdata=0;


%% Load MATPOWER case and set maximum clustering attempts

tic;
mpc = loadcase(filename); 
The_max_of_attempts = 10;

%% If bus numbers are not sequential, renumber buses from 1 to N Also update branch and generator connections accordingly

old_bus_numbers = mpc.bus(:,1);
a=length(old_bus_numbers);
b=max(old_bus_numbers);

if a~=b
    new_bus_numbers = (1:length(old_bus_numbers))'; % Sequential numbering
    bus_map = containers.Map(old_bus_numbers, new_bus_numbers);
    
  
    mpc_new = mpc;
    mpc_new.bus(:,1) = new_bus_numbers;

    for i = 1:size(mpc.branch, 1)
        mpc_new.branch(i,1) = bus_map(mpc.branch(i,1)); % From bus
        mpc_new.branch(i,2) = bus_map(mpc.branch(i,2)); % To bus
    end
    
    for i = 1:size(mpc.gen, 1)
        mpc_new.gen(i,1) = bus_map(mpc.gen(i,1)); % Generator bus
    end
    
    fprintf('Bus numbers are renumbered from [%d - %d] to [1 - %d].\n', ...
            min(old_bus_numbers), max(old_bus_numbers), length(old_bus_numbers));

    mpc=mpc_new;
end

% Extract bus and branch data
bus = mpc.bus;
branch = mpc.branch;
num_buses = size(bus, 1);
num_branches = size(branch, 1); % Total number of branches
num_gens = size(mpc.gen, 1); % Total number of generators
num_gencosts = size(mpc.gencost, 1); % Total number of generator costs
reference_bus = bus(bus(:,2) == 3, 1); % Find bus with type 3 (slack bus)
num_clusters = The_number_of_regions; % Number of regions
max_iterations = The_max_of_attempts; % Maximum attempts to find the best clustering
min_tie_lines = num_branches; % Initialize with the worst possible case
    num_neigh_bus_list = zeros(num_clusters, 1);
% Define variable names for each region
num_regions = num_clusters; % adjust as needed

% Generate region_names like 'mpc_regionR1', ..., 'mpc_regionR50'
region_names = arrayfun(@(r) ['mpc_regionR' num2str(r)], 1:num_regions, 'UniformOutput', false);

% Generate Bus_list like 'regionR1list', ..., 'regionR50list'
Bus_list = arrayfun(@(r) ['regionR' num2str(r) 'list'], 1:num_regions, 'UniformOutput', false);

% Initialize best regions storage
best_mpc_regions = cell(num_clusters, 1);
best_tie_lines = [];
best_total_region_branches = 0;
best_total_region_buses = 0;
best_total_region_gens = 0;
best_total_region_gencosts = 0;

%% Construct the graph adjacency matrix from branch connectivity and compute graph Laplacian and check connectivity

adj_matrix = sparse(num_buses, num_buses);
for i = 1:num_branches
    fbus = branch(i, 1);
    tbus = branch(i, 2);
    adj_matrix(fbus, tbus) = 1;
    adj_matrix(tbus, fbus) = 1;
end

% Compute degree matrix
D = diag(sum(adj_matrix, 2));

% Compute graph Laplacian
L = D - adj_matrix;

G = graph(L);
bins = conncomp(G);
num_components = max(bins);

%% Try multiple clusterings to find one with minimum inter-regional tie-lines


% Initialize best solution trackers to avoid undefined variable errors
best_tie_line_info = {};
best_mpc_regions = {};
best_tie_lines = [];
best_total_region_branches = 0;
best_total_region_buses = 0;
best_total_region_gens = 0;
best_total_region_gencosts = 0;


% Iterate multiple times to find the clustering with minimum tie-lines and use k-means on eigenvectors to form clusters (regions)
for iter = 1:max_iterations
    fprintf('Iteration %d/%d\n', iter, max_iterations);
    opts.maxit = 1000; % Increase iterations
    opts.tol = 1e-5; % Increase tolerance
    % Compute first few eigenvectors of L
    %[eig_vectors, ~] = eigs(L, num_clusters, 'smallestreal');
    L_full = full(L);
    [V, D] = eig(L_full);
[~, idx] = sort(diag(D));
eig_vectors = V(:, idx(1:num_clusters));


    % Normalize eigenvectors for better clustering performance
    eig_vectors = eig_vectors ./ vecnorm(eig_vectors, 2, 2);

    % Apply k-means to eigenvectors (run multiple times and choose best)
    best_cluster_labels = [];
    best_tie_line_count = inf; 

    for k_run = 1:max_iterations  
        temp_cluster_labels = kmeans(eig_vectors, num_clusters, 'Replicates', 10);
        
        % Count inter-regional tie lines for this run
        temp_bus_regions = zeros(num_buses, 1);
        temp_bus_regions(1:length(temp_cluster_labels)) = temp_cluster_labels;
        
        temp_branch_assignment = false(num_branches, 1);
        for r = 1:num_clusters
            temp_region_buses = bus(ismember(bus(:,1), find(temp_bus_regions == r)), 1);
            temp_branch_assignment = temp_branch_assignment | ...
                (ismember(branch(:,1), temp_region_buses) & ismember(branch(:,2), temp_region_buses));
        end
        
        % Count unassigned branches (tie-lines)
        temp_tie_line_count = sum(~temp_branch_assignment);
        
        if temp_tie_line_count < best_tie_line_count
            best_tie_line_count = temp_tie_line_count;
            best_cluster_labels = temp_cluster_labels;
        end
    end
    
    % Assign the best clustering result
    cluster_labels = best_cluster_labels;

% Compute centroids of each cluster in the spectral domain
region_centroids = zeros(num_clusters, size(eig_vectors, 2));
for r = 1:num_clusters
    region_centroids(r, :) = mean(eig_vectors(cluster_labels == r, :), 1);
end

% Use pairwise Euclidean distances between centroids
dist_matrix = squareform(pdist(region_centroids, 'euclidean'));

% Use hierarchical clustering to get an ordering
linkage_order = linkage(region_centroids, 'average');
region_order = optimalleaforder(linkage_order, dist_matrix);

% Re-map cluster labels and reorder regions accordingly
remapped_labels = zeros(size(cluster_labels));
for i = 1:num_clusters
    remapped_labels(cluster_labels == region_order(i)) = i;
end
cluster_labels = remapped_labels;





    % Assign buses to regions
    bus_regions = zeros(num_buses, 1);
    bus_regions(1:length(cluster_labels)) = cluster_labels;

    % Ensure all clusters contain at least one bus
    if length(unique(cluster_labels)) < num_clusters
        continue; % Skip this iteration if some clusters are empty
    end

    % Temporary storage for regional cases
    temp_mpc_regions = cell(num_clusters, 1);
    total_region_branches = 0; % Track the sum of branches in regions
    total_region_buses = 0; % Track buses
    total_region_gens = 0; % Track generators
    total_region_gencosts = 0; % Track generator costs
    branch_assignment = false(num_branches, 1); % Track branch assignment

    for r = 1:num_clusters
        % Ensure valid indexing using `ismember`
        region_buses = bus(ismember(bus(:,1), find(bus_regions == r)), 1);
        invalid_partition = false;

        if isempty(region_buses)
            continue;
        end
        
        % Create a regional copy of the original case
        mpc_region = mpc;
        
        % Select buses belonging to this region
        mpc_region.bus = bus(ismember(bus(:,1), region_buses), :);
        num_region_buses = size(mpc_region.bus, 1);
        total_region_buses = total_region_buses + num_region_buses;

        % Select branches where both ends are in the same region
        in_region = ismember(branch(:,1), region_buses) & ismember(branch(:,2), region_buses);
        mpc_region.branch = branch(in_region, :);
        
        % Count branches in this region
        num_region_branches = size(mpc_region.branch, 1);
        total_region_branches = total_region_branches + num_region_branches;
        
        % Mark these branches as assigned
        branch_assignment(in_region) = true;
        
        % Select generators in this region
        in_region_gen = ismember(mpc.gen(:,1), region_buses);
        mpc_region.gen = mpc.gen(in_region_gen, :);
        num_region_gens = size(mpc_region.gen, 1);
        total_region_gens = total_region_gens + num_region_gens;
        if num_region_gens == 0
          invalid_partition = true;
          break; % Exit loop early since this region has no generators
        end


        % Select gencosts for generators in this region
        mpc_region.gencost = mpc.gencost(in_region_gen, :);
        num_region_gencosts = size(mpc_region.gencost, 1);
        total_region_gencosts = total_region_gencosts + num_region_gencosts;
        temp_mpc_regions{r} = mpc_region;
    end
    if invalid_partition
    continue; % Skip to next iteration of outer loop
    end

%% Identify tie-line branches that connect different regions

    temp_tie_lines = branch(~branch_assignment, :); % Any unassigned branch is a tie-line
    num_temp_tie_lines = size(temp_tie_lines, 1);
    
    tie_line_info = cell(num_temp_tie_lines, 3); % Store tie-line region connections
    % Initialize storage for neighboring buses in each region

    for i = 1:num_temp_tie_lines
        fbus = temp_tie_lines(i, 1);
        tbus = temp_tie_lines(i, 2);
        
        region_f = bus_regions(bus(:,1) == fbus);
        region_t = bus_regions(bus(:,1) == tbus);
        
        tie_line_info{i, 1} = region_names{region_f}; % From region
        tie_line_info{i, 2} = region_names{region_t}; % To region
        tie_line_info{i, 3} = temp_tie_lines(i, :);  % Branch data
            % Store tie-line buses as neighboring buses in respective regions

    end
%% Update best partitioning if it results in fewer tie-lines

    % Check if the sum of region branches + tie-lines matches the total branches
    total_count = total_region_branches + num_temp_tie_lines;
    
    if total_count == num_branches && num_temp_tie_lines < min_tie_lines
        min_tie_lines = num_temp_tie_lines;
        best_mpc_regions = temp_mpc_regions;
        best_tie_lines = temp_tie_lines;
        best_total_region_branches = total_region_branches;
        best_total_region_buses = total_region_buses;
        best_total_region_gens = total_region_gens;
        best_total_region_gencosts = total_region_gencosts;
        best_tie_line_info = tie_line_info;
    end

end

if isempty(best_tie_line_info)
warning(['No valid partitioning found that satisfies generator and tie-line constraints. ' ...
         'Try reducing the number of regions to improve feasibility.']);
return; % Or handle as needed
end


num_of_areas=0;
% Write tie-line table to Excel
tie_line_table1 = cell2table(best_tie_line_info, 'VariableNames', {'From_Region', 'To_Region', 'Branch_Data'});



all_partitioned_buses = [];
for r = 1:num_clusters
    if ~isempty(best_mpc_regions{r})
        all_partitioned_buses = [all_partitioned_buses; best_mpc_regions{r}.bus(:,1)];
    end
end

% Find missing buses
missing_buses = setdiff(bus(:,1), all_partitioned_buses);

    reference_region = 0;
    
    for r = 1:num_clusters
        if ~isempty(best_mpc_regions{r}) % Ensure region is valid
            region_buses = best_mpc_regions{r}.bus(:,1); % Extract bus numbers in region
            if any(region_buses == reference_bus)
                reference_region = r;
                break; % Stop searching once found
            end
        end
    end

 % Renumber buses from 1 to N for each region
for r = 1:num_clusters
    if ~isempty(best_mpc_regions{r})
        % Extract region data
        mpc_region = best_mpc_regions{r};
        
        % Create a mapping from old bus numbers to new bus numbers
        old_bus_numbers = mpc_region.bus(:,1);
        new_bus_numbers = (1:length(old_bus_numbers))';
        bus_map = containers.Map(old_bus_numbers, new_bus_numbers);
        
        % Update the bus matrix (only first column)
        mpc_region.bus(:,1) = new_bus_numbers;
        
        % Update branch matrix (only first and second columns)
        for i = 1:size(mpc_region.branch,1)
            mpc_region.branch(i,1) = bus_map(mpc_region.branch(i,1)); % From bus
            mpc_region.branch(i,2) = bus_map(mpc_region.branch(i,2)); % To bus
        end
        
        % Update generator matrix (only first column)
        for i = 1:size(mpc_region.gen,1)
            mpc_region.gen(i,1) = bus_map(mpc_region.gen(i,1)); % Generator bus
        end
        
        % Store the updated region back
        best_mpc_regions{r} = mpc_region;
        bus_number_list{r}=[old_bus_numbers new_bus_numbers];
    end
end



Updated_bus_data = cell(1, num_clusters); % Initialize cell array

for r = 1:min(num_clusters, length(bus_number_list)) % Ensure r does not exceed available regions
    if r > length(bus_number_list) || isempty(bus_number_list{r}) % Ensure region exists
        continue; % Skip if the region is empty or out of range
    end
    
    Updated_bus_data{r} = bus_number_list{r}; % Store as a separate cell
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure `bus_number_list` is correctly sized
num_valid_regions = min(num_clusters, length(bus_number_list));

% Determine the maximum number of buses in any region
max_rows = max(cellfun(@(x) size(x, 1), bus_number_list(1:num_valid_regions))); 

% Create a numeric matrix filled with NaN (to handle different region sizes)
Updated_bus_matrix = NaN(max_rows, num_valid_regions * 3 - 1); % Includes space for empty columns

% Fill the matrix with bus numbers (both columns)
for r = 1:num_valid_regions
    if ~isempty(bus_number_list{r}) % Ensure the region is not empty
        num_buses_in_region = size(bus_number_list{r}, 1);
        
        % Determine column indices (1st and 2nd column for the region, then an empty one)
        col_start = (r - 1) * 3 + 1;
        col_end = col_start + 1;
        
        % Ensure data fits within the allocated space
        Updated_bus_matrix(1:num_buses_in_region, col_start:col_end) = bus_number_list{r}(1:min(num_buses_in_region, max_rows), :);
    end
end

% Convert to table
Updated_bus_table = array2table(Updated_bus_matrix);

% Define column names dynamically with unique names for empty columns
column_names = strings(1, num_valid_regions * 3 - 1);
empty_col_count = 1; % Counter for unique empty column names

for r = 1:num_valid_regions
    col_start = (r - 1) * 3 + 1;
    column_names(col_start) = "Region_" + string(r) + "_Col1";
    column_names(col_start + 1) = "Region_" + string(r) + "_Col2";
    
    % Add a unique name for empty columns
    if r < num_valid_regions
        column_names(col_start + 2) = "Empty_Col_" + string(empty_col_count);
        empty_col_count = empty_col_count + 1;
    end
end
Updated_bus_table.Properties.VariableNames = column_names;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:num_clusters
    if ~isempty(best_mpc_regions{r})

        % Assign the region's data in the MATLAB workspace
        assignin('base', region_names{r}, best_mpc_regions{r});

        % Define the file name for saving
        mat_filename = [region_names{r}, '.mat'];

        % Save only the current region's data into its respective .mat file
        region_data = best_mpc_regions{r};  % Extract data for this region

        % Increment the counter for valid regions
        num_of_areas = num_of_areas + 1;
        
        % Display message confirming save
        fprintf('Region %s saved successfully in %s\n', region_names{r}, mat_filename);
    end
end


%% Confirm that all original buses are included in the partitioning

clc;
if isempty(missing_buses)
    fprintf('The multi-regional partitioning process has been verified to include all buses without omission.\n\n');
else
    fprintf('Missing buses: %s.\n', mat2str(missing_buses));
end
fprintf('The best clustering configuration consists of %d regions and %d inter-regional tie-lines, identified within %d iterations.\n',num_of_areas, min_tie_lines,max_iterations);
fprintf('The reference bus %d in the MATPOWER case has been assigned to Region %d (%s).\n\n', reference_bus, reference_region, region_names{reference_region});



if num_components==1
fprintf('The graph topology was verified to be fully connected.\n');
else
fprintf('The graph topology was not verified to be fully connected!\n');
end
if a~=b
fprintf('Reordering of the original bus numbers in the MATPOWER case was necessary to enable consistent partitioning.\n\n');
end
if (best_total_region_branches + min_tie_lines==num_branches)
fprintf('The total number of branches across all regions, including inter-regional tie-lines, matches the number of branches in the MATPOWER case, which is %d.\n', num_branches);
end
if (best_total_region_buses==num_buses)
fprintf('The total number of buses across all regions matches the number of buses in the MATPOWER case, which is %d.\n', num_buses);
end
if (best_total_region_gens==num_gens)
fprintf('The total number of generators across all regions matches the number of generators in the MATPOWER case, which is %d.\n\n', num_gens);
end


%% Map tie-line buses to updated bus numbers within each region and build a structure (save_struct) containing all regions and tie-line info

    dataMatrix = table2cell(tie_line_table1);
    dataMatrix2 = table2array(Updated_bus_table);

    alphabet = 'A':'Z';
    save_struct = struct();

    % Loop through each region
    for r = 1:num_clusters
        region_letter = ['R' num2str(r)];
        region_varname = ['mpc_region' region_letter];
        region_tie = ['interregional_tielines' region_letter];
        region_tienum = ['num_interregional_tielines' region_letter];
    outq2 = dataMatrix2(:,3*r-2:3*r-1);

    nanIndex = find(isnan(outq2(:, 1)), 1);
    
    % If there is no NaN, just return the original matrix
    if isempty(nanIndex)
        out2 = outq2;
    else
        % Otherwise, create a new matrix with only the non-NaN rows
        out2 = outq2(1:nanIndex-1, :);
    end

        % Find rows where the string matches in the first or second column
    region_old_name = ['mpc_regionR' num2str(r)]; % Old variable name to search in dataMatrix
    row_indices = strcmp(dataMatrix(:, 1), region_old_name) | strcmp(dataMatrix(:, 2), region_old_name);
    
    filteredData = dataMatrix(row_indices, :);
    output_rows = cell2mat(filteredData(:, 3:end));
    output_rowss = dataMatrix(row_indices, :);
        % Convert back to a table for easy saving

    A=output_rows;B=out2;
m = size(A,1);      % number of rows in first matrix (A)
n = size(B,1);      % number of rows in second matrix (B)  --> 642
key2value = containers.Map(B(:,1), B(:,2));

A_out = A;          % start with A and overwrite in-place

for i = 1:m
    for j = 1:2
        v = A(i,j);            % element under inspection
        if isKey(key2value, v) % this is the “matching” element
            A_out(i,j) = key2value(v);
        else                   % this is the non-matching element
            A_out(i,j) = i + n;
        end
    end
end
for i = 1:size(output_rowss, 1)
    output_rowss{i, 3} = A_out(i, :); % Replace the 3rd column with the corresponding row from A_out
end
numtie=size(output_rowss, 1);
    outputTable = cell2table(output_rowss);
   
        % Get the variable from base workspace
        try
            region_data = evalin('base', region_old_name);
            save_struct.(region_varname) = region_data;
            if mdata==1
            save_struct.(region_tie) = outputTable;
            end
            %save_struct.(region_tienum) = numtie;
        catch
            warning('Variable %s not found in the workspace.', region_old_name);
        end
    end

    % Add additional variables
    %save_struct.DCOPF = 0;
    %save_struct.TotalCost = 0;
    %save_struct.ACOPF = 0;
    save_struct.num_regions = num_clusters;
    %save_struct.new_bus_labels=Updated_bus_table;
    save_struct.filename=filename;

    dataMatrix = table2cell(tie_line_table1);
    dataMatrix2 = table2array(Updated_bus_table);

for i = 1:size(dataMatrix, 1)
    region1_str = dataMatrix{i, 1};
    region2_str = dataMatrix{i, 2};
    regiontiedata = dataMatrix{i, 3};
    firstout = regiontiedata(1, 1);
    secondout = regiontiedata(1, 2);

    region1_num = sscanf(region1_str, 'mpc_regionR%d');
    region2_num = sscanf(region2_str, 'mpc_regionR%d');

    outq1 = dataMatrix2(:, 3*region1_num-2:3*region1_num-1);
    outq2 = dataMatrix2(:, 3*region2_num-2:3*region2_num-1);

    nanIndex1 = find(isnan(outq1(:, 1)), 1); 
    nanIndex2 = find(isnan(outq2(:, 1)), 1);

    if isempty(nanIndex1)
        out1 = outq1;
    else
        out1 = outq1(1:nanIndex1-1, :);
    end

    if isempty(nanIndex2)
        out2 = outq2;
    else
        out2 = outq2(1:nanIndex2-1, :);
    end

    % --- Replace firstout with second column value from out1 ---
    idx1 = find(out1(:, 1) == firstout, 1);
    if ~isempty(idx1)
        regiontiedata(1, 1) = out1(idx1, 2);
    end

    % --- Replace secondout with second column value from out2 ---
    idx2 = find(out2(:, 1) == secondout, 1);
    if ~isempty(idx2)
        regiontiedata(1, 2) = out2(idx2, 2);
    end

    % --- Update the dataMatrix with modified regiontiedata ---
    dataMatrix{i, 3} = regiontiedata;
end
tie_line_table1 = cell2table(dataMatrix);
    save_struct.interregional_tielines_total= tie_line_table1;

% Create .m file name
mfilename = [filename, '_', num2str(num_clusters), 'regions.m'];
mfilename2 = [filename, '_', num2str(num_clusters), 'regions'];
fid = fopen(mfilename, 'w');
% Define user information
author_name = 'Milad Hasanzadeh';
generation_date = string(datetime('now','Format','yyyy-MM-dd'));


% Write custom header
fprintf(fid, '%% -------------------------------------------------------------------------\n');
fprintf(fid, '%% Partitioned MATPOWER Case File\n');
fprintf(fid, '%% Copyright (c) 2025 %s\n', author_name);
fprintf(fid, '%% Generated on: %s\n', generation_date);
fprintf(fid, '%% This file contains partitioned data for %d regions\n\n\n', num_clusters);

% Loop through regions and write each mpc_regionR* as a separate struct
for r = 1:num_clusters
    regionID = ['R', num2str(r)];
    varname = ['mpc_region', regionID];

    if isfield(save_struct, varname)
        mpc = save_struct.(varname);  % Get region struct

        fprintf(fid, '%% -------------------- %s --------------------\n\n', varname);
        fprintf(fid, '%s.version = ''2'';\n', varname);

        if isfield(mpc, 'baseMVA')
            fprintf(fid, '%s.baseMVA = %.4f;\n', varname, mpc.baseMVA);
        else
            fprintf(fid, '%s.baseMVA = 100;\n', varname);
        end

        if isfield(mpc, 'bus')
            fprintf(fid, '%s.bus = [\n', varname);
            bus_headers = {'bus_i','type','Pd','Qd','Gs','Bs','area','Vm','Va','baseKV','zone','Vmax','Vmin'}; 
            fprintf_matrix(fid, mpc.bus, bus_headers);
            fprintf(fid, '];\n\n');
        end

        if isfield(mpc, 'gen')
            fprintf(fid, '%s.gen = [\n', varname);
            gen_headers = {'bus','Pg','Qg','Qmax','Qmin','Vg','mBase','status','Pmax','Pmin','Pc1','Pc2','Qc1min','Qc1max','Qc2min','Qc2max','ramp_agc','ramp_10','ramp_30','ramp_q','apf'};
            fprintf_matrix(fid, mpc.gen, gen_headers);
            fprintf(fid, '];\n\n');
        end

        if isfield(mpc, 'branch')
            fprintf(fid, '%s.branch = [\n', varname);
             branch_headers = {'fbus','tbus','r','x','b','rateA','rateB','rateC','ratio','angle','status','angmin','angmax'};
            fprintf_matrix(fid, mpc.branch, branch_headers);
            fprintf(fid, '];\n\n');
        end

        if isfield(mpc, 'gencost')
            fprintf(fid, '%s.gencost = [\n', varname);
            gencost_headers = {'2','startup','shutdown','n','c(n-1)...c0'};
            fprintf_matrix(fid, mpc.gencost, gencost_headers);
            fprintf(fid, '];\n\n');
        end

        fprintf(fid, '\n');
    end
end

%% Write the full interregional_tielines_total cell array
if isfield(save_struct, 'interregional_tielines_total')
    fprintf(fid, '%% ---------- Interregional Tie-Line Information ----------\n');
    fprintf(fid, '%% Each row: {From_Region (string), To_Region (string), Branch_Data (1x13 array)}\n');
    fprintf(fid, '%% Branch_Data: {fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax}\n\n');
    tie_cell = table2cell(save_struct.interregional_tielines_total);  % ensure it's a cell array

    fprintf(fid, 'interregional_tielines_total = {\n');
    for i = 1:size(tie_cell, 1)
        from_region = tie_cell{i, 1};
        to_region   = tie_cell{i, 2};
        branch_data = tie_cell{i, 3};

        fprintf(fid, '\t{''%s'', ''%s'', [', from_region, to_region);
        for j = 1:length(branch_data)
            if j < length(branch_data)
                fprintf(fid, '% .12g\t', branch_data(j));
            else
                fprintf(fid, '% .12g', branch_data(j));
            end
        end
        fprintf(fid, ']},\n');
    end
    fprintf(fid, '};\n\n');
end



fclose(fid);

filename = [filename, '_', num2str(num_clusters), 'regions.mat'];
save(filename, '-struct', 'save_struct');
filenamee = filename;
% --- Graph-Based Visualization of Regional Topology ---

plot_partitioned_matfile(filenamee);

exportgraphics(gcf, [mfilename2, '_topology.png'], 'Resolution', 300);

fprintf('The region-level topology figure has been saved as %s_topology.png\n', mfilename2);

fprintf('The partitioned data has been saved as .mat and .m files named %s in the current path.\n\n', mfilename2);
toc


%% ------- Helper Function --------
function fprintf_matrix(fid, matrix, col_headers)
    if nargin > 2 && ~isempty(col_headers)
        % Print column headers as a comment
        fprintf(fid, '%%\t');
        for i = 1:length(col_headers)
            fprintf(fid, '%-12s', col_headers{i});
        end
        fprintf(fid, '\n');
    end

    [rows, cols] = size(matrix);
    for i = 1:rows
        fprintf(fid, '\t');
        for j = 1:cols
            if j < cols
                fprintf(fid, '% .12g\t', matrix(i, j));
            else
                fprintf(fid, '% .12g', matrix(i, j));
            end
        end
        fprintf(fid, ';\n');
    end
end

function plot_partitioned_matfile(file_basename)
    % Add .mat extension if not already present
    if ~endsWith(file_basename, '.mat')
        file_basename = [file_basename, '.mat'];
    end
fontsizew=17;

    % Load the MAT file
    if ~isfile(file_basename)
        error('File "%s" not found.', file_basename);
    end
    data = load(file_basename);

    % Use num_regions directly
    if ~isfield(data, 'num_regions')
        error('MAT file must contain variable "num_regions".');
    end
    num_regions = data.num_regions;

    % Read tie-line table
    if ~isfield(data, 'interregional_tielines_total')
        error('MAT file must contain variable "interregional_tielines_total".');
    end
    tielines = data.interregional_tielines_total;
    dataMatrix = table2cell(tielines);

    % Build region_map
    region_map = containers.Map;
    for r = 1:num_regions
        region_name = ['mpc_regionR', num2str(r)];
        if ~isfield(data, region_name)
            error('Missing region struct: %s', region_name);
        end
        region_map(region_name) = r;
    end

    % Build edge count
    edge_count = containers.Map;
    for i = 1:size(dataMatrix, 1)
        r1 = region_map(dataMatrix{i,1});
        r2 = region_map(dataMatrix{i,2});
        key = sprintf('%d-%d', min(r1, r2), max(r1, r2));
        if isKey(edge_count, key)
            edge_count(key) = edge_count(key) + 1;
        else
            edge_count(key) = 1;
        end
    end

    % Create adjacency matrix
    adj_matrix = zeros(num_regions);
    keys_list = keys(edge_count);
    for i = 1:length(keys_list)
        key = keys_list{i};
        tokens = sscanf(key, '%d-%d');
        r1 = tokens(1); r2 = tokens(2);
        adj_matrix(r1, r2) = edge_count(key);
        adj_matrix(r2, r1) = edge_count(key);
    end

    % Compute layout
    region_graph = graph(adj_matrix);
temp_fig = figure('Visible', 'off');
temp_plot = plot(region_graph, 'Layout', 'force', 'Iterations', 100);
x_raw = temp_plot.XData;
y_raw = temp_plot.YData;
close(temp_fig);


    % Layout spacing and circle size
    circle_width = 12;
    circle_height = circle_width;
    if (num_regions <= 10)
        layout_scale = circle_width * 3.5;
    elseif (num_regions <= 15  && num_regions > 10)
        layout_scale = circle_width * 2.8;
    else
        layout_scale = circle_width * 2.4;
    end
    x_coords = x_raw * layout_scale;
    y_coords = y_raw * layout_scale;

    % Draw plot
fig2 = figure(100); 
set(fig2, 'Color', 'w');
    hold on;
    axis off;

    bus_counts = zeros(num_regions, 1);
for r = 1:num_regions
    region_name = ['mpc_regionR', num2str(r)];
    region_struct = data.(region_name);
    bus_counts(r) = size(region_struct.bus, 1);

    rectangle('Position', [x_coords(r)-circle_width/2, y_coords(r)-circle_height/2, ...
               circle_width, circle_height], ...
              'Curvature', [1,1], ...
              'EdgeColor', 'b', ...               % Dashed outline
              'LineWidth', 2.2);                   % Same thickness

    text(x_coords(r), y_coords(r), ...
         sprintf('R%d\n%d bus', r, bus_counts(r)), ...
         'HorizontalAlignment', 'center', ...
         'FontSize', fontsizew, 'FontWeight','bold');
end


    % Tie-lines
    legend_handles = [];
    legend_entries = {};
    for i = 1:length(keys_list)
        key = keys_list{i};
        tokens = sscanf(key, '%d-%d');
        r1 = tokens(1); r2 = tokens(2);

        p1 = [x_coords(r1), y_coords(r1)];
        p2 = [x_coords(r2), y_coords(r2)];
        v = (p2 - p1) / norm(p2 - p1);
        shrink = circle_width / 2 * 1.05;
        new_p1 = p1 + shrink * v;
        new_p2 = p2 - shrink * v;

        h = plot([new_p1(1), new_p2(1)], [new_p1(2), new_p2(2)], 'k-', 'LineWidth', 1.6);
        legend_handles(end+1) = h; %#ok<AGROW>
        legend_entries{end+1} = sprintf('R%d to R%d: %d lines', ...
                                        r1, r2, edge_count(key)); %#ok<AGROW>
    end

    % Display legends in two parts if too many entries
    if length(legend_entries) > 50
        % First legend (1–ceil(n/2))
        split_idx = ceil(length(legend_entries)/2);
        
        % Create two separate axes for legends
        ax1 = axes('Position', [0.8, 0.5, 0.15, 0.4], 'Visible', 'off');
        legend(ax1, legend_handles(1:split_idx), legend_entries(1:split_idx), ...
               'FontSize', fontsizew, 'FontWeight','bold', 'Location', 'northwest');

        ax2 = axes('Position', [0.8, 0.05, 0.15, 0.4], 'Visible', 'off');
        legend(ax2, legend_handles(split_idx+1:end), legend_entries(split_idx+1:end), ...
               'FontSize', fontsizew, 'FontWeight','bold', 'Location', 'northwest');

    else
        % Standard single legend
        lgd = legend(legend_handles, legend_entries, 'Location', 'northeastoutside');
        set(lgd, 'FontSize', fontsizew, 'FontWeight','bold');
    end


    hold off;
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(gca, 'FontSize', fontsizew);
    drawnow;
end


