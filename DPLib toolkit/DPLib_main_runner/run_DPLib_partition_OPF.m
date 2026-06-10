function [DPLibMainSummary, all_results] = run_DPLib_partition_OPF(case_list, mode, varargin)
%RUN_DPLIB_PARTITION_OPF
% -------------------------------------------------------------------------
% Purpose:
%   Function version of the DPLib main runner.
%
%   This function automatically:
%       1) partitions centralized PGLib/MATPOWER cases,
%       2) solves distributed DC OPF if requested,
%       3) solves distributed AC OPF if requested,
%       4) saves error/gap/rho figures,
%       5) returns and saves a summary table.
%
% Inputs:
%   case_list:
%       Cell array with two columns:
%           column 1: centralized case name
%           column 2: number of regions
%
%       Example:
%           case_list = {
%               'pglib_opf_case4917_goc', 10
%           };
%
%   mode:
%       'partition'  -> only partition
%       'dc'         -> partition + centralized DC + distributed DC
%       'ac'         -> partition + centralized AC + distributed AC
%       'both'       -> partition + DC + AC
%       'none'       -> same as partition
%
% Optional name-value inputs:
%   'force_partition'       : default true
%   'use_weighted'          : default false
%   'bus_region_map'        : default []
%
%   'thresholdDC'           : default 1e-3
%   'thresholdAC'           : default 1e-3
%   'maxIterDC'             : default 10000
%   'maxIterAC'             : default 10000
%   'rho_dc'                : default 1
%   'rho_ac'                : default 1
%   'holds'                 : default 1e6
%
%   'centralizedCostDC'     : default []
%   'centralizedCostAC'     : default []
%
%   'saveFigures'           : default true
%   'closeExtraFigures'     : default true
%   'saveSummary'           : default true
%   'summaryName'           : default 'DPLib_main_runner_summary'
%
% Outputs:
%   DPLibMainSummary:
%       Summary table of partitioning and OPF results.
%
%   all_results:
%       Cell array containing detailed outputs for each case.
%
% Example:
%   case_list = {
%       'pglib_opf_case4917_goc', 10
%   };
%
%   [summary, results] = run_DPLib_partition_OPF(case_list, 'both');
%
% Example with custom settings:
%   [summary, results] = run_DPLib_partition_OPF(case_list, 'both', ...
%       'rho_dc', 1, ...
%       'rho_ac', 1, ...
%       'thresholdDC', 1e-3, ...
%       'thresholdAC', 1e-3, ...
%       'maxIterDC', 10000, ...
%       'maxIterAC', 10000, ...
%       'holds', 1e6);

    %% ============================================================
    %  Basic checks
    % ============================================================
    if nargin < 1 || isempty(case_list)
        error('case_list is required.');
    end

    if nargin < 2 || isempty(mode)
        mode = 'both';
    end

    mode = lower(string(mode));

    valid_modes = ["partition", "none", "dc", "ac", "both"];

    if ~ismember(mode, valid_modes)
        error('Invalid mode. Use: partition, none, dc, ac, or both.');
    end

    %% ============================================================
    %  Parse optional inputs
    % ============================================================
    p = inputParser;

    addParameter(p, 'force_partition', true);
    addParameter(p, 'use_weighted', false);
    addParameter(p, 'bus_region_map', []);

    addParameter(p, 'thresholdDC', 1e-4);
    addParameter(p, 'thresholdAC', 1e-4);
    addParameter(p, 'maxIterDC', 3000);
    addParameter(p, 'maxIterAC', 3000);
    addParameter(p, 'rho_dc', 5000);
    addParameter(p, 'rho_ac', 5000);
    addParameter(p, 'holds', 1e8);

    addParameter(p, 'centralizedCostDC', []);
    addParameter(p, 'centralizedCostAC', []);

    addParameter(p, 'saveFigures', true);
    addParameter(p, 'closeExtraFigures', true);
    addParameter(p, 'saveSummary', true);
    addParameter(p, 'summaryName', 'DPLib_main_runner_summary');

    parse(p, varargin{:});
    opt = p.Results;

    %% ============================================================
    %  MATLAB display and figure settings
    % ============================================================
    format long g;

    set(groot, 'DefaultAxesFontSize', 9);
    set(groot, 'DefaultTextFontSize', 9);
    set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.0);
    set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.0);

    %% ============================================================
    %  Check required functions
    % ============================================================
    if exist('partitioning_code', 'file') ~= 2
        error('partitioning_code.m was not found. Add the DPLib partitioning folder to MATLAB path.');
    end

    if any(mode == ["dc", "both"]) && exist('DistDCOPF', 'file') ~= 2
        error('DistDCOPF.m was not found. Add the DPLib DC OPF solver folder to MATLAB path.');
    end

    if any(mode == ["ac", "both"]) && exist('DistACOPF', 'file') ~= 2
        error('DistACOPF.m was not found. Add the DPLib AC OPF solver folder to MATLAB path.');
    end

    %% ============================================================
    %  Preallocate summary fields
    % ============================================================
    num_cases = size(case_list, 1);

    CaseName = strings(num_cases, 1);
    PartitionedCaseName = strings(num_cases, 1);
    Regions = zeros(num_cases, 1);

    TieLines = NaN(num_cases, 1);
    RefBusRegion = strings(num_cases, 1);

    TotalBuses = NaN(num_cases, 1);
    TotalBranchesOriginal = NaN(num_cases, 1);
    TotalBranchesRegional = NaN(num_cases, 1);
    TotalGenerators = NaN(num_cases, 1);

    PartitionTimeSec = NaN(num_cases, 1);

    CentralizedDCCost = NaN(num_cases, 1);
    CentralizedACCost = NaN(num_cases, 1);

    RuntimeDC = NaN(num_cases, 1);
    RuntimeAC = NaN(num_cases, 1);

    StatusPartition = strings(num_cases, 1);
    StatusDC = strings(num_cases, 1);
    StatusAC = strings(num_cases, 1);
    Message = strings(num_cases, 1);

    all_results = cell(num_cases, 1);

    %% ============================================================
    %  Main loop
    % ============================================================
    for c = 1:num_cases

        case_name = case_list{c, 1};
        k_regions = case_list{c, 2};

        partitioned_case_name = sprintf('%s_%dregions', case_name, k_regions);
        generated_mat_file = sprintf('%s.mat', partitioned_case_name);

        CaseName(c) = string(case_name);
        PartitionedCaseName(c) = string(partitioned_case_name);
        Regions(c) = k_regions;

        fprintf('\n============================================================\n');
        fprintf('Running case %d/%d\n', c, num_cases);
        fprintf('Centralized case : %s\n', case_name);
        fprintf('Regions          : %d\n', k_regions);
        fprintf('Partitioned case : %s\n', partitioned_case_name);
        fprintf('Mode             : %s\n', mode);
        fprintf('============================================================\n');

        result_case = struct();
        result_case.case_name = case_name;
        result_case.k_regions = k_regions;
        result_case.partitioned_case_name = partitioned_case_name;

        try
            %% ------------------------------------------------------------
            %  Step 1: Partitioning
            % -------------------------------------------------------------
            if opt.force_partition || ~isfile(generated_mat_file)

                fprintf('\n------------------------------------------------------------\n');
                fprintf('[1] Running partitioning_code...\n');
                fprintf('------------------------------------------------------------\n');

                tic_partition = tic;

                partition_result = partitioning_code( ...
                    case_name, ...
                    k_regions, ...
                    opt.use_weighted, ...
                    opt.bus_region_map);

                partition_time = toc(tic_partition);

                result_case.partition_result = partition_result;
                PartitionTimeSec(c) = partition_time;
                StatusPartition(c) = "Success";

                fprintf('> Partitioning finished in %.4f seconds.\n', partition_time);

            else

                fprintf('\n------------------------------------------------------------\n');
                fprintf('[1] Partitioned file already exists. Skipping partitioning.\n');
                fprintf('------------------------------------------------------------\n');

                partition_result = struct();
                partition_result.message = 'Existing partitioned file used.';

                result_case.partition_result = partition_result;
                StatusPartition(c) = "Skipped";
                PartitionTimeSec(c) = 0;

            end

            if ~isfile(generated_mat_file)
                error('Generated .mat file was not found: %s', generated_mat_file);
            end

            fprintf('> Using partitioned file: %s\n', generated_mat_file);

            %% ------------------------------------------------------------
            %  Step 2: Extract partitioning information
            % -------------------------------------------------------------
            if isfield(partition_result, 'min_tie_lines')
                TieLines(c) = partition_result.min_tie_lines;
            end

            if isfield(partition_result, 'slack_region')
                RefBusRegion(c) = string(partition_result.slack_region);
            elseif isfield(partition_result, 'reference_bus_region')
                RefBusRegion(c) = string(partition_result.reference_bus_region);
            else
                RefBusRegion(c) = "";
            end

            %% ------------------------------------------------------------
            %  Step 3: Load generated distributed data
            % -------------------------------------------------------------
            generated_data = load(generated_mat_file);

            if isfield(generated_data, 'filename')
                centralized_case_from_file = generated_data.filename;
            else
                centralized_case_from_file = case_name;
            end

            fprintf('> Centralized case inside generated file: %s\n', centralized_case_from_file);

            if isfield(generated_data, 'num_regions')
                k_loaded = generated_data.num_regions;
            else
                k_loaded = k_regions;
            end

            total_region_buses = 0;
            total_region_branches = 0;
            total_region_gens = 0;

            for r = 1:k_loaded

                region_field = sprintf('mpc_regionR%d', r);

                if isfield(generated_data, region_field)

                    region_data = generated_data.(region_field);

                    if isfield(region_data, 'bus')
                        total_region_buses = total_region_buses + size(region_data.bus, 1);
                    end

                    if isfield(region_data, 'branch')
                        total_region_branches = total_region_branches + size(region_data.branch, 1);
                    end

                    if isfield(region_data, 'gen')
                        total_region_gens = total_region_gens + size(region_data.gen, 1);
                    end
                end
            end

            TotalBuses(c) = total_region_buses;
            TotalBranchesRegional(c) = total_region_branches;
            TotalGenerators(c) = total_region_gens;

            %% ------------------------------------------------------------
            %  Step 4: Original centralized branch count
            % -------------------------------------------------------------
            mpc_original = loadcase(case_name);
            TotalBranchesOriginal(c) = size(mpc_original.branch, 1);

            %% ------------------------------------------------------------
            %  If only partitioning is requested
            % -------------------------------------------------------------
            if mode == "partition" || mode == "none"

                fprintf('\nMode is "%s". OPF solvers are skipped.\n', mode);

                StatusDC(c) = "Skipped";
                StatusAC(c) = "Skipped";

                all_results{c} = result_case;
                continue;
            end

            %% ------------------------------------------------------------
            %  Step 5: DC OPF
            % -------------------------------------------------------------
            if mode == "dc" || mode == "both"

                fprintf('\n------------------------------------------------------------\n');
                fprintf('[2] Running DC OPF\n');
                fprintf('------------------------------------------------------------\n');

                if isempty(opt.centralizedCostDC)

                    fprintf('> Running Centralized DCOPF...\n');

                    [resdc, fvaldc] = run_dcopf_centralized(centralized_case_from_file);

                    result_case.centralizedDC.res = resdc;
                    result_case.centralizedDC.fval = fvaldc;

                    if isfield(resdc, 'status') && resdc.status ~= 0
                        fprintf('!! Centralized DCOPF failed with status = %d.\n', resdc.status);
                        fprintf('!! Setting centralized DC cost = 0.\n');
                        centralizedCostDC = 0;
                    else
                        centralizedCostDC = fvaldc;
                    end

                else

                    centralizedCostDC = opt.centralizedCostDC;
                    fprintf('> Using user-provided centralized DC cost: %.10g\n', centralizedCostDC);

                end

                CentralizedDCCost(c) = centralizedCostDC;

                fprintf('> Running Distributed DCOPF...\n');
                fprintf('  Case      : %s\n', partitioned_case_name);
                fprintf('  rho       : %.10g\n', opt.rho_dc);
                fprintf('  threshold : %.10g\n', opt.thresholdDC);
                fprintf('  maxIter   : %d\n', opt.maxIterDC);
                fprintf('  holds     : %.10g\n', opt.holds);
                fprintf('  central f : %.10g\n', centralizedCostDC);

                figs_before = findall(0, 'Type', 'figure');

                tic_dc = tic;

                resDistDC = DistDCOPF( ...
                    partitioned_case_name, ...
                    opt.rho_dc, ...
                    opt.thresholdDC, ...
                    opt.maxIterDC, ...
                    centralizedCostDC, ...
                    opt.holds);

                runtime_dc = toc(tic_dc);

                figs_after = findall(0, 'Type', 'figure');
                new_figs_dc = setdiff(figs_after, figs_before);

                RuntimeDC(c) = runtime_dc;
                StatusDC(c) = "Success";

                result_case.distributedDC.res = resDistDC;
                result_case.distributedDC.runtime = runtime_dc;

                fprintf('> Distributed DCOPF runtime: %.4f seconds\n', runtime_dc);

                if opt.saveFigures
                    save_dplib_figures(new_figs_dc, partitioned_case_name, 'DC', opt.closeExtraFigures);
                end

            else

                StatusDC(c) = "Skipped";
                fprintf('\n[2] DC OPF skipped.\n');

            end

            %% ------------------------------------------------------------
            %  Step 6: AC OPF
            % -------------------------------------------------------------
            if mode == "ac" || mode == "both"

                fprintf('\n------------------------------------------------------------\n');
                fprintf('[3] Running AC OPF\n');
                fprintf('------------------------------------------------------------\n');

                if isempty(opt.centralizedCostAC)

                    fprintf('> Running Centralized ACOPF...\n');

                    [resac, fvalac] = run_acopf_centralized(centralized_case_from_file);

                    result_case.centralizedAC.res = resac;
                    result_case.centralizedAC.fval = fvalac;

                    if isfield(resac, 'status') && resac.status ~= 0
                        fprintf('!! Centralized ACOPF failed with status = %d.\n', resac.status);
                        fprintf('!! Setting centralized AC cost = 0.\n');
                        centralizedCostAC = 0;
                    else
                        centralizedCostAC = fvalac;
                    end

                else

                    centralizedCostAC = opt.centralizedCostAC;
                    fprintf('> Using user-provided centralized AC cost: %.10g\n', centralizedCostAC);

                end

                CentralizedACCost(c) = centralizedCostAC;

                fprintf('> Running Distributed ACOPF...\n');
                fprintf('  Case      : %s\n', partitioned_case_name);
                fprintf('  rho       : %.10g\n', opt.rho_ac);
                fprintf('  threshold : %.10g\n', opt.thresholdAC);
                fprintf('  maxIter   : %d\n', opt.maxIterAC);
                fprintf('  holds     : %.10g\n', opt.holds);
                fprintf('  central f : %.10g\n', centralizedCostAC);

                figs_before = findall(0, 'Type', 'figure');

                tic_ac = tic;

                resDistAC = DistACOPF( ...
                    partitioned_case_name, ...
                    opt.rho_ac, ...
                    opt.thresholdAC, ...
                    opt.maxIterAC, ...
                    centralizedCostAC, ...
                    opt.holds);

                runtime_ac = toc(tic_ac);

                figs_after = findall(0, 'Type', 'figure');
                new_figs_ac = setdiff(figs_after, figs_before);

                RuntimeAC(c) = runtime_ac;
                StatusAC(c) = "Success";

                result_case.distributedAC.res = resDistAC;
                result_case.distributedAC.runtime = runtime_ac;

                fprintf('> Distributed ACOPF runtime: %.4f seconds\n', runtime_ac);

                if opt.saveFigures
                    save_dplib_figures(new_figs_ac, partitioned_case_name, 'AC', opt.closeExtraFigures);
                end

            else

                StatusAC(c) = "Skipped";
                fprintf('\n[3] AC OPF skipped.\n');

            end

            fprintf('\n> Finished case: %s\n', case_name);

        catch ME

            Message(c) = string(ME.message);

            if StatusPartition(c) == ""
                StatusPartition(c) = "Failed";
            end

            if StatusDC(c) == ""
                StatusDC(c) = "Failed";
            end

            if StatusAC(c) == ""
                StatusAC(c) = "Failed";
            end

            fprintf('\nFAILED for case %s:\n%s\n', case_name, ME.message);

            result_case.error = ME;

        end

        all_results{c} = result_case;

    end

    %% ============================================================
    %  Build summary table
    % ============================================================
    DPLibMainSummary = table( ...
        CaseName, ...
        PartitionedCaseName, ...
        Regions, ...
        TieLines, ...
        RefBusRegion, ...
        TotalBuses, ...
        TotalBranchesRegional, ...
        TotalBranchesOriginal, ...
        TotalGenerators, ...
        PartitionTimeSec, ...
        CentralizedDCCost, ...
        CentralizedACCost, ...
        RuntimeDC, ...
        RuntimeAC, ...
        StatusPartition, ...
        StatusDC, ...
        StatusAC, ...
        Message);

    disp(' ');
    disp('============================================================');
    disp('DPLib main runner summary');
    disp('============================================================');
    disp(DPLibMainSummary);

    %% ============================================================
    %  Save summary files
    % ============================================================
    if opt.saveSummary

        summary_name = opt.summaryName;

        save([summary_name '.mat'], ...
            'DPLibMainSummary', ...
            'all_results');

        writetable(DPLibMainSummary, [summary_name '.csv']);
        writetable(DPLibMainSummary, [summary_name '.xlsx']);

        fprintf('\nSaved summary files:\n');
        fprintf('  %s.mat\n', summary_name);
        fprintf('  %s.csv\n', summary_name);
        fprintf('  %s.xlsx\n', summary_name);

    end

    %% ============================================================
    %  Print compact LaTeX table
    % ============================================================
    print_latex_table(DPLibMainSummary);

    fprintf('\n=============================\n');
    fprintf('All cases completed.\n');
    fprintf('=============================\n');

end


function save_dplib_figures(new_figs, caseName, opfType, closeExtraFigures)
%SAVE_DPLIB_FIGURES
% Saves only DPLib error, optimality gap, and rho figures.

    case_tex = strrep(caseName, '_', '\_');

    for iFig = 1:numel(new_figs)

        f = new_figs(iFig);
        figure(f);
        ax = gca;

        ylbl = get(ax, 'YLabel');
        ystr = ylbl.String;

        if iscell(ystr)
            ystr = strjoin(ystr, ' ');
        elseif isstring(ystr)
            ystr = char(ystr);
        end

        if isempty(ystr)
            if closeExtraFigures
                close(f);
            end
            continue;
        end

        isRhoFig = contains(ystr, '\rho', 'IgnoreCase', true);
        isGapFig = contains(ystr, 'Optimality Gap', 'IgnoreCase', true);

        if strcmpi(opfType, 'DC')
            isErrFig = contains(ystr, 'Error', 'IgnoreCase', true) || ...
                       contains(ystr, 'Primal residual', 'IgnoreCase', true);
        else
            isErrFig = contains(ystr, 'Error', 'IgnoreCase', true);
        end

        if ~(isRhoFig || isGapFig || isErrFig)
            if closeExtraFigures
                close(f);
            end
            continue;
        end

        ax.FontSize        = 14;
        ax.XLabel.FontSize = 16;
        ax.YLabel.FontSize = 16;
        ax.Title.FontSize  = 14;
        ax.LineWidth       = 1.5;

        if isRhoFig
            title(sprintf('%s (%s)', case_tex, upper(opfType)), ...
                'Interpreter', 'latex');
            save_name = sprintf('%s_%s_rho.png', caseName, upper(opfType));

        elseif isGapFig
            title(sprintf('%s (%s)', case_tex, upper(opfType)), ...
                'Interpreter', 'latex');
            save_name = sprintf('%s_%s_gap.png', caseName, upper(opfType));

        else
            title(sprintf('%s (%s)', case_tex, upper(opfType)), ...
                'Interpreter', 'latex');
            save_name = sprintf('%s_%s_error.png', caseName, upper(opfType));
        end

        fprintf('  > Saving %s figure: %s\n', upper(opfType), save_name);
        saveas(f, save_name);

    end

end

function print_latex_table(T)
%PRINT_LATEX_TABLE
% Prints compact LaTeX table from the summary table.

    fprintf('\n\n%% ============================================================\n');
    fprintf('%% LaTeX table for generated and solved DPLib cases\n');
    fprintf('%% ============================================================\n\n');

    fprintf('\\begin{table*}[!t]\n');
    fprintf('\\centering\n');
    fprintf('{\\color{blue}\n');
    fprintf('\\caption{Generated DPLib benchmark cases and distributed OPF solution summary.}\n');
    fprintf('\\label{tab:dplib_main_runner_summary}\n');
    fprintf('\\footnotesize\n');
    fprintf('\\setlength{\\tabcolsep}{4pt}\n');
    fprintf('\\renewcommand{\\arraystretch}{1.08}\n');
    fprintf('\\begin{tabular}{lrrrrrrrr}\n');
    fprintf('\\toprule\n');
    fprintf('\\textbf{Case} & \\textbf{Buses} & \\textbf{Branches} & \\textbf{Gens.} & \\textbf{Regions} & \\textbf{Tie-Lines} & \\textbf{Part. Time} & \\textbf{DC Time} & \\textbf{AC Time} \\\\\n');
    fprintf('\\midrule\n');

    for i = 1:height(T)

        clean_case = strrep(char(T.CaseName(i)), '_', '\\_');

        if T.StatusPartition(i) == "Success" || T.StatusPartition(i) == "Skipped"

            fprintf('\\texttt{%s} & %.0f & %.0f & %.0f & %.0f & %.0f & %.2f & %.2f & %.2f \\\\\n', ...
                clean_case, ...
                T.TotalBuses(i), ...
                T.TotalBranchesOriginal(i), ...
                T.TotalGenerators(i), ...
                T.Regions(i), ...
                T.TieLines(i), ...
                T.PartitionTimeSec(i), ...
                T.RuntimeDC(i), ...
                T.RuntimeAC(i));

        else

            fprintf('\\texttt{%s} & -- & -- & -- & %.0f & -- & -- & -- & -- \\\\\n', ...
                clean_case, ...
                T.Regions(i));

        end
    end

    fprintf('\\bottomrule\n');
    fprintf('\\end{tabular}\n');
    fprintf('}\n');
    fprintf('\\end{table*}\n');

end

