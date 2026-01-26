clear; close all; clc;
format long g

%% ============================================================
%  Global defaults: smaller fonts everywhere
% ============================================================
set(groot, 'DefaultAxesFontSize', 9);          % tick labels
set(groot, 'DefaultTextFontSize', 9);
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.0);
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.0);

%% ============================================================
%  List of multi-region PGLib test cases
% ============================================================
case_list = {
    'pglib_opf_case9241_pegase_18regions', ...
};

%{
  'pglib_opf_case200_tamu_3regions', ...
  'pglib_opf_case300_ieee_3regions', ...
  'pglib_opf_case500_tamu_5regions', ...
  'pglib_opf_case1354_pegase_8regions', ...
  'pglib_opf_case2869_pegase_12regions', ...
  'pglib_opf_case4661_sdet_14regions', ...
%}

%% ============================================================
%%  Rho values for each case (same order as case_list)
%% ============================================================


thresholdAC = 1e-4;
thresholdDC = 1e-4;
maxIterDC = 10000;
maxIterAC = 10000;
holds = 1*10^(6);

results = struct();

for k = 1:length(case_list)

    Casename = case_list{k};
    rho_dc   = 1;
    rho_ac   = 1*10^(2);

    fprintf('\n=============================================\n');
    fprintf('Running Case (%d/%d): %s\n', k, length(case_list), Casename);
    fprintf('  -> Using rho = %.4g for DC OPF and rho = %.4g for AC OPF\n', rho_dc, rho_ac);
    fprintf('=============================================\n');

    %% ----------------------------------------------
    %% Load MATPOWER data
    %% ----------------------------------------------
    loadedData = load(Casename);        % .mat file with "filename" inside
    CasenameC  = loadedData.filename;   % MATPOWER .m case name

    % % ----------------------------------------------
    % %% CENTRALIZED DC-OPF
    % %% ----------------------------------------------
    % fprintf('> Running Centralized DCOPF...\n');
    % [resdc, fvaldc] = run_dcopf_centralized(CasenameC);
    % 
    % if resdc.status ~= 0
    %     fprintf('!! Centralized DCOPF FAILED (status = %d). Setting fvaldc = 0.\n', ...
    %         resdc.status);
    %     fvaldc = 0;    % no gap computation in DistDCOPF if this is 0
    % end
    % 
    % %% ----------------------------------------------
    % %% DISTRIBUTED DC-OPF
    % %% ----------------------------------------------
    % fprintf('> Running Distributed DCOPF (rho = %.4g)...\n', rho_dc);
    % 
    % figs_before = findall(0, 'Type', 'figure');
    % resDistDC   = DistDCOPF(Casename, rho_dc, thresholdDC, maxIterDC, fvaldc,holds);
    % figs_after  = findall(0, 'Type', 'figure');
    % 
    % new_figs_dc = setdiff(figs_after, figs_before);  % figures created by DistDCOPF
    % 
    % % TeX-safe case label
    % case_tex = strrep(Casename, '_', '\_');
    % 
    % % Only keep three figures: error, gap, rho
    % for iFig = 1:numel(new_figs_dc)
    %     f = new_figs_dc(iFig);
    %     figure(f);
    %     ax = gca;
    % 
    %     % --- Extract ylabel string robustly ---
    %     ylbl = get(ax, 'YLabel');
    %     ystr = ylbl.String;
    %     if iscell(ystr)
    %         ystr = strjoin(ystr, ' ');
    %     elseif isstring(ystr)
    %         ystr = char(ystr);
    %     end
    % 
    %     if isempty(ystr)
    %         % Nothing to identify: skip
    %         close(f);
    %         continue;
    %     end
    % 
    %     % --- Classify figure type ---
    %     isRhoFig  = contains(ystr, '\rho',           'IgnoreCase', true);
    %     isGapFig  = contains(ystr, 'Optimality Gap', 'IgnoreCase', true);
    %     % DC error plot label is "Primal residual (max error)" in DistDCOPF
    %     isErrFig  = contains(ystr, 'Error',          'IgnoreCase', true) || ...
    %                 contains(ystr, 'Primal residual', 'IgnoreCase', true);
    % 
    %     if ~(isRhoFig || isGapFig || isErrFig)
    %         % This is some auxiliary figure: don't save it
    %         close(f);
    %         continue;
    %     end
    % 
    %     % --- Nice, compact font sizes ---
    %     ax.FontSize        = 14;   % tick labels
    %     ax.XLabel.FontSize = 16;
    %     ax.YLabel.FontSize = 16;
    %     ax.Title.FontSize  = 14;
    %     ax.LineWidth       = 1.5;
    % 
    %     % --- Set LaTeX title & filename based on type ---
    %     if isRhoFig
    %         title(sprintf('%s (DC) ', case_tex), ...
    %               'Interpreter', 'latex');
    %         save_name = sprintf('%s_DC_rho.png', Casename);
    %     elseif isGapFig
    %         title(sprintf('%s (DC)', case_tex), ...
    %               'Interpreter', 'latex');
    %         save_name = sprintf('%s_DC_gap.png', Casename);
    %     else    % error figure
    %         title(sprintf('%s (DC)', case_tex), ...
    %               'Interpreter', 'latex');
    %         save_name = sprintf('%s_DC_error.png', Casename);
    %     end
    % 
    %     fprintf('  > Saving DC figure: %s\n', save_name);
    %     saveas(f, save_name);
    % end

    %% ----------------------------------------------
    %% CENTRALIZED AC-OPF
    %% ----------------------------------------------
    fprintf('> Running Centralized ACOPF...\n');
    %[resac, fvalac] = run_acopf_centralized(CasenameC);

    %% ----------------------------------------------
    %% DISTRIBUTED AC-OPF
    %% ----------------------------------------------
    fprintf('> Running Distributed ACOPF (rho = %.4g)...\n', rho_ac);

    figs_before = findall(0, 'Type', 'figure');
    resDistAC   = DistACOPF(Casename, rho_ac, thresholdAC, maxIterAC, 6243090,holds);
    figs_after  = findall(0, 'Type', 'figure');

    new_figs_ac = setdiff(figs_after, figs_before);  % figures created by DistACOPF

    case_tex = strrep(Casename, '_', '\_');

    for iFig = 1:numel(new_figs_ac)
        f = new_figs_ac(iFig);
        figure(f);
        ax = gca;

        % --- Extract ylabel string robustly ---
        ylbl = get(ax, 'YLabel');
        ystr = ylbl.String;
        if iscell(ystr)
            ystr = strjoin(ystr, ' ');
        elseif isstring(ystr)
            ystr = char(ystr);
        end

        if isempty(ystr)
            close(f);
            continue;
        end

        % In DistACOPF:
        %   - error plot y-label is "Error"
        %   - gap plot y-label is "Optimality Gap"
        %   - rho plot y-label has '\rho'
        isRhoFig = contains(ystr, '\rho',           'IgnoreCase', true);
        isGapFig = contains(ystr, 'Optimality Gap', 'IgnoreCase', true);
        isErrFig = contains(ystr, 'Error',          'IgnoreCase', true);

        if ~(isRhoFig || isGapFig || isErrFig)
            % Skip any other figures
            close(f);
            continue;
        end

        % --- Font sizes (same as DC) ---
        ax.FontSize        = 14;
        ax.XLabel.FontSize = 16;
        ax.YLabel.FontSize = 16;
        ax.Title.FontSize  = 14;
        ax.LineWidth       = 1.5;

        % --- LaTeX titles + filenames ---
        if isRhoFig
            title(sprintf('%s (AC)', case_tex), ...
                  'Interpreter', 'latex');
            save_name = sprintf('%s_AC_rho.png', Casename);
        elseif isGapFig
            title(sprintf('%s (AC) ', case_tex), ...
                  'Interpreter', 'latex');
            save_name = sprintf('%s_AC_gap.png', Casename);
        else    % error
            title(sprintf('%s (AC)', case_tex), ...
                  'Interpreter', 'latex');
            save_name = sprintf('%s_AC_error.png', Casename);
        end

        fprintf('  > Saving AC figure: %s\n', save_name);
        saveas(f, save_name);
    end

    fprintf('> Finished Case: %s\n', Casename);
end

fprintf('\n=============================\n');
fprintf('All cases completed.\n');
fprintf('=============================\n');
