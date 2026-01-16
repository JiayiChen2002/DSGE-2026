function run_shock_sensitivity()
% =============================================================
% Shock sensitivity runner for Annicchiarico & Di Dio (2015) replication
% Varies:
%   RHO_A  in rho_grid
%   stderr(eps_a) in sigma_grid
% Across regimes:
%   NO_POLICY, CAP_AND_TRADE, INTENSITY_TARGET, TAX_POLICY
%
% Output:
%   shock_sensitivity_irf_stats.csv
% =============================================================

    clc;

global oo_ M_ options_


    % ---- user settings ----
    templateFile = 'sens_template.mod';

    % Regimes you want to compare
    regimes = {'NO_POLICY','CAP_AND_TRADE','INTENSITY_TARGET','TAX_POLICY'};

    % Shock environment grid
    rho_grid   = [0.85, 0.95];
    sigma_grid = [0.005, 0.01, 0.02];

    % IRF extraction settings
    shockName = 'eps_a';
    vars = {'y','c','iv','l','mc','pie','z','u','pz','rnom'};
    T = 21; % must match stoch_simul(irf=21)

    % Calibration files (must exist in this folder)
    calib_no_policy = 'ann_dio_2015_calib_no_policy.inc';
    calib_policy    = 'ann_dio_2015_calib_environmental_policy.inc';

    % Output
    outCsv = 'shock_sensitivity_irf_stats.csv';

    % ---- sanity checks ----
    assert(exist(templateFile,'file')==2, 'Cannot find %s', templateFile);
    assert(exist(calib_no_policy,'file')==2, 'Cannot find %s', calib_no_policy);
    assert(exist(calib_policy,'file')==2, 'Cannot find %s', calib_policy);

    templateText = fileread(templateFile);

    % Preallocate a cell array for rows (fast enough for small grids)
    rows = {};
    rowIdx = 0;

    for r = 1:numel(regimes)
        regime = regimes{r};

        % Choose which calibration to load
        if strcmp(regime, 'NO_POLICY')
            calibFile = calib_no_policy;
        else
            calibFile = calib_policy;
        end

        % Build the regime macro block for Dynare preprocessor
        regimeBlock = make_regime_block(regime);

        for i = 1:numel(rho_grid)
            rhoA = rho_grid(i);

            for j = 1:numel(sigma_grid)
                sigmaA = sigma_grid(j);

                % Create a unique temp mod name (avoid dots in filenames)
                tmpMod = sprintf('tmp_sens_%s_rho%s_sig%s.mod', ...
                                 regime, num2str_safe(rhoA), num2str_safe(sigmaA));

                % Fill template
                modText = templateText;
                modText = strrep(modText, '%%REGIME_BLOCK%%', regimeBlock);
                modText = strrep(modText, '%%CALIB_FILE%%', calibFile);
                modText = strrep(modText, '%%RHO_A_VALUE%%', sprintf('%.6f', rhoA));
                modText = strrep(modText, '%%SIGMA_A_VALUE%%', sprintf('%.6f', sigmaA));

                % Write temp mod
                write_text(tmpMod, modText);

                % Clean previous Dynare globals to reduce cross-run contamination
                clear_global_dynare();

                % Run Dynare
                fprintf('\n=== Running: %s | RHO_A=%.3f | SIGMA_A=%.4f ===\n', regime, rhoA, sigmaA);
                dynare(tmpMod, 'noclearall');

                % Extract IRFs & summarize
                for v = 1:numel(vars)
                    varname = vars{v};
                    field = [varname '_' shockName];

                    if ~isfield(oo_.irfs, field)
                        error('IRF field oo_.irfs.%s not found. Check variable list or shock name.', field);
                    end

                    irf = oo_.irfs.(field);
                    irf = irf(1:min(T, numel(irf)));

                    s = summarize_irf(irf);

                    rowIdx = rowIdx + 1;
                    rows(rowIdx, :) = { ...
                        regime, rhoA, sigmaA, varname, ...
                        s.impact, s.peak_abs, s.t_peak, s.cum, s.halflife ...
                    };
                end

                % Optionally delete temp mod (keep if you want reproducibility)
                % delete(tmpMod);

            end
        end
    end

    % Build table and write CSV
    TBL = cell2table(rows, 'VariableNames', { ...
        'regime','rho_a','sigma_a','variable', ...
        'impact','peak_abs','t_peak','cumulative','halflife' ...
    });

    writetable(TBL, outCsv);
    fprintf('\nSaved results to %s\n', outCsv);
end

% ---------- helpers ----------

function block = make_regime_block(regime)
    % Creates Dynare macro definitions for the chosen regime
    switch regime
        case 'NO_POLICY'
            block = "@#define NO_POLICY = 1";
        case 'CAP_AND_TRADE'
            block = "@#define CAP_AND_TRADE = 1";
        case 'INTENSITY_TARGET'
            block = "@#define INTENSITY_TARGET = 1";
        case 'TAX_POLICY'
            block = "@#define TAX_POLICY = 1";
        otherwise
            error('Unknown regime: %s', regime);
    end
    block = sprintf('%s\n', block);
end

function s = summarize_irf(irf)
    % Returns summary statistics for an IRF vector
    % impact: first element
    % peak_abs: max absolute deviation
    % t_peak: index of peak (starting at 0 for first element)
    % cumulative: sum of IRF
    % halflife: periods after peak until abs(irf) falls below half peak

    s = struct();
    s.impact = irf(1);

    [pk, idx] = max(abs(irf));
    s.peak_abs = pk;
    s.t_peak = idx - 1; % treat first element as period 0

    s.cum = sum(irf);

    target = 0.5 * pk;
    after = abs(irf(idx:end));
    j = find(after <= target, 1, 'first');
    if isempty(j)
        s.halflife = NaN;
    else
        s.halflife = j - 1;
    end
end

function write_text(filename, txt)
    % Robust writer: works for char, string, cellstr across MATLAB versions
    if isstring(txt)
        txt = join(txt, sprintf('\n'));  % avoid newline() surprises
        txt = char(txt);
    elseif iscell(txt)
        txt = char(join(string(txt), sprintf('\n')));
    elseif ~ischar(txt)
        txt = char(string(txt));
    end

    % ---- IMPORTANT: strip non-ASCII chars that Dynare lexer cannot handle ----
    % keep: tab(9), LF(10), CR(13), and printable ASCII 32-126
    txt = regexprep(txt, '[^\x09\x0A\x0D\x20-\x7E]', '');

    fid = fopen(filename, 'w');
    if fid < 0
        error('Cannot open %s for writing.', filename);
    end

    fprintf(fid, '%s', txt);
    fclose(fid);
end

function s = num2str_safe(x)
    % convert number to a filename-safe string: 0.95 -> 0p95
    s = strrep(sprintf('%.3f', x), '.', 'p');
end

function clear_global_dynare()
    global M_ oo_ options_ ys0_ ex0_ it_
    clear global M_ oo_ options_ ys0_ ex0_ it_
end