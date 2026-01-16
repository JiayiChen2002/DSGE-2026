function analyze_sensitivity_results()
    % Analyze & visualize shock_sensitivity_irf_stats.csv
    clc;

    inCsv = 'shock_sensitivity_irf_stats.csv';
    assert(exist(inCsv,'file')==2, 'Cannot find %s', inCsv);

    T = readtable(inCsv);

    % ---- cleaning / conventions ----
    % Treat tiny numerical noise as zero (e.g. pz under NO_POLICY)
    tiny = 1e-12;
    T.impact(abs(T.impact)<tiny) = 0;
    T.peak_abs(abs(T.peak_abs)<tiny) = 0;
    T.cumulative(abs(T.cumulative)<tiny) = 0;

    % Mark variables that are not meaningful under NO_POLICY
    % (u and pz are mechanically fixed by regime definition)
    isNo = strcmp(T.regime,'NO_POLICY');
    isU  = strcmp(T.variable,'u');
    isPz = strcmp(T.variable,'pz');
    T(isNo & (isU | isPz), :) = [];

    % ---- choose baseline point ----
    rho0 = 0.95;
    sig0 = 0.01;

    % =========================
    % Figure 1: regime comparison at baseline
    % =========================
    vars_show = {'y','z','pz'}; % pz exists only for policy regimes
    M = T(T.rho_a==rho0 & T.sigma_a==sig0 & ismember(T.variable, vars_show), :);

    figure;
    % build grouped bars: rows=regime, cols=variable
    regimes = unique(M.regime,'stable');
    V = vars_show;
    A = NaN(numel(regimes), numel(V));
    for i=1:numel(regimes)
        for j=1:numel(V)
            ix = strcmp(M.regime, regimes{i}) & strcmp(M.variable, V{j});
            if any(ix)
                A(i,j) = M.peak_abs(find(ix,1));
            end
        end
    end

    bar(A);
    set(gca,'XTickLabel', regimes);
    legend(V,'Location','best');
    xlabel('Regime');
    ylabel('Peak |IRF|');
    title(sprintf('Peak IRF comparison (rho_a=%.2f, sigma_a=%.3f)', rho0, sig0));
    grid on;
% =========================
% Figure 2: heatmaps over (rho, sigma) for peak_abs( VAR )
% (with unified color scale across regimes for each variable)
% =========================

regimes_all = unique(T.regime,'stable');

% 你想分析的核心变量：按环境机制优先
vars_heat = {'y','z','pz','u'};   % 你也可以改成 {'mc','pie','rnom'} 等

for vv = 1:numel(vars_heat)
    vname = vars_heat{vv};

    % ---- NEW: compute global max for this variable across all regimes ----
    Tall = T(strcmp(T.variable, vname), :);
    vmax = max(Tall.peak_abs(~isnan(Tall.peak_abs)));
    if isempty(vmax) || vmax <= 0
        vmax = 1e-12;  % 防止全是0/空导致caxis报错
    end

    figure('Name', ['heatmap_peakabs_' vname]);

    for r = 1:numel(regimes_all)
        reg = regimes_all{r};

        % 取出该制度 + 该变量
        S = T(strcmp(T.regime,reg) & strcmp(T.variable,vname), :);

        subplot(2,2,r);

        % 如果该 regime 下该变量不存在/被你清理掉（如 NO_POLICY 的 pz/u）
        if isempty(S)
            axis off;
            title(sprintf('%s: peak\\_abs(%s) [N/A]', reg, vname), 'Interpreter','tex');
            continue;
        end

        rho_vals = unique(S.rho_a);
        sig_vals = unique(S.sigma_a);

        Z = NaN(numel(rho_vals), numel(sig_vals));
        for i=1:numel(rho_vals)
            for j=1:numel(sig_vals)
                ix = (S.rho_a==rho_vals(i) & S.sigma_a==sig_vals(j));
                if any(ix)
                    Z(i,j) = S.peak_abs(find(ix,1));
                end
            end
        end

        imagesc(sig_vals, rho_vals, Z);
        set(gca,'YDir','normal');
        colorbar;

        % ---- NEW: unified color axis (0 to global max) ----
        caxis([0, vmax]);

        xlabel('\sigma_a');
        ylabel('\rho_a');
        title(sprintf('%s: peak\\_abs(%s)', reg, vname), 'Interpreter','tex');
    end
end



    % =========================
    % Figure 3: scaling check (peak_abs(y) vs sigma_a)
    % =========================
    figure;
    hold on;
    for r=1:numel(regimes_all)
        reg = regimes_all{r};
        S = T(strcmp(T.regime,reg) & strcmp(T.variable,'y') & T.rho_a==rho0, :);
        [sig_sorted, idx] = sort(S.sigma_a);
        plot(sig_sorted, S.peak_abs(idx), '-o', 'DisplayName', reg);
    end
    hold off;
    xlabel('\sigma_a');
    ylabel('peak_abs(y)');
    title(sprintf('Linearity check at rho_a=%.2f', rho0));
    legend('Location','best');
    grid on;

    % ---- export figures (optional) ----
    % saveas(figure(1),'fig_baseline_peakabs.png');
    % saveas(figure(2),'fig_heatmaps_peakabs_y.png');
    % saveas(figure(3),'fig_linearity_peakabs_y.png');

    disp('Done: generated 3 figures.');
end


