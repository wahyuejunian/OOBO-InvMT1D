function PlotMTinv_syn(freq, app_data, phase_data, X_rho_best, X_thick_best, misfit_oobo, rsyn, tsyn)
%PLOTMTINV_SYN Visualizes inversion results for synthetic MT data
%
% INPUT:
%   freq         : Vector of frequency values
%   app_data     : Observed apparent resistivity values
%   phase_data   : Observed phase values (in degrees)
%   X_rho_best   : Matrix of inverted resistivity models (each row = one model)
%   X_thick_best : Matrix of inverted thickness models (each row = one model)
%   misfit_oobo  : Misfit values for each model
%   rsyn         : Resistivity values for the true (synthetic) model
%   tsyn         : Thickness values for the true (synthetic) model

% Find the best model (lowest misfit)
[~, idx] = min(misfit_oobo(:, end));
best_rho = X_rho_best(idx, :);
best_thc = X_thick_best(idx, :);

% Forward modeling of best model to get apparent resistivity and phase
[appres_best, phase_best] = MT1D(best_rho, best_thc, freq);

figure(1)

% Number of inversion runs
[row, ~] = size(X_rho_best);
run = row;

% --- Plot 1: Apparent Resistivity curves ---
for ipop = 1:run
    % Forward modeling for each inversion model
    L_rho = X_rho_best(ipop, :);
    L_thick = X_thick_best(ipop, :);
    [apparentResistivity, phase] = MT1D(L_rho, L_thick, freq);

    % Store results for all models
    rhoapp_cal(ipop, :) = apparentResistivity;
    phase_mod(ipop, :) = phase;

    % Plot calculated resistivity curves in gray
    subplot(3, 3, [1 5])
    loglog(1 ./ freq, rhoapp_cal, '-', 'Color', '#C0C0C0', 'LineWidth', 2);
end

% Plot observed and best-fit calculated apparent resistivity
subplot(3, 3, [1 5])
hold on;
a = loglog(1 ./ freq, app_data, 'ok', 'MarkerSize', 3.5, 'LineWidth', 3, 'DisplayName', 'AppRes Obs'); hold on
b = loglog(1 ./ freq, appres_best, '-r', 'MarkerSize', 3.5, 'LineWidth', 2.5, 'DisplayName', 'AppRes Cal');
axis([10^-4.5 10^3.5 10^0 10^4]);
xlabel('Periods (sec)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('App. Resistivity (\Omega·m)', 'FontSize', 12, 'FontWeight', 'Bold');
legend([a, b])
set(legend, 'Location', 'Northwest', 'fontsize', 8);
set(gca, 'LineWidth', 1.5);
hold on

% --- Plot 2: Phase curves ---
subplot(3, 3, [7 8])
c = loglog(1 ./ freq, phase_data, 'ok', 'MarkerSize', 3.5, 'LineWidth', 3, 'DisplayName', 'Phase Obs'); hold on
d = loglog(1 ./ freq, phase_best, '-r', 'MarkerSize', 3.5, 'LineWidth', 2.5, 'DisplayName', 'Phase Cal');
legend([c, d])
axis([10^-4.5 10^3.5 0 100]);
set(gca, 'YScale', 'linear');
xlabel('Periods (sec)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Phase (deg)', 'FontSize', 12, 'FontWeight', 'Bold');
set(legend, 'Location', 'Northwest', 'fontsize', 8);
set(gca, 'LineWidth', 1.5);

% --- Plot 3: Resistivity vs. Depth (Model Comparison) ---
subplot(3, 3, [3 9])
for irun = 1:run
    rr_inv = [0, X_rho_best(irun, :)];  % Resistivity
    tt_inv = [0, cumsum(X_thick_best(irun, :)), max(X_thick_best(irun, :)) * 500];  % Depth
    a = stairs(rr_inv, tt_inv, '-', 'Color', '#C0C0C0', 'LineWidth', 2, 'DisplayName', '30 inversion models');
    hold on;
end

% Plot true (synthetic) model
rr = [0, rsyn];
tt = [0, cumsum(tsyn), max(tsyn) * 10];
b = stairs(rr, tt, '-.k', 'LineWidth', 3, 'DisplayName', 'Synthetic Model');
hold on;

% Plot best model from inversion
rr_best = [0, best_rho];
tt_best = [0, cumsum(best_thc), max(best_thc) * 500];
c = stairs(rr_best, tt_best, '-r', 'LineWidth', 3, 'DisplayName', 'Best Model');

% Formatting
title('\bf Model');
axis([10^-0.5 10^5.5 0 3000])
xlabel('Resistivity (\Omega·m)', 'FontWeight', 'bold', 'FontSize', 10);
ylabel('Depth (m)', 'FontWeight', 'bold', 'FontSize', 10);
set(gca, 'XScale', 'log');
set(gca, 'YDir', 'reverse');
set(gca, 'LineWidth', 1.5);
legend([a, b, c]);
set(legend, 'Location', 'Southwest', 'fontsize', 8);
end
