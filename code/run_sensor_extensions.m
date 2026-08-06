%% run_sensor_extensions.m
% Sensor-oriented extensions for the target-certification paper.
% Produces sensor_extensions_results.csv and sensor_extensions.pdf/png.
clear; clc; close all;
rng(20260806,'twister');

Nmc = 5000;
M = 10;
K = 80;
deltaBar = 0.05;
rhoValues = 1.0:0.1:2.4;

heterogeneityNames = { ...
    'Homogeneous', ...
    'Mild heterogeneity', ...
    'Strong heterogeneity', ...
    'One degraded sensor'};

heterogeneityMultipliers = { ...
    ones(1,M), ...
    [0.8*ones(1,M/2), 1.2*ones(1,M/2)], ...
    [0.5*ones(1,M/2), 1.5*ones(1,M/2)], ...
    [ones(1,M-1), 1.5]};

distributionNames = { ...
    'Center', ...
    'Uniform disk', ...
    'Boundary', ...
    'Fixed magnitude'};

distributionCodes = {'center','uniform','boundary','fixed'};

nRows = numel(rhoValues) * ...
    (numel(heterogeneityNames) + numel(distributionNames));
experiment = strings(nRows,1);
caseName = strings(nRows,1);
rhoColumn = zeros(nRows,1);
Psys = zeros(nRows,1);
Tsys = nan(nRows,1);
row = 0;

% Experiment 1: heterogeneous sensor bounds.
for h = 1:numel(heterogeneityNames)
    multipliers = heterogeneityMultipliers{h};
    for r = 1:numel(rhoValues)
        row = row + 1;
        [Psys(row),Tsys(row)] = simulate_configuration( ...
            Nmc,M,K,deltaBar,rhoValues(r),multipliers,'uniform');
        experiment(row) = "heterogeneity";
        caseName(row) = heterogeneityNames{h};
        rhoColumn(row) = rhoValues(r);
    end
end

% Experiment 2: different error distributions under the same bound.
for d = 1:numel(distributionNames)
    for r = 1:numel(rhoValues)
        row = row + 1;
        [Psys(row),Tsys(row)] = simulate_configuration( ...
            Nmc,M,K,deltaBar,rhoValues(r),ones(1,M), ...
            distributionCodes{d});
        experiment(row) = "distribution";
        caseName(row) = distributionNames{d};
        rhoColumn(row) = rhoValues(r);
    end
end

T = table(experiment,caseName,rhoColumn,Psys,Tsys, ...
    'VariableNames',{'experiment','case_name','rho', ...
    'P_sys','T_sys_conditional'});
writetable(T,'sensor_extensions_results.csv');

plot_sensor_extensions(T,heterogeneityNames,distributionNames);
fprintf(['Generated sensor_extensions_results.csv, ', ...
    'sensor_extensions.pdf, and sensor_extensions.png.\n']);

function [Psys,Tsys] = simulate_configuration( ...
    Nmc,M,K,deltaBar,rho,multipliers,distributionCode)

    delta = deltaBar * multipliers(:)';
    deltaStar = rho * deltaBar;
    threshold = deltaStar - delta;

    % If one assigned sensor has a negative threshold, the complete task
    % cannot be certified with the proposed pairwise rule.
    if any(threshold < 0)
        Psys = 0;
        Tsys = NaN;
        return;
    end

    starts = -1 + 2*rand(Nmc,M,2);
    targets = -1 + 2*rand(Nmc,M,2);
    tau = inf(Nmc,M);
    active = true(Nmc,M);

    for k = 0:K
        if ~any(active,'all')
            break;
        end

        alpha = k/K;
        truePosition = starts + alpha*(targets-starts);
        U = rand(Nmc,M);
        theta = 2*pi*rand(Nmc,M);

        switch distributionCode
            case 'center'
                normalizedRadius = U.^2;
            case 'uniform'
                normalizedRadius = sqrt(U);
            case 'boundary'
                normalizedRadius = U.^(1/4);
            case 'fixed'
                normalizedRadius = ones(Nmc,M);
            otherwise
                error('Unknown distribution code.');
        end

        radius = normalizedRadius .* delta;
        errorVector = cat(3,radius.*cos(theta),radius.*sin(theta));
        measurement = truePosition + errorVector;
        distanceToTarget = vecnorm(measurement-targets,2,3);
        triggered = active & ...
            (distanceToTarget <= threshold + 1e-14);
        tau(triggered) = k;
        active(triggered) = false;
    end

    successfulTrial = all(isfinite(tau),2);
    Psys = mean(successfulTrial);

    if any(successfulTrial)
        Tsys = mean(max(tau(successfulTrial,:),[],2)/K);
    else
        Tsys = NaN;
    end
end

function plot_sensor_extensions(T,heterogeneityNames,distributionNames)
%PLOT_SENSOR_EXTENSIONS Create the two-panel publication figure using
%the same visual style adopted in numerical_tradeoffs.

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 29 8.5], ...
    'MenuBar','none', ...
    'ToolBar','none');

layout = tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Same family of line and marker styles used in the previous figure.
styles = {'-o','--s','-.^',':d'};

% Same color palette used in the previous numerical figure.
% Blue, orange, yellow, and purple are combined with distinct
% line styles and markers for accessibility and grayscale printing.
colors = [ ...
    0.0000 0.4470 0.7410; ... % blue
    0.8500 0.3250 0.0980; ... % orange
    0.9290 0.6940 0.1250; ... % yellow
    0.4940 0.1840 0.5560];    % purple

boundaryColor = [0.30 0.30 0.30];

% ---------------------------------------------------------------
% Panel (a): heterogeneous position-sensing bounds
% ---------------------------------------------------------------
ax1 = nexttile(layout);
hold(ax1,'on');

for c = 1:numel(heterogeneityNames)
    rows = T.experiment=="heterogeneity" & ...
        T.case_name==heterogeneityNames{c};

    plot(ax1,T.rho(rows),T.P_sys(rows),styles{c}, ...
        'Color',colors(c,:), ...
        'LineWidth',1.2, ...
        'MarkerSize',4, ...
        'DisplayName',heterogeneityNames{c});
end

% Boundary shown without an in-plot annotation and included in the legend.
xline(ax1,2,'--', ...
    'Color',boundaryColor, ...
    'LineWidth',1.2, ...
    'DisplayName','$\rho=2$');

xlabel(ax1,'$\rho=\delta^*/\bar{\delta}$', ...
    'Interpreter','latex');
ylabel(ax1,'$P_{\mathrm{sys}}$', ...
    'Interpreter','latex');
title(ax1,'(a) Sensor heterogeneity');

xlim(ax1,[1 2.4]);
ylim(ax1,[0 1.03]);
grid(ax1,'on');
box(ax1,'on');
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 9;
ax1.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax1);

legend(ax1, ...
    'Location','northwest', ...
    'Interpreter','latex', ...
    'FontSize',8, ...
    'Box','on');

% ---------------------------------------------------------------
% Panel (b): bounded-error radial distributions
% ---------------------------------------------------------------
ax2 = nexttile(layout);
hold(ax2,'on');

for c = 1:numel(distributionNames)
    rows = T.experiment=="distribution" & ...
        T.case_name==distributionNames{c};

    plot(ax2,T.rho(rows),T.P_sys(rows),styles{c}, ...
        'Color',colors(c,:), ...
        'LineWidth',1.2, ...
        'MarkerSize',4, ...
        'DisplayName',distributionNames{c});
end

xline(ax2,2,'--', ...
    'Color',boundaryColor, ...
    'LineWidth',1.2, ...
    'DisplayName','$\rho=2$');

xlabel(ax2,'$\rho=\delta^*/\delta$', ...
    'Interpreter','latex');
ylabel(ax2,'$P_{\mathrm{sys}}$', ...
    'Interpreter','latex');
title(ax2,'(b) Bounded-error distributions');

xlim(ax2,[1 2.4]);
ylim(ax2,[0 1.03]);
grid(ax2,'on');
box(ax2,'on');
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = 9;
ax2.Toolbar.Visible = 'off';
disableDefaultInteractivity(ax2);

legend(ax2, ...
    'Location','northwest', ...
    'Interpreter','latex', ...
    'FontSize',8, ...
    'Box','on');

% Publication-ready exports.
drawnow;
exportgraphics(fig,'sensor_extensions.pdf', ...
    'ContentType','vector');
exportgraphics(fig,'sensor_extensions.png', ...
    'Resolution',600);
end
