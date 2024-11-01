% Estimated data from the image
x = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]; % x-axis points
Cp_upper = [-1.5, -1.2, -1.1, -0.9, -0.7, -0.6, -0.6, -0.7, -0.8, -0.9, -1.0]; % Upper surface Cp
Cp_lower = [0.5, 0.4, 0.3, 0.2, 0.0, -0.1, -0.2, -0.3, -0.35, -0.4, -0.45]; % Lower surface Cp

% Plot upper surface
figure;
plot(x, Cp_upper, 'k-', 'LineWidth', 1.5); % Present computation (upper)
hold on;

% Plot lower surface
plot(x, Cp_lower, 'k--', 'LineWidth', 1.5); % Present computation (lower)

% Experimental data points (estimated)
x_exp = [0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9];
Cp_exp_upper = [-1.5, -1.3, -1.0, -0.8, -0.7, -0.85, -0.9];
Cp_exp_lower = [0.5, 0.4, 0.2, 0.1, -0.1, -0.25, -0.3];

% Plot experimental data for upper and lower surfaces
plot(x_exp, Cp_exp_upper, 'ko', 'MarkerSize', 6, 'DisplayName', 'Upper surface');
plot(x_exp, Cp_exp_lower, 'ks', 'MarkerSize', 6, 'DisplayName', 'Lower surface');

% Add labels and legend
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$C_p$', 'Interpreter', 'latex', 'FontSize', 14);
legend({'Present comp.', 'Full boundary condition', 'Upper surface Exp.', 'Lower surface Exp.'}, ...
    'Location', 'northeast');

% Add phase and angle of attack text
text(0.05, 0.8, ['$\alpha(t) = 2.34^\circ$' newline 'Phase = 67.8$^\circ$'], ...
    'Interpreter', 'latex', 'FontSize', 12);

% Adjust plot limits and appearance
xlim([0 1]);
ylim([-2 1]);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

grid on;
hold off;
