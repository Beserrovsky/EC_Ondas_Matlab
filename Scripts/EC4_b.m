clear; clc; close all;

%% Medium parameters
u0 = 4*pi*1e-7;
e0 = 8.854e-12;

u1 = u0;       e1 = 3.702*e0;
u2 = u0;       e2 = e0;

n1 = sqrt(u1/e1);
n2 = sqrt(u2/e2);

%% θ1 from 0° to 90°
theta1_degs = 0:0.1:90;
theta1_rads = deg2rad(theta1_degs);

%% Transmitted angle θ2
theta2_rads = arrayfun(@(a) getTransmittedAngle(a, e1, e2), theta1_degs);

%% Load impedance Zl = n2 cos(theta2)
Zl = n2 .* cos(theta2_rads);
Zl_real = real(Zl);
Zl_imag = imag(Zl);

%% Angles to highlight
angles_mark = [0, 31.4, 90];
theta2_mark = asin((n1/n2).*sin(deg2rad(angles_mark)));
Zl_mark = n2 .* cos(theta2_mark);

Zl_mark_real = real(Zl_mark);
Zl_mark_imag = imag(Zl_mark);

%% -------------------------------------------------
%                     PLOT
%% -------------------------------------------------
figure('Color','w');
hold on; grid on;

% Real part curve (blue)
plot(theta1_degs, Zl_real, 'b-', 'LineWidth', 2);

% Imag part curve (red dashed)
plot(theta1_degs, Zl_imag, 'r--', 'LineWidth', 2);

xlabel('\theta_1 (graus)');
ylabel('Z_l');
title('Partes Real e Imaginária de Z_l = n_2 cos(\theta_2)');

legend('Re\{Z_l\}', 'Im\{Z_l\}', 'Location','Best');

%% Highlight the selected angles (Real + Imag) with complex labels

for i = 1:length(angles_mark)
    th = angles_mark(i);
    % index corresponding to this angle in the sampled vector
    idx = find(abs(theta1_degs - th) < 1e-6);

    % Real and imaginary values
    zr = Zl_real(idx);
    zi = Zl_imag(idx);

    % Plot markers on both curves
    plot(th, zr, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);  % real
    plot(th, zi, 'ks', 'MarkerFaceColor', 'c', 'MarkerSize', 8);  % imag

    % Build label string: Zl = a + bi
    txt = sprintf('Z_l(%.2f°) = %.4f + %.4fi', th, zr, zi);

    % Offset labels slightly to the right for readability
    text(th + 1, zr + 0.02*max(Zl_real), txt, ...
        'FontWeight','bold', 'Color','k');
end

% Labels
% text(0 + 1, Zl_mark_real(1), ' 0°', 'FontWeight','bold');
% text(31.31 + 1, Zl_mark_real(2), ' 31.31°', 'FontWeight','bold');
% text(90 - 12, Zl_mark_real(3), ' 90°', 'FontWeight','bold');

hold off;
