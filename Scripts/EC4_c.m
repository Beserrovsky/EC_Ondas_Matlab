% Poynting-z vs theta1 (incident, reflected, transmitted)
clear; clc; close all;

% -------------------------
% Media and wave parameters
% -------------------------
u0 = 4*pi*1e-7;
e0 = 8.854e-12;

u1 = u0;        e1 = 3.702*e0;    % medium 1 (dielectric)
u2 = u0;        e2 = e0;          % medium 2 (air)

% intrinsic impedances (SI)
eta1 = sqrt(u1 / e1);
eta2 = sqrt(u2 / e2);

% refractive indices (for Brewster/ Snell)
n1_refr = sqrt(e1 / e0);
n2_refr = sqrt(e2 / e0);

H0 = 1; % incident magnitude

% angle vector
theta1_deg = 0:0.1:90;
N = numel(theta1_deg);

% preallocate
Sz_inc = zeros(1,N);
Sz_ref = zeros(1,N);
Sz_tr  = zeros(1,N);

for k = 1:N
    th1 = deg2rad(theta1_deg(k));
    cos1 = cos(th1);

    % Snell
    sin_th2 = (n1_refr / n2_refr) * sin(th1);

    if abs(sin_th2) <= 1
        cos2 = sqrt(1 - sin_th2^2);
    else
        cos2 = 1i * sqrt(sin_th2^2 - 1); % TIR
    end

    % TE coef (Hy polarization)
    denom = eta2*cos1 + eta1*cos2;
    rho_TE = (eta2*cos1 - eta1*cos2) ./ denom;
    tH = (2*eta2*cos1) ./ denom;

    % Time-average Poynting components
    Sz_inc(k) = 0.5 * eta1 * cos1 * abs(H0)^2;
    Sz_ref(k) = 0.5 * eta1 * cos1 * abs(rho_TE)^2 * abs(H0)^2;
    Sz_tr(k)  = 0.5 * eta2 * real(cos2) * abs(tH)^2 * abs(H0)^2;
end

% Plot
figure('Color','w','Units','normalized','Position',[0.1 0.15 0.6 0.55]);
hold on; grid on; box on;

plot(theta1_deg, abs(Sz_inc), 'b-', 'LineWidth', 1.8);
plot(theta1_deg, abs(Sz_ref), 'r--', 'LineWidth', 1.8);
plot(theta1_deg, abs(Sz_tr),  'g-.', 'LineWidth', 1.8);

xlabel('$\theta_1\ (^\circ)$','Interpreter','latex','FontSize',14);
ylabel('$|S_z|\ (W/m^2)$','Interpreter','latex','FontSize',14);

title({'\bf Componentes do vetor de Poynting (TE)', ...
       '$S_z^{(1,+)}$, $S_z^{(1,-)}$, $S_z^{(2,+)}$'}, ...
       'Interpreter','latex','FontSize',15);

% create legend with LaTeX entries
lg = legend({'$|S_z^{(1,+)}|$', '$|S_z^{(1,-)}|$', '$|S_z^{(2,+)}|$'}, ...
    'Interpreter','latex', 'FontSize',14, 'Location','northeast');

% turn AutoUpdate off so new plots won't change the legend
lg.AutoUpdate = 'off';



% Brewster angle (TM)
theta_p_deg = 28;

% Points to highlight
angles_req = [0, theta_p_deg, 40];

for i = 1:length(angles_req)
    th = angles_req(i);
    [~, idx] = min(abs(theta1_deg - th));

    % Plot points
    plot(th, abs(Sz_inc(idx)), 'ko', 'MarkerFaceColor','b', 'MarkerSize',8);
    plot(th, abs(Sz_ref(idx)), 'ko', 'MarkerFaceColor','r', 'MarkerSize',8);
    plot(th, abs(Sz_tr(idx)),  'ko', 'MarkerFaceColor','g', 'MarkerSize',8);

    % LaTeX annotation
    txt = sprintf(['$\\theta_1 = %.3f^{\\circ}$\n' ...
                   '$|S_z^{(1,+)}| = %.3f$ W/m$^2$\n' ...
                   '$|S_z^{(1,-)}| = %.3f$ W/m$^2$\n' ...
                   '$|S_z^{(2,+)}| = %.3f$ W/m$^2$'], ...
                   theta1_deg(idx), ...
                   abs(Sz_inc(idx)), ...
                   abs(Sz_ref(idx)), ...
                   abs(Sz_tr(idx)));

    text(th + 0.8, ...
         max([abs(Sz_inc(idx)), abs(Sz_ref(idx)), abs(Sz_tr(idx))]) * 0.92, ...
         txt, ...
         'Interpreter','latex', ...
         'FontSize',11, ...
         'BackgroundColor','w', ...
         'EdgeColor','k');
end

% Print on console
fprintf('\n--- Numerical Results (LaTeX version) ---\n');
fprintf('Brewster angle (theta_p) = %.6f°\n\n', theta_p_deg);
for th = angles_req
    [~, idx] = min(abs(theta1_deg - th));
    fprintf(['theta1 = %7.3f°:\n' ...
             '   |S_inc| = %.3f W/m^2\n' ...
             '   |S_ref| = %.3f W/m^2\n' ...
             '   |S_tr | = %.3f W/m^2\n\n'], ...
             theta1_deg(idx), abs(Sz_inc(idx)), abs(Sz_ref(idx)), abs(Sz_tr(idx)));
end
