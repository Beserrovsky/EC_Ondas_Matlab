% Medium Constraints
% Dados (como no seu script)
u0 = 4*pi*1e-7;
e0 = 8.854e-12;
u1 = u0; e1 = 3.702*e0;
u2 = u0; e2 = e0;

theta1_degs = 0:0.1:90;

% Calcular theta2 em radianos (getTransmittedAngle deve devolver rad)
theta2_rads = arrayfun(@(a) getTransmittedAngle(a, e1, e2), theta1_degs);
theta2_degs = real(rad2deg(theta2_rads));

% Encontrar primeiro índice onde theta2 >= 90° (com tolerância)
tol = 1e-6;
idx = find(theta2_degs >= 90 - tol, 1, 'first');

% Plot básico
plot(theta1_degs, theta2_degs, 'b-', 'LineWidth',1.2)
xlabel('Ângulo da onda incidente  \theta_1 (degrees)')
ylabel('Ângulo da onda transmitida \theta_2 (degrees)')
title('Variação de ângulo na interface vidro-ar')
grid on
hold on

% Se encontrado, marcar e anotar
if ~isempty(idx)
    theta1_crit = theta1_degs(idx);
    theta2_at_crit = theta2_degs(idx);
    plot(theta1_crit, theta2_at_crit, 'ro', 'MarkerSize',8, 'LineWidth',1.5)
    txt = sprintf('\\theta_c = %.2f°', theta1_crit);
    text(theta1_crit + 1, theta2_at_crit - 5, txt, 'Color','r', 'FontWeight','bold')
else
    % opcional: informar que não atinge 90°
    disp('Nenhum theta1 no intervalo produz theta2 = 90° (dentro da tolerância).')
end

hold off
