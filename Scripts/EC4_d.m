clear; clc;

%% PARAMETERS
theta1_deg = 60;
u0 = 4*pi*1e-7;
e0 = 8.854e-12;

u1 = u0; 
e1 = 3.702*e0;

u2 = u0;
e2 = e0;

f = 1e9; 
omega = 2*pi*f;

theta1 = deg2rad(theta1_deg);
n1 = sqrt(u1/e1);
n2 = sqrt(u2/e2);

%% Snell
theta2 = asin((n1/n2) * sin(theta1));

%% Waves
k1 = omega*sqrt(u1*e1);
k2 = omega*sqrt(u2*e2);

beta_x1 = k1*sin(theta1);
beta_z1 = k1*cos(theta1);
beta_z2 = k2*cos(theta2);

%% Reflection and transmission coefficients for H (TE)
rho0 = (n2*cos(theta1) - n1*cos(theta2)) / (n2*cos(theta1) + n1*cos(theta2));
tH   = (2*n2*cos(theta1)) / (n2*cos(theta1) + n1*cos(theta2));   % CORRECT ONE

%% Incident amplitude
H_y1 = 1;

%% z points
z_cm = -20:1:10;
z = z_cm/100;
Nz = numel(z);

Hy = zeros(Nz,1);

for iz = 1:Nz
    zi = z(iz);
    x = tand(theta1_deg) * zi;

    if zi < 0
        Hy(iz) = H_y1 * ( ...
              exp(-1i*(beta_x1*x + beta_z1*zi)) ...
            + rho0 * exp(-1i*(beta_x1*x - beta_z1*zi)) );
    elseif zi > 0
        Hy(iz) = tH * H_y1 * exp(-1i*(beta_x1*x + beta_z2*zi));
    else
        Hy(iz) = H_y1 * (1 + rho0); % total field at interface
    end
end

%% -----------------------------
% PLOT
% -----------------------------
figure; hold on; grid on;
plot(z_cm, abs(Hy), 'LineWidth', 2);
xlabel('z (cm)'); ylabel('|H_y(z)|');
title('|H_y| Amplitude do campo magnético na distância');
set(gca, 'FontSize', 12);

%% Highlight z=0 and z=5 cm
H0  = abs(Hy(z_cm == 0));
H5  = abs(Hy(z_cm == 5));

plot(0,  H0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
text(0+0.5, H0*1.05, sprintf('|H| at 0 cm = %.4f', H0), ...
    'Color','r','FontWeight','bold');

plot(5,  H5, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b');
text(5+0.5, H5*1.05, sprintf('|H| at 5 cm = %.4f', H5), ...
    'Color','b','FontWeight','bold');

hold off;

%% Display numerical results
fprintf('\n============================\n');
fprintf(' |H_y| at z = 0 cm  = %.6f\n', H0);
fprintf(' |H_y| at z = 5 cm  = %.6f\n', H5);
fprintf('============================\n\n');
