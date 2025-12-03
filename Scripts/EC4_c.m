% VARIABLES.
theta1_deg = 40;

% Medium Constraints
u0 = 4*pi*10^(-7);
e0 = 8.854*10^(-12);

u1 = u0;
e1 = 3.702*e0;

u2 = u0;
e2 = e0;

% Wave params
f = 10^9; % frequency [Hz]

% Incidency Angles
theta1_rad = deg2rad(theta1_deg);

theta2_deg = getTransmittedAngle(theta1_deg, e1, e2);
theta2_rad = deg2rad(theta2_deg);

% Electric and Magnetic Vectors
    % (SI units) x; y; z; --> Complex values in each slot.

H_y1 = 1e-3;

% !!! Onda incidente!!! 
%
% H_1Plus = [0; H_y1 * (1 - rho0); 0];       % A/m < -- Given with y+
% E_1Plus = [0; 0; 0];       % V/m
% N_1Plus = [0; 0; 0];       % W/m^2
%
H_1Plus = [0; H_y1 * (1 - rho0); 0];
E_1Plus = zeros(size(H_1Plus)); 
%N_1Plus = zeros(size(H_1Plus)); 

% !!! Onda refletida!!! 
%
% H_1Minus = [0; 0; 0];    % A/m
% E_1Minus = [0; 0; 0];    % V/m
% N_1Minus = [0; 0; 0];       % W/m^2
%
H_1Minus = zeros(size(H_1Plus)); 
E_1Minus = zeros(size(H_1Plus)); 


% !!! Onda transmitida!!! 
%
% H_2Plus = [0; 0; 0];    % A/m
% E_2Plus = [0; 0; 0];    % V/m
% N_2Plus = [0; 0; 0];       % W/m^2
%
H_2Plus = zeros(size(H_1Plus)); 
E_2Plus = zeros(size(H_1Plus)); 

% Get Eletric from Magnetic.
n1 = sqrt(u1/e1);
n2 = sqrt(u2/e2);

% SIMPLIFICATION:
% Everything for x = 0; z = 0;
%
% k1 = 2 * pi * f * sqrt(u1*e1);
% k2 = 2 * pi * f * sqrt(u2*e2);
% 
% beta_x1 = k1 * sin(theta1_rad);
% beta_z1 = k1 * cos(theta1_rad);

Zl = n2 * cos(theta2_rad); % has a complex angle.
Zz1 = n1 * cos(theta1_rad);

rho0 = (Zl - Zz1) / (Zl + Zz1);

% Ex_1Plus = n1 * cos(theta1_rad) * H_1Minus * exp^(-1i*beta_x1*x) * (exp^(-1i*beta_z1*z) + (rho0 * exp^(1i*beta_z1*z)));
    % Evertyhing for x = 0; z = 0;
    % Ex_1Plus = n1 * cos(theta1_rad) * H_1Minus * exp^0 * (exp^0 + (rho0 * exp^0));

% Incidente
E_1Plus(1) = n1 * cos(theta1_rad) * H_y1 * (1 + rho0); % Ex
E_1Plus(3) = -n1 * sin(theta1_rad) * H_y1 * (1 - rho0); % Ez

N_1Plus = cross(E_1Plus, H_1Plus);     % W/m^2

% Transmitida
H_2Plus(2) = H_1Plus(2);
E_2Plus(1) = n2 * cos(theta2_rad) * H_y1 * (1 - rho0); % Ex
E_2Plus(3) = -n2 * sin(theta2_rad) * H_y1 * (1 - rho0); % Ez

N_2Plus = cross(E_2Plus, H_2Plus);     % W/m^2

% Refletida - Errado.
H_1Minus(2) = H_1Plus(2);
E_1Minus(1) = n1 * cos(theta1_rad) * H_y1 * (-rho0); % Ex - Negativo
E_1Minus(3) = -n1 * sin(theta1_rad) * H_y1 * (-rho0); % Ez

N_1Minus = cross(E_1Minus, H_1Minus);     % W/m^2

% Poynting vectors (instantaneous for real-valued fields; for phasors use 0.5*real(cross(E,conj(H))))
S_1Plus  = cross(E_1Plus,  H_1Plus);    % incident
S_1Minus = cross(E_1Minus, H_1Minus);   % reflected
S_2Plus  = cross(E_2Plus,  H_2Plus);    % transmitted

% Prepare 3D plot
figure('Color','w','Name','E, H and Poynting Vectors');
hold on; grid on; axis equal;
view(38,22);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');

% Draw interface plane z = 0
[xp, yp] = meshgrid(linspace(-0.6,0.6,2), linspace(-0.6,0.6,2));
zp = zeros(size(xp));
surf(xp, yp, zp, 'FaceAlpha',0.15, 'EdgeColor','none', 'FaceColor',[0.8 0.8 1]);

% Origin
orig = [0;0;0];

% Plot vectors with scale factor for visibility
scale = 1; % adjust if vectors too small/large

quiver3(orig(1), orig(2), orig(3), real(E_1Plus(1))*scale, real(E_1Plus(2))*scale, real(E_1Plus(3))*scale, 'r', 'LineWidth',1.8, 'MaxHeadSize',0.6);
quiver3(orig(1), orig(2), orig(3), real(H_1Plus(1))*scale, real(H_1Plus(2))*scale, real(H_1Plus(3))*scale, 'm', 'LineWidth',1.8, 'MaxHeadSize',0.6);
quiver3(orig(1), orig(2), orig(3), real(S_1Plus(1))*scale, real(S_1Plus(2))*scale, real(S_1Plus(3))*scale, 'g', 'LineWidth',2, 'MaxHeadSize',0.8);

quiver3(orig(1), orig(2), orig(3), real(E_1Minus(1))*scale, real(E_1Minus(2))*scale, real(E_1Minus(3))*scale, 'r--', 'LineWidth',1.4, 'MaxHeadSize',0.6);
quiver3(orig(1), orig(2), orig(3), real(H_1Minus(1))*scale, real(H_1Minus(2))*scale, real(H_1Minus(3))*scale, 'm--', 'LineWidth',1.4, 'MaxHeadSize',0.6);
quiver3(orig(1), orig(2), orig(3), real(S_1Minus(1))*scale, real(S_1Minus(2))*scale, real(S_1Minus(3))*scale, 'g--', 'LineWidth',1.6, 'MaxHeadSize',0.8);

quiver3(orig(1), orig(2), orig(3), real(E_2Plus(1))*scale, real(E_2Plus(2))*scale, real(E_2Plus(3))*scale, 'c', 'LineWidth',1.8, 'MaxHeadSize',0.6);
quiver3(orig(1), orig(2), orig(3), real(H_2Plus(1))*scale, real(H_2Plus(2))*scale, real(H_2Plus(3))*scale, 'b', 'LineWidth',1.8, 'MaxHeadSize',0.6);
quiver3(orig(1), orig(2), orig(3), real(S_2Plus(1))*scale, real(S_2Plus(2))*scale, real(S_2Plus(3))*scale, 'k', 'LineWidth',2, 'MaxHeadSize',0.8);

% Annotations and legend entries with magnitudes
mE1 = vecnorm(real(E_1Plus).');
mH1 = vecnorm(real(H_1Plus).');
mS1 = vecnorm(real(S_1Plus).');

mE1m = vecnorm(real(E_1Minus).');
mH1m = vecnorm(real(H_1Minus).');
mS1m = vecnorm(real(S_1Minus).');

mE2 = vecnorm(real(E_2Plus).');
mH2 = vecnorm(real(H_2Plus).');
mS2 = vecnorm(real(S_2Plus).');

legend({
  sprintf('E_inc (|E|=%.3g V/m)', mE1), ...
  sprintf('H_inc (|H|=%.3g A/m)', mH1), ...
  sprintf('S_inc (|S|=%.3g W/m^2)', mS1), ...
  sprintf('E_ref (|E|=%.3g V/m)', mE1m), ...
  sprintf('H_ref (|H|=%.3g A/m)', mH1m), ...
  sprintf('S_ref (|S|=%.3g W/m^2)', mS1m), ...
  % sprintf('E_tr  (|  sprintf('H_tr  (|H|=%.3g A/m)', mH2), ...'
  sprintf('S_tr  (|S|=%.3g W/m^2)', mS2) }, 'Location','northeastoutside');

title(sprintf('Fields at interface (\\theta_1 = %g°). rho0 = %g%+gi', theta1_deg, real(rho0), imag(rho0)));
axis([-0.6 0.6 -0.6 0.6 -0.6 0.6]);

% Add small text labels near arrow tips
txtOffset = 0.02;
text(real(E_1Plus(1))+txtOffset, real(E_1Plus(2)), real(E_1Plus(3)), 'E_{inc}', 'Color','r');
text(real(H_1Plus(1))+txtOffset, real(H_1Plus(2)), real(H_1Plus(3)), 'H_{inc}', 'Color','m');
text(real(S_1Plus(1))+txtOffset, real(S_1Plus(2)), real(S_1Plus(3)), 'S_{inc}', 'Color','g');

text(real(E_1Minus(1))+txtOffset, real(E_1Minus(2)), real(E_1Minus(3)), 'E_{ref}', 'Color','r');
text(real(H_1Minus(1))+txtOffset, real(H_1Minus(2)), real(H_1Minus(3)), 'H_{ref}', 'Color','m');
text(real(S_1Minus(1))+txtOffset, real(S_1Minus(2)), real(S_1Minus(3)), 'S_{ref}', 'Color','g');

text(real(E_2Plus(1))+txtOffset, real(E_2Plus(2)), real(E_2Plus(3)), 'E_{tr}', 'Color','c');
text(real(H_2Plus(1))+txtOffset, real(H_2Plus(2)), real(H_2Plus(3)), 'H_{tr}', 'Color','b');
text(real(S_2Plus(1))+txtOffset, real(S_2Plus(2)), real(S_2Plus(3)), 'S_{tr}', 'Color','k');

hold off;
