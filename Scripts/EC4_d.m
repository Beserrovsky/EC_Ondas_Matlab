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


% Get Eletric from Magnetic.
n1 = sqrt(u1/e1);
n2 = sqrt(u2/e2);

k1 = 2 * pi * f * sqrt(u1*e1);
k2 = 2 * pi * f * sqrt(u2*e2);

beta_x = k1 * sin(theta1_rad);
beta_z1 = k1 * cos(theta1_rad);
beta_z2 = k2 * cos(theta2_rad);

z = -20:1:10;
Nz = numel(z);

H3 = zeros(3, 1, Nz); 

H_z0 = [0; H_y1 * (1 - rho0); 0];

H = zeros(size(H_z0)); 

for iz = 1:Nz
    x = tand(theta1_deg) * z;

    % compute components depending on z(iz)
    Hx = 0;
    if (z < 0 )
        Hy = H_y1 * exp(-1i*beta_x1*x)*(exp(-1i*beta_z1*z) - rho0*exp(1i*beta_z1*z));
    end
    if (z > 0 )
        Hy = (1 - rho0) * H_y1 * exp(-1i*beta_x1*x)* exp(-1i*beta_z2*z);
    end
    if (z == 0)
        Hy = H_y1 * (1 - rho0);
    end
    Hz = 0;
    H3(:,1,iz) = [Hx; Hy; Hz];
end
