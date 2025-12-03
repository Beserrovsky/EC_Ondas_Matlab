% Medium Constraints
u0 = 4*pi*10^(-7);
e0 = 8.854*10^(-12);

u1 = u0;
e1 = 3.702*e0;

u2 = u0;
e2 = e0;

% Create Angles - 0° --> 90°; 0.1° step
theta1_degs = 0:0.1:90;

% Find theta2 angles
theta2_degs = arrayfun(@(a) getTransmittedAngle(a, e1, e2), theta1_degs);

% Find Zl for every angle
Zl_out = arrayfun(@(a) getLoadImpedance(a, u2, e2), theta2_degs);

plot(theta2_degs, Zl_out)
