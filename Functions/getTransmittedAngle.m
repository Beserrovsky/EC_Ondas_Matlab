function theta2_deg = getTransmittedAngle(theta1_deg, e1, e2)
%GETTRANSMITTEDANGLE   Computes Snell's transmitted angle, including complex case.

arguments
    theta1_deg (1,1) double
    e1 (1,1) double
    e2 (1,1) double
end

theta1_rad = theta1_deg * pi/180;

snell = sin(theta1_rad) * sqrt(e1/e2);

if snell <= 1 && snell >= -1
    theta2_rad = asin(snell);               % real angle (radians)
else
    theta2_rad = (pi/2) + 1i * acosh(snell); % complex angle (radians)
end

theta2_deg = real(rad2deg(theta2_rad));         % MATLAB handles complex fine

% Ignores imaginary
end
