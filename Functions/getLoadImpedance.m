function Zl = getLoadImpedance(theta2_deg, u2, e2)
%GETTRANSMITTEDANGLE   Computes Impedance of TM signal

arguments
    theta2_deg (1,1) double
    u2 (1,1) double
    e2 (1,1) double
end

theta2 = theta2_deg * pi/180;

intrinsic = sqrt(u2/e2);

Zl = intrinsic * cos(theta2);

end
