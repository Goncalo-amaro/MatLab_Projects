function [Bm] = DipoleMagModel(P_v , t_p)

%% This function makes a simple approximation of the earth's magnetic field 
% using the dipole theory, with a coordinate changes a way of impletmenting 
% the magnetic tilted axis.

%% variables
 
b_0      = 8.0e15; % magnetic dipole moment (Tm^3)

theta    = deg2rad(169.74); % angle between the magnetic North pole and the geocentric North pole ( ?) 

alpha_0  = deg2rad(10); % angle between RAAN of the mangetic and geocentric models ( ?)

omega_0  = deg2rad(4.1667e-3); % Earth's rotation angular rate ( ?/s)

alpha_m  = alpha_0 + omega_0 * t_p; % magnetic dipole inertial rotation angle 

r_d      = norm(P_v);

r_v      = [ (P_v(1) / r_d) ; (P_v(2) / r_d) ; (P_v(3) / r_d) ]; % spacecraft position vector

N_v      = [ sin(theta) * cos(alpha_m);...
             sin(theta) * sin(alpha_m) ;...    % dipole unit vector
             cos(theta) ]; 



%% expression for the earth magnetic field model

Bm = (3 * dot(N_v,r_v) * r_v/r_d^2  - N_v); 

Bm = (b_0 / r_d^3) * Bm;

end











