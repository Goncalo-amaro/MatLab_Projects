function [L_f,L_f_norm,L_f_tang] = LorentzForce ( q, m, R, V, B)

 omega_e = [0; 0; deg2rad(4.1667e-3)]; 
 V_rel = V - cross(omega_e,R);

 
% Lorentz force calculation

 L_f      = cross( ((q/m) * V_rel) , B); 
 
%  L_f_norm = L_f(1,1); % component that results from the normal between the magnetic field line and the radial vector ( x-axis )
%  
%  L_f_tang = L_f(2:3,1); % component that results from the tangent between the magnetic field line and the co-planar position components ( y & z axis )
 
 
end 