function total_rot_mat = Rtot(Omega,offs,phi,time)
% Returns the rotational matrix associated with a point of a linear chirp
% 
% Inputs: pulse parameters
%     - Omega, point B1 field strength
%     - offs, point offset
%     - phi, point phase
%     - time, point time
%
% Output:
%     - total_rot_mat: rotational matrix point associated with the 
%     pulse point


% angular frequency of the effective field Beff
Omega_eff = sqrt((2 * pi * Omega)^2 + (2 * pi * offs)^2);

% angle between Beff and B1
theta = atan2(Omega, offs);

%flip angle alpha 
alpha = time * Omega_eff;

% total rotational matrix associated with the linear chirp point
total_rot_mat = Rz(phi) * Ry(theta) * Rz(alpha) * Ry(-theta) * Rz(-phi);

end