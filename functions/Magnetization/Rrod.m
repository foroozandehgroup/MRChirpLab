function total_rot_mat = Rrod(Cx_t, Cy_t, offs, delta_t)
% Returns the rotational matrix associated with a point of a linear chirp
% (or segment of the chirp)
% 
% Uses Rodrigues formula
%
% Input:
%     - Omega, point B1 field strength
%     - offs, point offset
%     - phi, point phase
%     - delta_t, point time (duration of the segment)
%
% Output:
%     - total_rot_mat: rotational matrix point associated with the 
%     pulse point
%


offs = offs + eps; % avoid singularity at offs = 0, Cx_t = 0 and Cy_t = 0

R = [0 -offs Cy_t; offs 0 -Cx_t;-Cy_t Cx_t 0] * delta_t;

gamma = delta_t * norm([Cx_t Cy_t offs]);

total_rot_mat = eye(3) + (sin(gamma)/gamma) * R + ((1-cos(gamma))/(gamma^2)) * R * R;


end