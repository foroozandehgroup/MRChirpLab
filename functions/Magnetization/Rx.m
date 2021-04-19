function rot_mat_x = Rx(phi)
% Returns the rotational matrix for an angle phi around the x-axis
%
% Input:
%   - phi, angle of rotation
%
% Output:
%   - rot_mat_x, rotational matrix around x-axis


rot_mat_x = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];

end
