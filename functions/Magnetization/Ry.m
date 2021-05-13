function rot_mat_y = Ry(phi)
% Returns the rotational matrix for an angle phi around the y-axis
%
% Input:
%   - phi, angle of rotation
%
% Output:
%   - rot_mat_y, rotational matrix around y-axis


rot_mat_y = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];

end