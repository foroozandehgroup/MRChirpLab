function rot_mat_z = Rz(phi)
% Returns the rotational matrix for an angle phi around the z-axis
%
% Input:
%   - phi, angle of rotation
%
% Output:
%   - rot_mat_z, rotational matrix around z-axis


rot_mat_z = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

end