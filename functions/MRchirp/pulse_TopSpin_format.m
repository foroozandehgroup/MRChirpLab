function shape = pulse_TopSpin_format(p)
% Creates a matrix in the format of pulses used by Bruker from a pulse
% structure
%
% Input:
%   - p, pulse to extract the shape from
%
% Ouput
%   - shape, a matrix witht the shape of the pulse in TopSpin format
%   (polar coordinates with r normalized from 0 to 100 and phi from
%	0 to 360°)


shape(:,1) = 100 * p.Pr / max(p.Pr);

% does not work with Matlab 2019
% shape(:,1) = 100 * rescale(p.Pr, 'minmax');

shape(:,2) = rad2deg(wrapTo2Pi(p.Pph));

end