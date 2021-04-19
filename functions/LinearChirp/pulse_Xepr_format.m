function shape = pulse_Xepr_format(pulse)
% Creates a matrix in the format of pulses used by Xepr from a pulse
% structure
%
% Input:
%   - pulse, pulse to extract the shape from
%
% Ouput
%   - shape, a matrix witht the shape of the pulse in Xepr format
%   (cartesian coordinates normalized from -1 to 1)


shape(:,1) = rescale(pulse.Cx, 'maxabs');
shape(:,2) = rescale(pulse.Cy, 'maxabs');

% only working from Matlab 2020a
% shape(:,1) = normalize(pulse.Cx, 'range', [-1 1]);
% shape(:,2) = normalize(pulse.Cy, 'range', [-1 1]);

end