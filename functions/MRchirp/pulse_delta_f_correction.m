function frequency_offset_pulse = pulse_delta_f_correction(p, delta_f_corr)
% Adds a frequency offset delta_f to a pulse p
%
% Input:
%   - p, pulse
%   - delta_f, frequency offset to be added
%
% Output:
%   - frequency_offset_pulse, pulse p with added frequency offset delta_f

p.Pph = p.Pph + 2 * pi * delta_f_corr * (p.t - p.delta_t);
p.delta_f = p.delta_f + delta_f_corr;

% recomputing cartesian coordinates
[p.Cx, p.Cy] = pol2cart(p.Pph, p.Pr);

frequency_offset_pulse = p;

end