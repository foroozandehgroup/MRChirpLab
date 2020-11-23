function compensated_pulse = pulse_resonator_easyspin(p, f, H_f, nu)
% [EPR] Makes the pulse compensate for the resonator effects with the
% resonator function from EasySpin
% 
% Input:
%   - p, pulse to be compensated
%   - f, frequency axis of the transfer function (GHz)
%   - H_f, transfer function of the resonator
%   - nu, center frequency of the resonator compensation (GHz)
%
% Output:
%   - compensated_pulse, pulse treated with the resonator function
%
% The transfer function of the resonator H_f is used to change the pulse 
% shape to compensate for resonator effects.
%
% Easyspin is required: http://www.easyspin.org


t = 1e6 * (p.t-min(p.t)); % time needs to be in us and starting from 0
y_t = p.Cx + 1i*p.Cy;

% extending resonator transfer function
f2 = linspace(min(f), max(f), 8*length(H_f));
[H_f2,index] = unique(H_f); % avoid repeated values
H_f2 = interp1(f(index), H_f2, f2, 'spline');

% resonator compensation
[t2, y_t2] = resonator(t, y_t, nu, f2, H_f2, 'compensate');

% interpolation of the compensated pulse on original time grid
y_t2 = interp1(t2, y_t2, t, 'nearest');

% recomputing pulse coordinates
p.Cx = real(y_t2);
p.Cy = imag(y_t2);
[p.Pr, p.Pph] = cart2pol(p.Cx, p.Cy);

compensated_pulse = p;

end

