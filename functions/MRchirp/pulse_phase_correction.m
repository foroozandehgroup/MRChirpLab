function phased_pulse = pulse_phase_correction(p, ph_corr)
% Adds a phase ph_corr to a pulse p
%
% Input:
%   - p, pulse
%   - ph_corr, phase correction to be added. Can be input as one phase or a
%   vector of phase which will be applied to the phase poins of the pulse
%   after interpolation in case the two vectors present different number of
%   points.
%
% Output:
%   - phased_pulse, pulse p with added phase correction ph_corr


grumble(p, ph_corr)

if length(ph_corr) > 1
    if length(p.Pph) ~= length(ph_corr)

        x_p_Ph = linspace(1, length(p.Pph), p.np); 
        ph_corr = interp1(ph_corr, x_p_Ph);

    end
end

p.Pph = p.Pph + ph_corr;
p.phi0 = p.phi0 + ph_corr;

% recomputing cartesian coordinates
[p.Cx, p.Cy] = pol2cart(p.Pph, p.Pr);

phased_pulse = p;

end


function grumble(p, ph_corr)

if length(ph_corr) > 1
    if size(ph_corr, 2) ~= size(p.Pph, 2)
       error('2nd dimension missmatch between ph_corr and p.Pph')
    end
end

end