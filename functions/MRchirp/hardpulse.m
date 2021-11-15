function p = hardpulse(param)
% Generates a hard pulse of duration tp and flip angle alpha
%
% Input (properties in the structure param):
%   Required:
%   - tp, pulse duration (s)
%   - alpha, pulse flip angle (rad)
%   Optional:
%   - phi0, the phase of the pulse, 0 by default (rad)
%   - delta_t, the position of the pulse, tp/2 by default.
%
% Computed properties:
%   - np, number of points (2)
%   - tres, time resolution of the pulse (s)
%   - t, time frame of the pulse (s),1*np array
%   - Pr (Hz) and PPh (rad), waveform polar coordinates, 1*np arrays
%   - Cx (Hz) and Cy (Hz), waveform cartesian coordinates, 1*np arrays
%
% Output:
%   - p, a stucture with all the input and computed properties of the 
%   hard pulse as fields.

grumble(param)

% flip angle
p.alpha = param.alpha;

% duration
p.tp = param.tp;

% position
if isfield(param, 'delta_t')
    p.delta_t = param.delta_t;
else
    p.delta_t = p.tp / 2;
end

% amplitude
p.w1 = p.alpha/ p.tp;
p.w1 = p.w1/(2*pi); % rad to Hz

% time resolution
if isfield(param, 'tres')
    p.tres = param.tres;
else
    p.tres = p.tp;
end

% number of points
p.np = floor(p.tp / p.tres);

% adjusted tres
p.tres = p.tp / p.np;

% time frame of the pulse
pulse_start = p.delta_t - p.tp / 2;
p.t = pulse_start + linspace(p.tres/2, p.tres*p.np-p.tres/2, p.np);

% overall phase
if isfield(param, 'phi0')
    p.phi0 = param.phi0;
else
    p.phi0 = 0;
end

% amplitude vector
p.Pr = p.w1 * ones(1, length(p.t));

% phase vector
p.Pph = 0*p.Pr + p.phi0;

% cartesian coordinates
[p.Cx, p.Cy] = pol2cart(p.Pph, p.Pr);

end


function grumble(param)

if ~isfield(param, 'tp')
    error('param must contain the pulse duration tp (s)')
elseif ~isreal(param.tp) || param.tp < 0
    error('tp must be a positive real number')
end

if isfield(param,'delta_t')
    if ~isreal(param.delta_t) || param.delta_t < 0
    error('delta_t must be a positive real number')
    end
end

if isfield(param,'phi0')
    if ~isreal(param.phi0)
    error('phi0 must be a real number')
    end
end

end