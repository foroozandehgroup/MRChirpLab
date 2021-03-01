function p = LinearChirp(param)
% Returns a structure containing the different properties of a linear chirp
% pulse
%
% Input (properties in the structure param):
%   Required (3 of them):
%   - bw, bandwidth (Hz)
%   - w1, excitation field value (Hz)
%   - tp, pulse duration (s)
%   - Q, adiabacitiy factor
%   Optional:
%   - delta_t, time offset or the pulse center position (s) 
%   - phi0, phase offset (rad)
%   - tres, time resolution of the pulse (s)
%   - type, pulse type - either superGaussian (default) or sinsmoothed
%   - delta_f, freqency offset (Hz) - not defined for sinsmoothed
%   - n, smoothing enveloppe exponent for superGaussian
%   - sm, smoothing percentage of the chirp - required for sinsmoothed
%
% Computed properties:
%   - one of bw, w1, tp or Q
%   - np, number of point
%   - TBP, time bandwidth product of the chirp
%   - t, time frame of the pulse (s),1*np array
%   - Pr and PPh, polar coordinates of the waveform, 1*np array
%   - Cx and Cy, cartesian coordinates of the waveform, 1*np array
%
% Output:
%   - p, a stucture with all the input and computed properties of the 
%   linear chirp pulse
%
% 3 parameters defintion can be made with bw, w1, tp. These can also be 
% configured by using 2 of them and a specified adiabaticity factor Q.
% 6 parameters definition can be made by adding delta_t, delta_f, phi0
%
% The pulse p is configured as a superGaussian by default. An option
% sinsmoothed is also availbale to create a rectangular shaped pulse whose
% sides are apodized with a sine function.


grumble(param); % input check

if isfield(param, 'type')
    p.type = param.type;
else
    p.type = "superGaussian"; % default type
end

% calculating missing parameter out of the input tp, w1, bw and Q
[calculated_param_value, calculated_param_name] = tp_w1_bw_Q(param);

if calculated_param_name == "tp"
    param.tp = calculated_param_value;
end

if calculated_param_name == "w1"
    param.w1 = calculated_param_value;
end

if calculated_param_name == "bw"
    param.bw = calculated_param_value;
end

if calculated_param_name == "Q"
    param.Q = calculated_param_value;
end

p.bw = param.bw;
p.w1 = param.w1;
p.tp = param.tp;
p.Q = param.Q;

% negative time equivalent to reversing the frequency sweep
tpa = abs(param.tp); % tpa used for computation

% 3 additional parametes for 6 parameters definition
if isfield(param,'delta_t')
    p.delta_t = param.delta_t;
else
    p.delta_t = tpa / 2;
end

% overall phase
if isfield(param,'phi0')
    p.phi0 = param.phi0;
else
    p.phi0 = 0;
end

% time resolution
if isfield(param,'tres')
    p.tres = param.tres;
else
    p.tres = 5e-7;
end

% TBP factor
p.TBP = p.bw * tpa;

% number of points
p.np = floor(tpa / p.tres);

if p.np > 50000 % check to avoid overloading memory
    
    disp("Caution, high number of points (" + num2str(p.np)+ " points)")
    prompt = input('Keep going ? (Y to continue) ','s');
    
    if prompt ~= 'Y'
        error('Execution aborted')
    end
    disp("")
end

% adjusted tres
p.tres = tpa / p.np;

% time frame of the pulse
pulse_start = p.delta_t - p.tp / 2;
p.t = pulse_start + linspace(0,p.tres*p.np,p.np);

if p.type == "superGaussian"
    
    % smoothing controlled by superGaussian index n
    if isfield(param,'n')
        p.n = param.n;
    else
        p.n = 40;
    end
    
    % delta_f only available for superGaussian type
    if isfield(param,'delta_f')
        p.delta_f = param.delta_f;
    else
        p.delta_f = 0;
    end

    % polar coordinates Pr and Pph
    p.Pr = p.w1 * exp( ...
             -(2^(p.n + 2)) * ((p.t - p.delta_t) / p.tp).^p.n);

    p.Pph = p.phi0 + ...
              pi * p.bw * (p.t - p.delta_t).^2 / p.tp + ...
              2 * pi * p.delta_f * (p.t - p.delta_t);
          
    % smoothing percentage sm
    i_sm = 1; % unsmoothed part beginning index
    while p.Pr(i_sm) < 0.99 * p.w1
        i_sm = i_sm+1;
    end
    p.sm = 100 * i_sm / length(p.Pr);

elseif p.type == "sinsmoothed"
    
    % smoothing percentage
    if isfield(param, 'sm')
        p.sm = param.sm;
    else
        p.sm = 10; % default sinsmoothed sm value
    end
    
    n_sm = floor((p.np * p.sm) / 100); % number of points smoothed
    n_unsm = p.np - (2 * n_sm); % number of points unsmoothed
    
    % amplitude apodized with a sine function taken from 0 to pi/2
    unsmoothed_middle = p.w1 * ones(1, n_unsm);
    smoothed_side = p.w1 * (sin(linspace(0, pi/2, n_sm)));

    % phase calculated from instantaneous frequency
    % d(phase)/dt = sweep_rate * t + f0 (= instant. freq.)
    f0 = -p.bw / 2;
    sweep_rate = p.bw / p.tp;
    t = linspace(0, p.tp, p.np);
    integral_instant_freq = (sweep_rate * t.^2) / 2 + f0 * t;
    
    % polar coordinates Pr and Pph
    p.Pr = [smoothed_side unsmoothed_middle flip(smoothed_side)];
    p.Pph = 2 * pi * integral_instant_freq + p.phi0;

end

% continuity of the phase
p.Pph = unwrap(wrapTo2Pi(p.Pph));

% cartesian coordinates
[p.Cx, p.Cy] = pol2cart(p.Pph, p.Pr);

end

function grumble(param)

if isfield(param,'delta_t')
    if ~isreal(param.delta_t) || param.delta_t < 0
    error('delta_t must be a positive real number')
    end
end

if isfield(param,'delta_f')
    if ~isreal(param.delta_f)
    error('delta_f must be a real number')
    end
end

if isfield(param,'phi0')
    if ~isreal(param.phi0)
    error('phi0 must be a real number')
    end
end

if isfield(param,'tres')
    if ~isreal(param.tres) || param.tres < 0
    error('tres must be a positive real number')
    end
end

if isfield(param,'n') 
    if ~isreal(param.n) || rem(param.n,2)~=0 || floor(param.n)~=param.n || param.n <= 0

    error(['n, the factor determining the smoothing enveloppe must be ' ...
          'an even positive integer.'])
    end
end

if isfield(param,'type')
    if param.type ~= "superGaussian" && param.type ~= "sinsmoothed"
        error(['The type of the linear chirp must take one of the ' ...
               'following values: superGaussian, sinsmoothed.'])
    end
    
    if param.type == "sinsmoothed"
        if isfield(param, 'sm')
            if ~isreal(param.sm) || param.sm < 0 || param.sm > 100
                error(['sm, the smoothing percetnage must be a real ...'
                       'number between 0 and 100 '])
            end
        end
        if isfield(param, 'n')
            warning(['Gaussian index n not required for a sinsmoothed ' ...
                     'linear chirp, n neglected during creation.'])
        end
    end
end

end