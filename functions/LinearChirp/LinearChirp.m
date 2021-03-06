function p = LinearChirp(param)
% Returns a structure containing the different properties of a linear chirp
% pulse
%
% Input (properties in the structure param):
%   Required (3 of them):
%   - bw, bandwidth (Hz) - if negative, reversed sweep
%   - w1, excitation field value (Hz)
%   - tp, pulse duration (s)
%   - Q, adiabacitiy factor
%   Optional:
%   - delta_t, time offset or the pulse center position (s) 
%   - phi0, phase offset (rad)
%   - tres, time resolution of the pulse (s)
%   - type, pulse type - either superGaussian (default), sinsmoothed or
%   WURST
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
% sides are apodized with a sine function. A WURST pulse shape is also
% available.


grumble(param); % input check

if isfield(param, 'type')
    p.type = param.type;
else
    p.type = "superGaussian"; % default type
end

% reverse sweep
if isfield(param, 'bw')
    if param.bw < 0
        bw0 = param.bw; % negative value for calculation
        param.bw = -bw0;
    end
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

if exist('bw0','var') % reverse sweep case
    p.bw = bw0;
else
    p.bw = param.bw;
end
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
p.TBP = param.bw * tpa;

% number of points
p.np = floor(tpa / p.tres);

if p.np > 50000 % check to avoid overloading memory
    
    warning("Caution, high number of points (" + num2str(p.np)+ " points)")
    prompt = input('Keep going ? (Y to continue) ','s');
    
    if prompt ~= 'Y'
        error('Execution aborted')
    end
    disp("")
end

% adjusted tres
p.tres = tpa / p.np;

% time frame of the pulse: center of each segment (representation of an
% analytical function)
pulse_start = p.delta_t - p.tp / 2;
p.t = pulse_start + linspace(p.tres/2, p.tres*p.np-p.tres/2, p.np);

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
    t = p.t - p.delta_t + p.tp/2;
    integral_instant_freq = (sweep_rate * t.^2) / 2 + f0 * t;
    
    % polar coordinates Pr and Pph
    p.Pr = [smoothed_side unsmoothed_middle flip(smoothed_side)];
    p.Pph = 2 * pi * integral_instant_freq + p.phi0;

elseif p.type == "WURST"

        % smoothing percentage
    if isfield(param, 'n')
        p.n = param.n;
    else
        p.n = 80; % default n values
    end
    
    % phase calculated from instantaneous frequency
    % d(phase)/dt = sweep_rate * t + f0 (= instant. freq.)
    f0 = -p.bw / 2;
    sweep_rate = p.bw / p.tp;
    t = p.t - p.delta_t + p.tp/2;
    integral_instant_freq = (sweep_rate * t.^2) / 2 + f0 * t;
    
    % polar coordinates Pr and Pph
    p.Pr = p.w1 * (1-abs(sin((pi * (t-p.tp/2))/p.tp)).^p.n);
    p.Pph = 2 * pi * integral_instant_freq + p.phi0;

    % smoothing percentage sm
    try
        i_sm = 1; % unsmoothed part beginning index

        while p.Pr(i_sm) < 0.99 * p.w1
            i_sm = i_sm+1;
        end
        p.sm = 100 * i_sm / length(p.Pr);
    catch
        warning('Smoothing sm could not be computed')
    end  
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
    
    if isfield(param,'type')
        
        if param.type == "WURST"
            if ~isreal(param.n) || ~isreal(param.n) <= 0
                error(['n, the factor determining the smoothing enveloppe must be ' ...
                       'a positive real number for WURST pulse.'])
            end
        elseif param.type == "sinsmoothed"
            warning(['Index n not required for a sinsmoothed ' ...
                     'linear chirp, n neglected during creation.'])
            
        elseif param.type == "superGaussian"
            if ~isreal(param.n) || rem(param.n,2)~=0 || floor(param.n)~=param.n || param.n <= 0
                error(['n, the factor determining the smoothing enveloppe must be ' ...
                       'an even positive integer for a superGaussian pulse.'])
            end
        end
    elseif ~isreal(param.n) || rem(param.n,2)~=0 || floor(param.n)~=param.n || param.n <= 0
        error(['n, the factor determining the smoothing enveloppe must be ' ...
               'an even positive integer for a superGaussian pulse.'])
    end
end

if isfield(param,'type')
    if param.type ~= "superGaussian" && param.type ~= "sinsmoothed" && param.type ~= "WURST"
        error(['The type of the linear chirp must take one of the ' ...
               'following values: superGaussian, sinsmoothed, WURST.'])
    end
    
    if param.type == "sinsmoothed"
        if isfield(param, 'sm')
            if ~isreal(param.sm) || param.sm < 0 || param.sm > 100
                error(['sm, the smoothing percetnage must be a real ...'
                       'number between 0 and 100 '])
            end
        end
    end
end

end