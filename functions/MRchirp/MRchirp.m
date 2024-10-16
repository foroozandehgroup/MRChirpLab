function p = MRchirp(par)
% Returns a structure containing the different properties of a 
% Frequency-swept (FS) pulse
%
% Input (properties in the structure param):
%   Required:
%   - tres, time resolution of the pulse (s)
%   choose
%   - bw, bandwidth (Hz) - if negative, reversed sweep
%   - w1, excitation field value (Hz)
%   - tp, pulse duration (s) - can be replaced by B for tanh phase
%   or 2 of the above and:
%   - Q, adiabacitiy factor
%   or:
%   - [amp = "custom"] "P.r", amplitude polar coordinate of the pulse
%   - [phase = "custom"] "P.ph", phase polar coordinate of the pulse
%
%   Optional:
%   - delta_t, time offset or the pulse center position (s)
%   - phi0, phase offset (rad)
%   - delta_f, freqency offset (Hz)
%   - amp, string for the amplitude modulation function of the pulse:
%     "superGaussian" (default), "sinsmoothed", "WURST", "sech", "custom"
%   - phase, string for the phase modulation function of the pulse:
%     "chirp" (default), "tanh", "custom"
%   - [amp = "WURST"] n, smoothing enveloppe exponent
%     (= 20 by default)
%   - [amp = "superGaussian"] n, smoothing enveloppe exponent
%     (= 40 by default)
%   - [amp = "sinsmoothed"] sm, smoothing percentage of the chirp
%     (= 10 by default)
%   - [amp = "sech" or phase = "tanh"] B, smoothing parameter 
%     (= k / tp)
%   - [amp = "sech" or phase = "tanh"] k, smoothing parameter
%     (10.6 by default)
%
% Computed properties:
%   - one of bw, w1, tp or Q (except for custom pulse)
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
% Configure standard chirps with phase = "chirp" and 
% amp = "sinsmoothed", "superGaussian", "WURST"
% Configure standard Hyperbolic Sechant (HS) pulses with phase = "tanh"
% amp = "sech"
% 
% The pulse p is configured as a superGaussian chirp by default.


%% initialization 

grumble(par) % input check

% pulse structure
p = par;

% default amplitude and phase type
if ~isfield(p, 'amp') && ~isfield(p, 'phase')
    
    p.amp = "superGaussian";
    p.phase = "chirp";
    
elseif ~isfield(p, 'amp')
   
    if p.phase == "chirp"
       p.amp = "superGaussian";
   elseif p.phase == "tanh"
       p.amp = "sech";
   end
   
elseif ~isfield(p, 'phase')
    
   if p.amp == "sinsmoothed" || p.amp == "superGaussian" ...
		|| p.amp == "WURST" || p.amp == "linearsmoothed"
       p.phase = "chirp";
   elseif p.phase == "sech"
       p.phase = "tanh";
   end

end

%% required input: a combination of bw-w1-tp +Q(chirp)-k(HS)

% reverse sweep: negative bandwidth value par.bw
if isfield(p, 'bw')
   p.bw = abs(p.bw); % positive value required for computation
end

% calculating last parameter
if p.phase == "chirp" || p.phase == "superGaussian"
    
    if ~isfield(p, 'w1')
        p.w1 = sqrt(p.bw * p.Q / (2 * pi * p.tp));
    elseif ~isfield(p, 'tp')
        p.tp = p.bw * p.Q / (2 * pi * p.w1^2);
    elseif ~isfield(p, 'bw')
        p.bw = p.w1^2 * 2 * pi * p.tp / p.Q;
    elseif ~isfield(p, 'Q')
        p.Q = p.w1^2 * 2 * pi * p.tp / p.bw;
    end
    
elseif p.phase == "tanh"

    if ~isfield(p, 'k')
        p.k = 10.6;
    end
    
    if ~isfield(p, 'B') && isfield(par, 'tp')
        p.B = p.k/p.tp;
    end
    
    if ~isfield(p, 'w1')
        p.w1 = sqrt(p.bw * p.Q *p.B / (4 * pi));
    elseif ~isfield(p, 'B')
        p.B = 4 * pi * p.w1^2 / (p.Q * p.bw);
    elseif ~isfield(p, 'bw')
        p.bw = 4 * pi * p.w1^2 / (p.Q * p.B);
    elseif ~isfield(p, 'Q')
        p.Q = 4 * pi * p.w1^2 / (p.bw * p.B);
    end
    
    if ~isfield(p, 'tp')
        p.tp = p.k/p.B;
    end

    
elseif p.phase == "custom"
    
    % no additional parameter calculation
    
end

%% optional input default values 

% time offset (position)
if ~isfield(p,'delta_t')
    p.delta_t = p.tp / 2;
end

% frequency offset
if ~isfield(p,'delta_f')
    p.delta_f = 0;
end

% overall phase
if ~isfield(p,'phi0')
    p.phi0 = 0;
end

%% parameters computation: TBP, np, tres (adjusted), t

% TBP factor
if isfield(p, 'bw')
    p.TBP = p.bw * p.tp;
end

% number of points
p.np = floor(p.tp / p.tres);

% check to avoid overloading memory
if p.np > 50000 
    
    warning("Caution, high number of points (" + num2str(p.np)+ " points)")
    prompt = input('Keep going ? (Y to continue) ','s');
    
    if prompt ~= 'Y'
        error('Execution aborted')
    end
    disp("")
end

% adjusted tres
p.tres = p.tp / p.np;

% time frame of the pulse
% center of each segment (representation of an analytical function)
pulse_start = p.delta_t - p.tp / 2;
p.t = pulse_start + linspace(p.tres/2, p.tres*p.np-p.tres/2, p.np);

%% polar coordinate Pph: phase-modulation (frequency-modulation)

if p.phase == "chirp"
    
    % phase calculated from instantaneous frequency
    % d(phase)/dt = sweep_rate * t + f0 (= instant. phase.)
    % sweep_rate = p.bw / p.tp;     
    % instant_phase_integral = (sweep_rate * t.^2) / 2 + f0 * t;
    
    p.Pph = p.phi0 + ...
          pi * p.bw * (p.t - p.delta_t).^2 / p.tp + ...
          2 * pi * p.delta_f * (p.t - p.delta_t);
    
elseif p.phase == "tanh"
    
    p.Pph = 0;

    % phase calculated from instantaneous frequency integral
    % instant_freq = 0.5 * p.bw * tanh(p.B*t);
    % instant_phase_integral = 0.5 * p.bw * (1/B) * log(cosh(B * (t))); 

    p.Pph = p.phi0 + ...
            pi * p.bw * (1/p.B) * log(cosh(p.B * (p.t - p.delta_t))) +...
            2 * pi * p.delta_f * (p.t - p.delta_t);

elseif p.phase == "custom"
    
    % Pph input by the user
    
    if length(par.Pph) ~= p.np
        error("length(Pph) does not match tp and tres.")
    end
    
end

%% polar coordinates Ppr: amplitude-modulation function

if p.amp == "WURST"
    
    % smoothing factor n default value
    if ~isfield(p, 'n')
        p.n = 80;
    end
    
    p.Pr = p.w1 * ...
           (1 - abs(sin((pi * (p.t - p.delta_t)) / p.tp)).^p.n);

elseif p.amp == "superGaussian"
    
    % smoothing factor: superGaussian index n
    if ~isfield(p,'n')
        p.n = 40;
    end
        
    p.Pr = p.w1 * ...
           exp(-(2^(p.n + 2)) * ((p.t - p.delta_t) / p.tp).^p.n);

elseif p.amp == "sinsmoothed" || p.amp == "linearsmoothed"
    
    % default smoothing percentage value
    if ~isfield(p, 'sm')
        p.sm = 10; 
    end
    
    n_sm = floor((p.np * p.sm) / 100); % number of points smoothed
    n_unsm = p.np - (2 * n_sm); % number of points unsmoothed

    if p.amp == "sinsmoothed"
        % amplitude apodized with a sine function taken from 0 to pi/2
        unsmoothed_middle = p.w1 * ones(1, n_unsm);
        smoothed_side = p.w1 * (sin(linspace(0, pi/2, n_sm)));
    elseif p.amp == "linearsmoothed"
        % simple linear ramp
        unsmoothed_middle = p.w1 * ones(1, n_unsm);
        smoothed_side = p.w1 * linspace(0, 1, n_sm);
    end

    p.Pr = [smoothed_side unsmoothed_middle flip(smoothed_side)];

elseif p.amp == "sech"
    
    if ~isfield(p,'B')
        p.B = 10.6/p.tp;
    end
    
    p.Pr = p.w1 * sech(p.B*(p.t-p.delta_t));

elseif p.amp == "custom"
    
    % Pr input by the user
    
    if length(par.Pr) ~= p.np
        error("length(Pr) does not match tp and tres.")
    end
    
end

% estimation of smoothing percentage sm computed if sm not input
if ~isfield(par, 'sm')
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

%% Cartesian coordinates

% cartesian coordinates
[p.Cx, p.Cy] = pol2cart(p.Pph, p.Pr);

% reverse sweep case
if isfield(par, 'bw')
   if par.bw < 0
      p.Cy = -p.Cy;
      [p.Pph, p.Pr] = cart2pol(p.Cx, p.Cy);
      p.bw = -p.bw;
   end
end

% continuity of the phase
p.Pph = unwrap(wrapTo2Pi(p.Pph));

end

function grumble(par)

%% required input checks

if ~isfield(par,'tres')
    error('tres must be input')
elseif ~isreal(par.tres) || par.tres < 0
    error('tres must be a positive real number')
end

if isfield(par,'tp')
    if ~isreal(par.tp)
        error('tp must be a real number')
    end
    if par.tp < par.tres
        error('tp must be superior or equal to tres')
    end
end

if isfield(par,'w1')
    if ~isreal(par.w1) || par.w1 < 0
        error('w1 must be a positive real number')
    end
end

if isfield(par,'bw')
    if ~isreal(par.bw)
        error('bw must be a real number')
    end
end

if isfield(par, 'B')
    if ~isreal(par.B) || par.B <= 0
        error('B must be a positive real number for sech pulse.')
    end
end

%% pulse types input checks

if isfield(par, "phase")
    
    if par.phase == "chirp"

        if sum([isfield(par,'bw') isfield(par,'w1') ...
                isfield(par,'tp') isfield(par,'Q')]) ~= 3
            error(['3 of the following parameters must be defined: ' ...
                   'bw w1 tp Q'])
        end

        if isfield(par,'Q')
            if ~isreal(par.Q) || par.Q < 0
                error('Q must be a positive real number')
            end
        end
        
    elseif par.phase == "tanh"
       
       if isfield(par, 'tp') || isfield(par, 'B')
           if isfield(par, 'tp') && isfield(par, 'B')
               error('Only one of tp and B should be input.')
           end
           if sum([isfield(par,'bw') isfield(par,'w1') ...
                   isfield(par,'Q')]) ~= 2
                error(['Only 2 of the following parameters must be ' ...
                      'defined in addtion to tp or B: bw w1 Q'])
           end           
       elseif ~isfield(par,'bw') && ~isfield(par,'w1') && ~isfield(par,'Q')
            error(['Without tp or B, the following parameters ' ...
                   'must be defined: bw w1 Q'])
       end
       
       if isfield(par, 'k')
           if ~isreal(par.k) || par.k <= 0
               error('k must be a positive real number for sech pulse.')
           end
       end
       
    elseif par.phase == "custom"

        if ~isfield(par, 'Pph')
            error('Pph required for custom pulse')
        end
        
    else
        error(['phase must be a string which takes one of the ' ...
               'following value: chirp, superGaussian, tanh, custom'])
    end
    
    
end

if isfield(par, "amp")
    
    if par.amp == "WURST"

        if isfield(par, 'n')
            if ~isreal(par.n) || par.n <= 0
                error(['n, the factor determining the smoothing enveloppe must be ' ...
                       'a positive real number for WURST pulse.'])
            end
        end

    elseif par.amp == "superGaussian"

        if isfield(par, 'n')
            if ~isreal(par.n) || rem(par.n,2)~=0 || floor(par.n)~=par.n || par.n <= 0
                error(['n, the factor determining the smoothing enveloppe must ' ...
                       'be an even positive integer for a superGaussian pulse.'])
            end
        end

    elseif par.amp == "sinsmoothed" ||  par.amp == "linearsmoothed"

        if isfield(par, 'sm')
            if ~isreal(par.sm) || par.sm < 0 || par.sm > 100
                error(['sm, the smoothing percetnage must be a real ...'
                       'number between 0 and 100 '])
            end
        end

    elseif par.amp == "sech"

        if isfield(par, 'B')
            if ~isreal(par.B) || par.B <= 0
                error(['B, the factor determining the smoothing ' ...
                       'enveloppe must be a positive real number ' ...
                       'for sech pulse.'])
            end
        end

    elseif par.amp == "custom"

        if ~isfield(par, 'Pr')
            error('Pr required for custom pulse')
        end
        
    else
        error(['amp must be a string which must take one of the ' ...
               'following values: WURST, superGaussian, sinsmoothed, ' ...
               'linearsmoothed, sech, custom'])
    end
end

%% optional input checks

if isfield(par, 'delta_t')
    if ~isreal(par.delta_t) || par.delta_t < 0
        error('delta_t must be a positive real number')
    end
end

if isfield(par, 'delta_f')
    if ~isreal(par.delta_f)
        error('delta_f must be a real number')
    end
end

if isfield(par, 'phi0')
    if ~isreal(par.phi0)
        error('phi0 must be a real number')
    end
end

%% unexpected input checks

% checking for unexpected input
input_par = fieldnames(par);
for i = 1:length(input_par)
    
    if ~ismember(input_par{i}, ...
                 ["tres", "bw", "tp", "w1", "phi0", "delta_t", ...
                  "delta_f" , "amp", "phase", "Q", "sm", "n", "B", "k", ...
                  "Pph", "Pr"])
        
        warning(['Careful, ' input_par{i} ' is not a standard parameter.'])
    end
end

end