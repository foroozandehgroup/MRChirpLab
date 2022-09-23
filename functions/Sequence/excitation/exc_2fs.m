function seq = exc_2fs(param)
% Creates a double Frequency-Swept (FS) pulse sequence for excitation
% (Kunz scheme, Bohlen-Bodenhausen scheme)
% 
% Input:
%   - param, structure containing the parameters allowing to defined the
%   sequence
%
% Output:
%   - seq, structure containing the sequence information
%
% Required fields for param:
%   - bw, target bandwidth of the pulse sequence (Hz)
%   - tres, time resolution of the pulse sequence (s)
%   - one of the following:
%       - w1max, maximum amplitude of the excitation field B1
%       - t90min, minimum pulse duration for a 90 degree pulse(s)
%       - t180min, minimum pulse duration for a 180 degree pulse (s)
%
% Optional fields for parma:
%   - TBPmin, minimum time bandwidth product if w1max is used
%   - Q90 and Q180, adiabaticity factors of the 90 degree and 180
%   degrees pulses respectively   
%   - pulse_param, structure containing desired MRchirp parameters
%   - phase_polynomial_fitting, a boolean to launch a magnetizaiton
%   computation which leads to a phase correction of the sequence (set to
%   fault by default). When set to true, these additional parameters are
%   available (cf. polyfit_ph documentation for more details):
%       - polyfit_degree
%       - polyfit_start
%       - polyfit_stop
%   - display_result, boolean which allows to display the sequence and the 
%   results of simulation/calculation (set to false by default)
%
% Fields contained in seq:
%   - all the fields mentionned above (with input/default values)
%   - tau, vector containing the duration of the pulses and delays in
%   order (s)
%   - pulses, cell array containing the pulse structures (MRchirp)
%   - total_time, total time of the pulse sequence (s)
%   - pc, a proposed phase cycle - used for possible simulations in the
%   function


grumble(param)

% default values for optional parameters
if ~isfield(param, 'Q90')
    param.Q90 = (2 / pi) * log(2 / (cosd(90) + 1));
end

if ~isfield(param, 'Q180')
    param.Q180 = 5;
end

if isfield(param, 'w1max') && ~isfield(param, 'TBPmin')
    param.TBPmin = 100;
end

if ~isfield(param, 'phase_polynomial_fitting')
    param.phase_polynomial_fitting = false;
end

if ~isfield(param, 'display_result')
    param.display_result = false;
end

if isfield(param, 'pulse_param')
    pulse_param = param.pulse_param;
end

% common parameters for each pulse
pulse_param.bw = param.bw;
pulse_param.tres = param.tres;

% values for sequence required paramters
if isfield(param, 'w1max')      
    pulse_param.Q = param.Q90;
    t90min = w1max_TBPmin_MRChirp(pulse_param, param.w1max, param.TBPmin);
    pulse_param.Q = param.Q180;
    t180min = w1max_TBPmin_MRChirp(pulse_param, param.w1max, param.TBPmin);
elseif isfield(param, 't90min')
    t90min = param.t90min;
    t180min = 0;
elseif isfield(param, 't180min')
    t180min = param.t180min;
    t90min = 0;
end

% rounding the pulse duration values
rounding_power = ceil(log10(param.tres));
t90min = 10^rounding_power * ceil(10^(-rounding_power) * t90min);
t180min = 10^rounding_power * ceil(10^(-rounding_power) * t180min);

% defining the sequence with the most critical pulse
if t90min > 2*t180min
    disp("pi/2 defined") % high RF power and high TBPmin
    tau = [t90min 0.5*t90min 0.5*t90min];
    param.TBPmin = t90min * param.bw;
else    
    disp("pi defined")
    tau = [t180min*2 t180min t180min];
    param.TBPmin = t180min * param.bw;
end

% pulse 1: pi/2 pulse
pulse_param.tp = tau(1);
pulse_param.delta_t = tau(1)/2;
pulse_param.Q = param.Q90;
p1 = MRchirp(pulse_param);

% pulse 2: pi pulse
pulse_param.tp = tau(2);
pulse_param.delta_t = tau(1) + tau(2)/2;
pulse_param.Q = param.Q180;
p2 = MRchirp(pulse_param);

% double chirp pulse sequence
seq = param; % saving all the parameters
seq.tau = tau;
seq.pulses = {p1, p2};
seq.total_time = sum(tau);

% proposed phase cycling
ph1 = [0 0 0 0];
ph2 = [0 1 2 3];

CTP = [-1 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph1; ph2], CTP);

seq.pc = pi/2 * [ph1; ph2; phrec];

% magnetization calculation for polynomial fitting/display
if param.phase_polynomial_fitting == true || param.display_result == true
    
    % offsets
    off = linspace(-seq.bw/2+seq.pulses{1}.delta_f, ...
                   seq.bw/2+seq.pulses{1}.delta_f, 101);
    
    opt.pc = seq.pc;
    
    disp('Magnetization computation...')
    final_magn_1 = magn_calc_rot(seq.pulses, seq.total_time, off, opt);
    
    if param.display_result == true
        plot_magn(final_magn_1, off)
    end
end

%  polynomial fitting for phase correction
if param.phase_polynomial_fitting == true
    

    % polyfit options
    if ~isfield(param, 'polyfit_degree')
        param.polyfit_degree = 5;
    end
    
    if isfield(param, 'polyfit_start')
        polyfit_options.start = param.polyfit_start;
    end
    
    if isfield(param, 'polyfit_stop')
        polyfit_options.stop = param.polyfit_stop;
    end
    
    polyfit_options.polyfit_degree = param.polyfit_degree;
    polyfit_options.display_result = param.display_result;
    
    % phase retrieval
    ph = magn_phase(final_magn_1);
    
    % polyfit
    ph_corr = polyfit_ph(p1, ph, polyfit_options);

    % pulse 1 phase correction
    seq.pulses{1} = pulse_phase_correction(seq.pulses{1}, -ph_corr);
    
    if param.display_result == true
        
        disp('Magnetization computation...')
        final_magn_2 = magn_calc_rot(seq.pulses, seq.total_time, off, opt);
        plot_magn(final_magn_2, off)
        
    end
end

if param.display_result == true
    plot_seq(seq);
end

end

function grumble(param)

if ~isfield(param, 'tres')
    error('seq_param must contain the time resolution tres (s)')
end

if ~isfield(param, 'bw')
    error('seq_param must contain the bandwidth bw (Hz)')
elseif ~isreal(param.bw) || param.bw <= 0
    error('seq_param.bw must be a real positive number')
end

% definiton check for t90min and t180min
if sum([isfield(param, 't90min') isfield(param, 't180min') isfield(param, 'w1max')]) ~= 1
    error(['seq_param needs to contain one and one only of the' ...
           ' following parameters: t90min, t180min, w1max.'])
end

if isfield(param, 't90min')
    if ~isreal(param.t90min)
        error('t90min must be real')
    end
    if isfield(param, 'TBPmin')
        error('t90min and t180min already input, TBPmin cannot be input')
    end
    if param.t90min < param.tres
        error('t90min must be superior or equal to tres')
    end
elseif isfield(param, 't180min')
    if ~isreal(param.t180min)
        error('t180min must be real')
    end
    if isfield(param, 'TBPmin')
        error('t90min and t180min already input, TBPmin cannot be input')
    end
    if param.t180min < param.tres
        error('t180min must be superior or equal to tres')
    end 
elseif isfield(param, 'w1max')
    if ~isreal(param.w1max) || param.w1max <= 0
        error('w1max must be a real positive number')
    end
    if isfield(param, 'TBPmin')
        if ~isreal(param.TBPmin) || param.TBPmin <= 0
            error('TBPmin must be a real positive number')
        end
    end
end

if isfield(param, 'Q90')
    if ~isreal(param.Q90) || param.Q90 <= 0
        error('Q90 must be a real positive number')
    end
end

if isfield(param, 'Q180')
    if ~isreal(param.Q180) || param.Q180 <= 0
        error('Q180 must be a real positive number')
    end
end

if isfield(param, 'phase_polynomial_fitting')
    if ~islogical(param.phase_polynomial_fitting)
        error('phase_polynomial_fitting must be a boolean')
    end
end

if isfield(param, 'display_result')
    if ~islogical(param.display_result)
        error('display_result must be a boolean')
    end
end

if isfield(param, 'pulse_param')
    if ~isstruct(param.pulse_param)
        error('pulse_param must be a structure')
    end
end

% checking for unexpected input
input_param = fieldnames(param);
for i = 1:length(input_param)
    if ~ismember(input_param{i},["bw", "tres", "TBPmin", "w1max", ...
            "t90min", "t180min", "Q90" , "Q180", "pulse_param", ...
            "display_result", "phase_polynomial_fitting", ...
            "polyfit_degree", "polyfit_start", "polyfit_stop",])
        warning(['Careful, ' input_param{i} ' is not a standard parameter.'])
    end
end


end