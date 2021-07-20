function seq = ref_2fs(param)
% Creates a refocusing pulse sequence with 2 FS pulses
% 
% Input:
%   - param, a structure containing the parameters allowing to defined the
%   sequence
%
% Output:
%   - seq, a structure containing the sequence information
%
% Required fields for param:
%   - bw, the target bandwidth of the pulse sequence (Hz)
%   - tres, the time resolution of the pulse sequence (s)
%   - one of the following combination):
%       - TBPmin, the minimum time bandwidth product and w1max, the maximum
%       amplitude of the excitation field B1
%       - t180min, the minimum pulse duration for a 90 degree pulse and a 
%       180 degree pulse respectively (s)
%
% Optional fields for param:
%   - Q180 the adiabaticity factor of the 180 degree pulse (by default 5)
%   - pulse_param, a structure containing desired LinearChirp parameters
%   - display_result, a boolean which allows to display the sequence and
%   the results of simulation/calculation (set to false by default)
%   - t_delay to add a delay between the 180deg pulses and at the start/end
%   of the sequence (s)
%
% Fields contained in seq:
%   - all the field mentionned above (with input/default value)
%   - tau, a vector containing the duration of the pulses and delays in
%   order (s)
%   - pulses, a cell array containing the pulse structures (LinearChirp)
%   - total_time, the total time of the pulse sequence (s)
%   - ph_cy, a proposed phase cycle - used for possible simulations in the
%   function


grumble(param)

% default values for optional parameters

if ~isfield(param, 'Q180')
    param.Q180 = 5;
end

if isfield(param, 'w1') && ~isfield(param, 'TBPmin')
    param.TBPmin = 100;
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
if isfield(param, 't180min')
    t180min = param.t180min;
    param.TBPmin = t180min * param.bw;
end

if isfield(param, 'w1max')
    t180min = w1max_TBP_compression(param.w1max, param.TBPmin, ...
                                    param.Q180, param.bw);
end

% rounding the pulse duration values
rounding_power = ceil(log10(param.tres));
t180min = 10^rounding_power * ceil(10^(-rounding_power) * t180min);

% potential delay at the end end of the sequence added to the delay between
% the 180 pulses
if isfield(param, 't_delay')
    t_delay = param.t_delay;
else
    t_delay = 0;
end

% sequence timing;
tau = [t180min t180min];

% pulse 1: pi/2 pulse
pulse_param.tp = tau(1);
pulse_param.delta_t = tau(1)/2;
pulse_param.Q = param.Q180;
p1 = MRchirp(pulse_param);

% pulse 2: pi pulse
pulse_param.tp = tau(2);
pulse_param.delta_t = tau(1) + tau(2)/2;
pulse_param.Q = param.Q180;
p2 = MRchirp(pulse_param);

% sequence
seq = param; % saving all the parameters
seq.tau = tau;
seq.pulses = {p1, p2};
seq.total_time = p2.delta_t + p2.tp / 2;

% suggested phase cycling
ph2 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph3 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];

CTP = [-2 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph2; ph3], CTP);

seq.pc = pi/2 * [ph2; ph3; phrec];

if isfield(param, 't_delay')
    seq = seq_add_delay(seq, t_delay/2, 3);
    seq = seq_add_delay(seq, t_delay, 2);
    seq = seq_add_delay(seq, t_delay/2, 1);
end
if param.display_result == true
    
    plot_seq(seq);
    seq_pulses_disp(seq);
    
    % offsets
    n_off = 101;
    off = linspace(-seq.bw/2, seq.bw/2, n_off);
    
    opt.magn_init = repmat([0,1,0]', 1, n_off);
    opt.pc = seq.pc;
    disp('Magnetization computation...')
    final_magn = magn_calc_rot(seq.pulses, seq.total_time, off, opt);
    
    if param.display_result == true
        plot_magn(final_magn, off)
    end
end

end

function grumble(param)

if ~isfield(param, 'tres')
    error('param must contain the time resolution tres (s)')
end

if ~isfield(param, 'bw')
    error('param must contain the bandwidth bw (Hz)')
elseif ~isreal(param.bw) || param.bw <= 0
    error('param.bw must be a real positive number')
end

% definiton check for t90min and t180min
if ~isfield(param, 't180min') && ~isfield(param, 'w1max')
    error('param must contain either t180min or w1max.')
elseif isfield(param, 't180min')
    if isfield(param, 'w1max')
        error('w1max already input, w1max cannot be input')
    end
    if ~isreal(param.t180min)
        error('t180min must be real')
    end
    if isfield(param, 'TBPmin')
        error('t90min and t180min already input, TBPmin cannot be input')
    end
elseif isfield(param, 'w1max')
    if isfield(param, 't180min')
        error('t180min already input, w1max cannot be input')
    end
    if ~isreal(param.w1max) || param.w1max <= 0
        error('w1max must be a real positive number')
    end
    if isfield(param, 'TBPmin')
        if ~isreal(param.TBPmin) || param.TBPmin <= 0
            error('TBPmin must be a real positive number')
        end
    end
end

if isfield(param, 'Q180')
    if ~isreal(param.Q180) || param.Q180 <= 0
        error('Q180 must be a real positive number')
    end
end

if isfield(param, 'display_result')
    if ~islogical(param.display_result)
        error('display_result must be a boolean')
    end
end

% checking for unexpected input
input_param = fieldnames(param);
for i = 1:length(input_param)
    if ~ismember(input_param{i},["bw", "tres", "TBPmin", "w1max", ...
                                 "t180min", "Q180", "pulse_param", ...
                                 "display_result", "t_delay"])
        warning(['Careful, ' input_param{i} ' is not a standard parameter.'])
    end
end


end