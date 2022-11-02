function seq = prochorus(param)
% Creates a perfect echo chorus (prochorus) pulse sequence
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
%   - one of the following combination:
%       - TBPmin and w1max, the minimum time bandwidth product and the
%       maximum amplitude of the excitation field B1
%       - t90min and t180min, the minimum pulse duration for a 90 degree 
%       pulse and a 180 degree pulse respectively (s)
%
% Optional fields for param:
%   - Q90 and Q180 the adiabaticity factors of the 90 degree and 180
%   degrees pulses respectively   
%   - pulse_param, a structure containing desired MRchirp parameters
%   - compression, a boolen to use the compressed version of prochorus
%   - display_result, a boolean which allows to display the sequence and
%   the results of simulation/calculation (set to false by default)
%   - coeff,a custom set of coefficients to adjust the length of the pulses
%   - phase_ph_opt, a boolean to launch an optimization which leads to a 
%   phase correction of the sequence (set to faulse by default).
%   Optional fields:
%       - ph_opt_parnb, the number of parameters to optimize - can be either
%       4 for optimization with the 3rd and 6th pulse or 8 to add the 1st
%       and 4th pulse (set to 4 by default)
%       - ph_opt_alpha, the weigth of the phase variance in the cost
%       function
%       - ph_opt_x, used as optimization solution (no optimization run)
%
% Fields contained in seq:
%   - all the field mentionned above (with input/default value)
%   - tau, a vector containing the duration of the pulses and delays in
%   order (s)
%   - pulses, a structure containing the pulse structures (MRchirp)
%   - total_time, the total time of the pulse sequence (s)
%
% ! Phase correction optimization not satifying?
% Try:
%   - using the phase variance by setting ph_opt_alpha = 1
%   - varying the weight of the phase variance ph_opt_alpha
%   - increasing the number of parameters to optimize (will make the
%   optimization slower)
% Other solutions (at the expense of longer duration/lower bandwidth) 
% include:
%   - increasing the TBPmin (only give minor improvements > 150 for chirps)
%   - increasing the Q factor of the 180deg Q180 (only give minor
%   improvements > 5 for chirps)
%   - increasing the smoothing factor of the pulse sm (loss of effective
%   bandwidth)

grumble(param)

%%  default values for optional parameters

if ~isfield(param, 'Q90')
    param.Q90 = (2 / pi) * log(2 / (cosd(90) + 1));
end

if ~isfield(param, 'Q180')
    param.Q180 = 5;
end

if isfield(param, 'w1') && ~isfield(param, 'TBPmin')
    param.TBPmin = 100;
end

if isfield(param, 'pulse_param')
    pulse_param = param.pulse_param;
end

% possible compression or use of custom coefficients
if isfield(param, 'coeff')
    coeff = param.coeff;
elseif isfield(param,'compression') && param.compression == true
    coeff = [1 0 1 0 1 0.5 0.5 0; 0 1 0 1 0 1 0 1]';
else
    % uncompressed (also default value if no compression or coeff field)
    coeff = [2 1/2 3/2 0 1 3/2 1/2 1; 0 1 0 1 0 1 0 1]';
end

% possible phase correction parameters
if isfield(param, 'phase_ph_opt') && param.phase_ph_opt == true
    if ~isfield(param, 'ph_opt_parnb')
        param.ph_opt_parnb = 4;
    end
    if ~isfield(param, 'ph_opt_alpha')
        param.ph_opt_alpha = 0;
    end
else    
    param.phase_ph_opt = false;
end

if ~isfield(param, 'display_result')
    param.display_result = false;
end

%% values for sequence required paramters

if isfield(param, 't90min') && isfield(param, 't180min')
    t90min = param.t90min;
    t180min = param.t180min;

    param.TBPmin = t90min * param.bw;
end

if isfield(param, 'w1max')
    pulse_param.Q = param.Q90;
    t90min = w1max_TBPmin_MRChirp(pulse_param, param.w1max, param.TBPmin);
    pulse_param.Q = param.Q180;
    t180min = w1max_TBPmin_MRChirp(pulse_param, param.w1max, param.TBPmin);
end

% rounding the pulse duration values
rounding_power = ceil(log10(param.tres));
t90min = 10^rounding_power * ceil(10^(-rounding_power) * t90min);
t180min = 10^rounding_power * ceil(10^(-rounding_power) * t180min);

%% sequence structure inialization
% excitation - inversion into 90y - refocusing
% 90 -180-delay-180-90y -180-delay-180

seq = param; % saving all parameters

% adiabaticity factor (Q = 0 indicates a delay)
Q90 = param.Q90;
Q180 = param.Q180;
seq.Q = [Q90 Q180 0 Q180 Q90 Q180 0 Q180];
   
seq.tau = coeff * [t90min t180min]'; % sequence timing
seq.tau = seq.tau';

seq.total_time = sum(seq.tau);

%% sequence pulses

% common parameters for each pulse
pulse_param.bw = param.bw;
pulse_param.tres = param.tres;

% generating the sequence pulses into a structure
ip = 1;
for i = 1:length(seq.tau)
    
    if seq.Q(i) > 0 % indicates presence of a pulse
        
        pulse_param.Q = seq.Q(i);
        pulse_param.tp = seq.tau(i);
        pulse_param.delta_t = seq.tau(i)/2; % pulses positions
        
        if i>1
            pulse_param.delta_t = pulse_param.delta_t + sum(seq.tau(1:i-1));
        end
        
        seq.pulses{ip} = MRchirp(pulse_param);
        ip = ip+1;
    end
end

% phase shift of the 2nd 90deg pulse
seq.pulses{4} = pulse_phase_correction(seq.pulses{4}, pi/2);

%% possible phase correction with superGaussian pulses

if param.phase_ph_opt == true
   
    if isfield(param, 'ph_opt_x')
        seq.ph_opt_x = param.ph_opt_x;
    else
        % initial values of parameters    
        if param.ph_opt_parnb == 4
            X0 = [0 0 0 0];
        elseif param.ph_opt_parnb == 8
            X0 = [0 0 0 0 0 0 0 0];
        end
        
        % optimizer options
        options = struct('Display', 'iter', 'HessUpdate', 'bfgs', ...
                     'GoalsExactAchieve', 1, 'GradConstr', true, ...
                     'TolX', 1e-8, 'TolFun', 1e-8, 'GradObj', 'off', ...
                     'MaxIter', 80, 'MaxFunEvals', 300*numel(X0)-1, ...
                     'DiffMaxChange', Inf, 'DiffMinChange', 1e-9, ...
                     'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, ...
                     'tau1',3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN',20);
        
        % optimizer
        [seq.ph_opt_x,~] = fminlbfgs(@mynestedfun, X0, options); 
    end
    
    % pulse correction with optimization results 
    x = seq.ph_opt_x;
    
    if length(x) == 4
        seq.pulses{3} = pulse_phase_correction(seq.pulses{3}, x(1));
        seq.pulses{6} = pulse_phase_correction(seq.pulses{6}, x(2));
        seq.pulses{3} = pulse_delta_f_correction(seq.pulses{3}, x(3));
        seq.pulses{6} = pulse_delta_f_correction(seq.pulses{6}, x(4));
    elseif length(x) == 8
        seq.pulses{1} = pulse_phase_correction(seq.pulses{1}, x(1));
        seq.pulses{3} = pulse_phase_correction(seq.pulses{3}, x(2));
        seq.pulses{4} = pulse_phase_correction(seq.pulses{4}, x(3));
        seq.pulses{6} = pulse_phase_correction(seq.pulses{6}, x(4));
        seq.pulses{1} = pulse_delta_f_correction(seq.pulses{1}, x(5));
        seq.pulses{3} = pulse_delta_f_correction(seq.pulses{3}, x(6));
        seq.pulses{4} = pulse_delta_f_correction(seq.pulses{4}, x(7));
        seq.pulses{6} = pulse_delta_f_correction(seq.pulses{6}, x(8));
    end

end

%% proposed phase cycle

% 64 steps phase cycling for 4 pi pulses
ph1 = zeros(1, 64);
ph2 = repmat([0 1 2 3],1,16);
ph3 = repmat(repelem([0 1 2 3],4),1,4);
ph4 = zeros(1, 64);
ph5 = repmat(repelem([0 1 2 3],8),1,2);
ph6 = repelem([0 1 2 3],16);

CTP = [+1 -2 +2 0 -2 +2]; % coherence transfer pathway

phrec = phase_cycle_receiver([ph1; ph2; ph3; ph4; ph5; ph6], CTP);

seq.pc = pi/2 * [ph1; ph2; ph3; ph4; ph5; ph6; phrec];

%% add potential delay

if isfield(param, 't_delay')
    
    seq = seq_add_delay(seq, seq.t_delay/2,9);
    seq = seq_add_delay(seq, seq.t_delay,7);
    seq = seq_add_delay(seq, seq.t_delay/2,6);
    
    seq = seq_add_delay(seq, seq.t_delay/2,5);
    seq = seq_add_delay(seq, seq.t_delay,3);
    seq = seq_add_delay(seq, seq.t_delay/2,2);
    
end

if param.display_result == true
    seq_pulses_disp(seq);
end

%% possible display of the sequence and its result (simulation)

if param.display_result == true
    
    plot_seq(seq)
    
    off = linspace(-seq.bw/2, seq.bw/2, 101);
    
    opt.pc = seq.pc;
    final_magn = magn_calc_rot(seq.pulses, seq.total_time, off, opt);
    
    plot_magn(final_magn, off)
    
end

%% optimization function

    function diff = mynestedfun(par)
  
    seqO = seq; % sequence to be modified for optimization

    % offsets
    offO = linspace(-seqO.bw/2, seqO.bw/2, 50);

    % set up tested parameters
    if param.ph_opt_parnb == 4
        seqO.pulses{3} = pulse_phase_correction(seqO.pulses{3}, par(1));
        seqO.pulses{6} = pulse_phase_correction(seqO.pulses{6}, par(2));
        seqO.pulses{3} = pulse_delta_f_correction(seqO.pulses{3}, par(3));
        seqO.pulses{6} = pulse_delta_f_correction(seqO.pulses{6}, par(4));
    elseif param.ph_opt_parnb == 8
        seqO.pulses{1} = pulse_phase_correction(seqO.pulses{1}, par(1));
        seqO.pulses{3} = pulse_phase_correction(seqO.pulses{3}, par(2));
        seqO.pulses{4} = pulse_phase_correction(seqO.pulses{4}, par(3));
        seqO.pulses{6} = pulse_phase_correction(seqO.pulses{6}, par(4));
        seqO.pulses{1} = pulse_delta_f_correction(seqO.pulses{1}, par(5));
        seqO.pulses{3} = pulse_delta_f_correction(seqO.pulses{3}, par(6));
        seqO.pulses{4} = pulse_delta_f_correction(seqO.pulses{4}, par(7));
        seqO.pulses{6} = pulse_delta_f_correction(seqO.pulses{6}, par(8));
    end
    
    % magnetization simulation
    magn = magn_calc_rot(seqO.pulses, seqO.total_time, offO);

    % target unsmoothed part for optimization

    % unsmoothed part beginning index for the magnetization phase vector
    i_sm_magn = floor(length(magn) * seqO.pulses{1}.sm/100); % 5 /100
    % beginning and end of smoothing
    begin_sm_magn = i_sm_magn;
    end_sm_magn = length(magn) - i_sm_magn;
    % unsmoothed part selection
    magn = magn(:, begin_sm_magn:end_sm_magn);
    magn = magn / norm(magn);

    % target - +y
    target = zeros(3,length(magn(1,:)));
    target(2,:) = 1;
    target = target/norm(target);

    % cost function with variance of the phase
    % alpha: variance of the phase weight (set to 0 to avoid using it)
    diff = -trace(target'*magn) + seq.ph_opt_alpha * ...
           var(abs(atan2(magn(2,:),magn(1,:))));

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
if ~isfield(param, 't90min') && ~isfield(param, 'w1max')
    error('seq_param must contain either t90min and t180min, or w1max.')
elseif isfield(param, 't90min') || isfield(param, 't180min')
    if ~isfield(param, 't90min')
        error('seq_param must contain t90min if t180min is input')
    end
    if ~isfield(param, 't180min')
        error('seq_param must contain t180min if t90min is input')
    end
    if ~isreal(param.t90min)
        error('t90min must be real')
    end
    if ~isreal(param.t180min)
        error('t180min must be real')
    end
    if isfield(param, 'w1max')
        error('t90min and t180min already input, w1max cannot be input')
    end
    if isfield(param, 'TBPmin')
        error('t90min and t180min already input, TBPmin cannot be input')
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

if isfield(param, 'phase_ph_opt')
    if ~islogical(param.phase_ph_opt)
        error('phase_ph_opt must be a boolean')
    end
else
    if isfield(param, "ph_opt_x")
        error("ph_opt_x unused, please set phase_ph_opt to true.")
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
            "t90min", "t180min", "Q90" , "Q180", "pulse_param", ...
            "compression", "display_result" , "coeff", "phase_ph_opt", ...
            "ph_opt_parnb", "ph_opt_alpha", "ph_opt_x", "t_delay"])
        disp(['Careful, ' input_param{i} ' is not a standard parameter.'])
    end
end

    
end
