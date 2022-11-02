function seq = ref_3chirps_sum(param, param_sum)
% Creates an overlapped version of ref_3fs
% 
% Input:
%   - param, a structure containing the parameters allowing to defined the
%   sequence
%   - pulse_param, a structure containing the parameters of the overlap:
%
% Required fields for param:
%   - see help of non-overlapped sequence function
%
% Required fields for param_sum
%   - param_sum.overlap, the overlap duration (s)
%
% Optional fields for param_sum
%   - phase_ph_opt, a boolean to launch an optimization which leads to a 
%   phase correction of the sequence (set to faulse by default).
%   Optional fields:
%       - phase_opt_alpha, the weigth of the phase variance in the cost
%       function
%   - display_result, a boolean which allows to display the sequence and
%   the results of simulation/calculation (set to false by default)
%
% Fields contained in seq:
%   - all the field mentionned above (with input/default value)
%   - all the fields from the non-overlapped sequence (see help of 
%     non-overlapped sequence function)
%   - p, the waveform of the overlapped pulse

%% initialization with 3 frequency-swept pulses refocusing scheme

seq = ref_3fs(param);

%% overlap

seq.pulses{2} = pulse_modif(seq.pulses{2}, "delta_t", ...
                            seq.pulses{2}.delta_t-param_sum.overlap);
seq.pulses{3} = pulse_modif(seq.pulses{3}, "delta_t", ....
                            seq.pulses{3}.delta_t-param_sum.overlap*2);
                        
seq.total_time = seq.total_time - 2 * param_sum.overlap;

%% options

if isfield(param_sum, 'display_result')
    display_result = param_sum.display_result;
else
    display_result = true;
end

if isfield(param_sum, 'phase_opt')
    seq.phase_opt = param_sum.phase_opt;   
else
    seq.phase_opt = true;
end

%% 0th order phase correction with optimization
if seq.phase_opt == true
    
    % initial value of parameters
    X0 = [0];

    options = struct('Display', 'iter', 'HessUpdate', 'bfgs', ...
                     'GoalsExactAchieve', 1, 'GradConstr', true, ...
                     'TolX', 1e-8, 'TolFun', 1e-8, 'GradObj', 'off', ...
                     'MaxIter', 80, 'MaxFunEvals', 300*numel(X0)-1, ...
                     'DiffMaxChange', Inf, 'DiffMinChange', 1e-9, ...
                     'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, ...
                     'tau1',3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN',20);
        
    [x,~]=fminlbfgs(@mynestedfun, X0, options);

    seq.pulses{1} = pulse_phase_correction(seq.pulses{1}, x(1));
   
    seq.opt_x = x;
    
end

seq.p = pulse_sum(seq.pulses{1}, seq.pulses{2});
seq.p = pulse_sum(seq.p, seq.pulses{3});

if display_result == true
    
    plot_seq(seq)
    plot_pulse(seq.p)

    % offsets
    n_off = 101;
    off = linspace(-seq.bw/2, seq.bw/2, n_off);

    opt.magn_init = repmat([0,1,0]', 1, n_off);
    
    magn = magn_calc_rot_sum(seq.pulses, seq.total_time, off, opt);
    plot_magn(magn, off)
    
end

    function diff = mynestedfun(para)

    seqO = seq;

    % offsets
    n_offs = 51;
    offs = linspace(-seqO.bw/2, seqO.bw/2, n_offs);

    % correction
    seqO.pulses{1} = pulse_phase_correction(seqO.pulses{1}, para(1));
    
    % pulse sum
    pO = pulse_sum(seqO.pulses{1}, seqO.pulses{2});
    pO = pulse_sum(pO, seqO.pulses{3});
    
    % NB: faster than magn_calc_rot_sum when delays are in the sequence
    optO.magn_init = repmat([0,1,0]', 1, n_offs);
    magnO = magn_calc_rot({pO}, seqO.total_time, offs, optO);

    % only optimize on unsmoothed part

    % unsmoothed part beginning index for the magnetization phase vector
    if seqO.pulses{1}.amp == "sech"
        seqO.pulses{1}.sm = 5;
    end
    
    i_sm_magn = floor(length(magnO) * seqO.pulses{1}.sm /100);
    % beginning and end of smoothing
    begin_sm_magn = i_sm_magn;
    end_sm_magn = length(magnO) - i_sm_magn;
    % unsmoothed part selection
    magnO = magnO(:, begin_sm_magn:end_sm_magn);
    magnO = magnO / norm(magnO);

    % target - +y
    target = zeros(3,length(magnO(1,:)));
    target(2,:) = 1;
    target = target/norm(target);

    % cost function with variance of the phase
    % alpha: variance of the phase weight (set to 0 to avoid using it)
    diff = -trace(target'*magnO);

    end

end