function seq = chorus_sum(param, param_sum)


%% initialization with chorus

seq = exc_3fs(param);

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

if isfield(param_sum, 'phase_polynomial_fitting')
    seq.phase_polynomial_fitting = param_sum.phase_polynomial_fitting;
else
    seq.phase_polynomial_fitting = true;
end

if isfield(param_sum, 'phase_opt')
    
    seq.phase_opt = param_sum.phase_opt;
    
    if ~isfield(param_sum, 'phase_opt_alpha')
        seq.phase_opt_alpha = 1;
    end
    
else
    seq.phase_opt = true;
    seq.phase_opt_alpha = 1;
end

%% simulation
if display_result == true || seq.phase_polynomial_fitting == true
    
    off = linspace(-param.bw/2, param.bw/2, 201);
    opt.pc = seq.pc;
    magn = magn_calc_rot_sum(seq.pulses, seq.total_time, off, opt);
    
    if display_result == true
        plot_magn(magn,off)
    end
end

%% 1st phase correction with polynomial fitting
if seq.phase_polynomial_fitting == true
    
    if isfield(param_sum, 'polyfit_options')
        polyfit_options = param_sum.polyfit_options;
    else
        polyfit_options = struct();
    end
    
    % polyfit options
    if ~isfield(polyfit_options, 'polyfit_degree')
        param_sum.polyfit_degree = 10;
    end
    
    polyfit_options.display_result = param_sum.display_result;
    
    % phase retrieval
    ph = magn_phase(magn);
    
    % polyfit
    ph_corr = polyfit_ph(seq.pulses{1}, ph, polyfit_options);

    % pulse 1 phase correction
    seq.pulses{1} = pulse_phase_correction(seq.pulses{1}, ph_corr);
    
    if display_result == true
        magn = magn_calc_rot_sum(seq.pulses, seq.total_time, off, opt);
        plot_magn(magn, off)
    end
end

%% 2nd phase correction with optimization
if seq.phase_opt == true
    
    % initial value of parameters
    X0 = [0 0 0 0 0 0];

    options = struct('Display', 'iter', 'HessUpdate', 'bfgs', ...
                     'GoalsExactAchieve', 1, 'GradConstr', true, ...
                     'TolX', 1e-8, 'TolFun', 1e-8, 'GradObj', 'off', ...
                     'MaxIter', 80, 'MaxFunEvals', 300*numel(X0)-1, ...
                     'DiffMaxChange', Inf, 'DiffMinChange', 1e-9, ...
                     'OutputFcn', [], 'rho', 0.0100, 'sigma', 0.900, ...
                     'tau1',3, 'tau2', 0.1, 'tau3', 0.5, 'StoreN',20);
        
    [x,~]=fminlbfgs(@mynestedfun, X0, options);

    seq.pulses{1} = pulse_phase_correction(seq.pulses{1}, x(1));
    seq.pulses{2} = pulse_phase_correction(seq.pulses{2}, x(2));
    seq.pulses{3} = pulse_phase_correction(seq.pulses{3}, x(3));
    seq.pulses{1} = pulse_delta_f_correction(seq.pulses{1}, x(4));
    seq.pulses{2} = pulse_delta_f_correction(seq.pulses{2}, x(5));
    seq.pulses{3} = pulse_delta_f_correction(seq.pulses{3}, x(6));
    
    seq.opt_x = x;
    
end

seq.p = pulse_sum(seq.pulses{1}, seq.pulses{2});
seq.p = pulse_sum(seq.p, seq.pulses{3});

if display_result == true
    
    plot_seq(seq)
    plot_pulse(seq.p)
    
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
    seqO.pulses{2} = pulse_phase_correction(seqO.pulses{2}, para(2));
    seqO.pulses{3} = pulse_phase_correction(seqO.pulses{3}, para(3));
    seqO.pulses{1} = pulse_delta_f_correction(seqO.pulses{1}, para(4));
    seqO.pulses{2} = pulse_delta_f_correction(seqO.pulses{2}, para(5));
    seqO.pulses{3} = pulse_delta_f_correction(seqO.pulses{3}, para(6));
    
    % pulse sum
    pO = pulse_sum(seqO.pulses{1}, seqO.pulses{2});
    pO = pulse_sum(pO, seqO.pulses{3});
    
    % NB: faster than magn_calc_rot_sum when delays are in the sequence
    magnO = magn_calc_rot({pO}, seqO.total_time, offs);

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
    diff = -trace(target'*magnO) + ...
           seqO.phase_opt_alpha * var(abs(atan2(magnO(2,:),magnO(1,:))));

    end

end