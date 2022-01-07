function pulse_sum_test()

%% partially overlapping pulses: Kunz-Bodenhausen

overlap = 450e-6;

par.tres = 0.5e-6;
par.bw = 300e3;
par.t90min = 2000e-6;
par.pulse_param.type = "superGaussian";

seq = exc_2fs(par);

plot_seq(seq)

seq.pulses{2} = pulse_modif(seq.pulses{2}, "delta_t", seq.pulses{2}.delta_t-overlap);

plot_seq(seq)

% summing the sequence pulses
p = pulse_sum(seq.pulses{1}, seq.pulses{2});

% sequence structure
seq1.total_time = p.tp + seq.tau(3) - overlap; % end delay
seq1.pulses = {p};

plot_seq(seq1)

% simulation

off = linspace(-150e3, 150e3, 100);

magn = magn_calc_rot(seq1.pulses, seq1.total_time, off);
plot_magn(magn,off);

%% fully overlapping pulses

p1 = pulse_modif(seq.pulses{2}, "delta_t", seq.pulses{2}.tp/2);
p2 = pulse_modif(p1, "delta_t", p1.tp + p1.tp/2);

p = pulse_sum(p1, p2);
p = pulse_sum(p, seq.pulses{1});

plot_pulse(p);

%% non overlapping pulses: chorus

par.t90min = 250e-6;
par.t180min = 500e-6;
par.phase_polynomial_fitting = true;

seq = exc_3fs(par);

% summing the sequence pulses
p = pulse_sum(seq.pulses{1}, seq.pulses{2});
p = pulse_sum(p, seq.pulses{3});

plot_pulse(p)

% simulation: should give the same result

magn = magn_calc_rot(seq.pulses, seq.total_time, off);
plot_magn(magn,off)

magn = magn_calc_rot({p}, p.tp, off);
plot_magn(magn,off)

end