function pulse_resonator_easyspin_test()

param.bw = 500e6;
param.tp = 100e-9;
param.Q = 5;
param.tres = 0.625e-9;

p = LinearChirp(param);

f = load('resonator_profile_f.mat').f;
H_f = load('resonator_profile_H_f.mat').H_f;
p_compensated = pulse_resonator_easyspin(p, f, H_f, 9.4);

pulses = [p, p_compensated];
titles = ["Uncompensated pulse", "Compensated pulse"];
plot_pulse(pulses, '',titles)
plot_pulse(pulses, 'polar',titles)


end