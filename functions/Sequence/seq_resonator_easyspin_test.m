function seq_resonator_easyspin_test()

param.tres = 0.625e-9;
param.bw = 500e6;
param.t90min = 64e-9;
param.t180min = 64e-9;

chorus_uncompensated = chorus(param);

f = load('resonator_profile_f.mat').f;
H_f = load('resonator_profile_H_f.mat').H_f;
chorus_compensated = seq_resonator_easyspin(chorus_uncompensated, f, H_f, 9.55);

plot_seq(chorus_uncompensated)
plot_seq(chorus_compensated)

% phase offset for resonator compensated pulse
% pulses for phase cycling need to be created with the appropriate phase
% offset before correcting them
chorus_090 = pulse_phase_correction(chorus_compensated.pulses{1}, pi/2);

plot_pulse({chorus_090 chorus_compensated.pulses{1}})

chorus_simulated = seq_resonator_easyspin(chorus_compensated, f, H_f, 9.55, 'simulate');
plot_seq(chorus_simulated)

end