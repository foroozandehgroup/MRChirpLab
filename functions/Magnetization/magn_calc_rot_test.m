function magn_calc_rot_test()

% Example of pulse generation for EPR

clear all
close all

% sequence required parameters
param.t90min = 64e-9; % 90deg pulse minimum duration
param.t180min = 128e-9; % 180deg pulse minimum duration
param.bw = 500e6; % bandwidth
param.tres = 0.625e-9;

% sequence optional parameters
param.display_result = false;

% chorus pulse sequence
seq = chorus(param);

% phase cycling
ph1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ph2 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph3 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];

CTP = [1 -2 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph1; ph2; ph3], CTP);

ph_cy = pi/2 * [ph1; ph2; ph3; phrec];

% offsets
n_offs = 101;
offs = linspace(-seq.bw/2, seq.bw/2, n_offs);

tic
opt.pc = ph_cy;
final_magn1 = magn_calc_rot(seq.pulses, seq.total_time, offs, opt);
toc
plot_magn(final_magn1, offs)


end