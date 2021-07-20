% Example of pulse generation for EPR with mapping of the y-magnetization
% as a function time and offsets

clear all
close all

% sequence required parameters
param.t90min = 80e-9; % 90deg pulse minimum duration
param.t180min = 160e-9; % 180deg pulse minimum duration
param.t_delay = 90e-9;
param.bw = 500e6; % bandwidth
param.tres = 0.625e-9;

% sequence optional parameters
param.phase_polynomial_fitting = true;
param.polyfit_degree = 5;
param.display_result = true;

% sequence optional pulses parameters
param.pulse_param.amp = "sinsmoothed";
param.pulse_param.sm = 12.5;

% chorus pulse sequence
seq = chorus(param);

% phase cycling
ph1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ph2 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph3 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];

CTP = [1 -2 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph1; ph2; ph3], CTP);

ph_cy = pi/2 * [ph1; ph2; ph3; phrec];

% offsets and time axis
n_offs = 201;
offs = linspace(-seq.bw/2, seq.bw/2, n_offs);
t = 0:param.tres:800e-9;

% -------------------------------------------------------------------------
% CHORUS magnetization computation and plot
% -------------------------------------------------------------------------

% without phase cycling
final_magnetization = magn_calc_rot_full(seq.pulses, t, [0 0 0 0]', offs);

y_traj = transpose(squeeze(final_magnetization(2,:,:)));

figure('Position', [0 0 1000 400])
imagesc(1e9*t, 1e-6*offs, y_traj)
xlh = xlabel('$t$ (ns)','Interpreter','latex','FontSize',18);
ylabel('$\Omega$ / $2 \pi$ (MHz)','Interpreter','latex','FontSize',18);
h = gca;
set(h,'FontSize',18)
set(h, 'xtick', [0 200 400 600 800]);
set(h, 'ytick', [-250 0 250]);

% % with phase cycling
% final_magnetization2 = magn_calc_rot_full(seq.pulses, t, ph_cy, offs);
% 
% y_traj2 = transpose(squeeze(final_magnetization2(2,:,:)));
% 
% figure('Position', [0 0 1000 400])
% imagesc(1e9*t, 1e-6*offs, y_traj2)
% xlh = xlabel('$t$ (ns)','Interpreter','latex','FontSize',18);
% ylabel('$\Omega$ / $2 \pi$ (MHz)','Interpreter','latex','FontSize',18);
% h = gca;
% set(h,'FontSize',18)
% set(h, 'xtick', [0 200 400 600 800]);
% set(h, 'ytick', [-250 0 250]);












