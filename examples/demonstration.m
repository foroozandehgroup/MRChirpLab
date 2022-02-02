%% Demonstration file on use of MRChirLab
% can be run section by section

clear all; close all;


%% Create pulses with MRchirp()
% creation of an HS1 pulse for inversion in NMR

% pulse parameters
par = struct('bw', 65e3, 'tp', 1000e-6, 'Q', 0.441, 'tres', 1e-6, ...
             'amp', 'sech','phase', 'tanh');

% create pulse
p1 = MRchirp(par);

% look at the pulse
plot_pulse(p1)


%% Manipulate pulses with functions from the MRchirp folder
% modify the pulse adiabaticity factor

p1 = pulse_modif(p1, "Q", 5);
plot_pulse(p1)


%% Simulate pulses with functions from the Magnetization folder
% simulate the effect of a single pulse

offsets = linspace(-p1.bw/2, p1.bw/2, 50);
magnetization = magn_calc_rot({p1}, p1.tp, offsets);
plot_magn(magnetization, offsets)


%% Use several pulses to make sequences
% double pulse for refocusing

% new position
par.delta_t = p1.delta_t + p1.tp;

% update the Q factor
par.Q = 5;

p2 = MRchirp(par);

plot_pulse(p2)
xlim([0 p2.t(end)])

% phase cycling
ph1 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph2 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];
CTP = [-2 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph1; ph2], CTP);

% simulate with added option to start with magnetization on +y 
% and phase cyle
opt.pc = pi/2 * [ph1; ph2; phrec];
opt.magn_init = repmat([0,1,0]', 1, 50);
magnetization = magn_calc_rot({p1, p2}, p2.delta_t+p2.tp/2, offsets, opt);
plot_magn(magnetization, offsets)


%% Create sequences directly (folders in Sequence folder)
% % the previous section would be equivalent to:
% par = struct('t180min', p1.tp, 'bw', p1.bw, ... 
%              'tres', 1e-6, 'display_result', true, ...
%              'pulse_param', struct('amp', 'sech','phase', 'tanh'));
% 
% ref = ref_2fs(par);

% example with compressed ABSTRUSE

% sequence parameters
par = struct('t90min', p1.tp/2, 't180min', p1.tp, 'bw', p1.bw, ... 
             'tres', 1e-6, ...
             'pulse_param', struct('amp', 'sech','phase', 'tanh'));

% sequence creation
abstruse = exc_3fs(par);

% plot the sequence
plot_seq(abstruse)


%% Use additional options from the sequence creation functions
% polynomial fitting phase correction

% new parameters
par.phase_polynomial_fitting = true;
par.polyfit_start = 10;
par.polyfit_stop = 10;
par.display_result = true;

tic
abstruse = exc_3fs(par);
toc


%% Manipulate sequences with functions from the Sequence folder

% add a pulse at the end of the sequence
abstruse = seq_add_pulse(abstruse, p1);
plot_seq(abstruse)


%% Export a pulse shape to be used on the instrument

% select current directory
path = [pwd '\'];

pulse_TopSpin_file(p1, 'demo_pulse', path)


%% Export a whole sequence

close all; % avoid figures export by seq_TopSpin_shape_files
seq_TopSpin_shape_files(abstruse, 'demo_seq', path)

