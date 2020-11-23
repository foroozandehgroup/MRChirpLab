% JBV 18.11.2020
% Sequence pulses generation for broadband EPR at Q-band

clear all
close all

tic

% sequence required parameters
param.bw = 350e6;      % bandwidth
param.tres = 0.625e-9; % time resolution
% param.w1max = 37.5e6;
% param.TBPmin = 32;
param.t180min = 120e-9; % 180deg pulse minimum duration

% sequence optional parameters
param.display_result = true;
param.Q180 = 3;

% sequence optional pulses parameters
param.pulse_param.type = "sinsmoothed";
param.pulse_param.sm = 15;

DoubleChirp = doublechirp(param);

path = [pwd '\'];

% resonator profile information
resonator.f = load('resonator_profile_f_201001.mat').resonator_profile_f_201001;
resonator.H_f = load('resonator_profile_H_f_201001.mat').resonator_profile_H_f_201001;
resonator.nu = 33.8;

seq_Xepr_shape_files(DoubleChirp, '2chirps', '2313', path, resonator)
close all

toc