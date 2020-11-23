% JBV 18.11.2020
% Sequence pulses generation for article on use of linear chirp for 
% broadband excitation in EPR

clear all
close all

tic

% -------------------------------------------------------------------------
% common parameters for both sequences
% -------------------------------------------------------------------------

% sequence required parameters
param.bw = 400e6;      % bandwidth
param.tres = 0.625e-9; % time resolution
% param.w1max = 40e6;
% param.TBPmin = 32;

% sequence optional parameters
param.display_result = true;
param.Q180 = 4;

% sequence optional pulses parameters
param.pulse_param.type = "sinsmoothed";
param.pulse_param.sm = 12.5;
param.phase_polynomial_fitting = true;
param.polyfit_degree = 6;
param.display_results = true;

% resonator profile information
resonator.f = load('resonator_profile_f_200731.mat').f;
resonator.H_f = load('resonator_profile_H_f_200731.mat').H_f;
resonator.nu = 9.45;

path = [pwd '\'];

% -------------------------------------------------------------------------
% Kunz scheme without/with resonator compensation
% -------------------------------------------------------------------------

param.t180min = 160e-9; % 90deg pulse minimum duration

DoubleChirp = doublechirp(param);

% pulse shape Xepr files
seq_Xepr_shape_files(DoubleChirp, 'TwoChirp', '4201', path)
close all
seq_Xepr_shape_files(DoubleChirp, 'TwoChirp', '4209', path, resonator)
close all

% -------------------------------------------------------------------------
% CHORUS without/with resonator compensation
% -------------------------------------------------------------------------

param.t90min = 80e-9; % 90deg pulse minimum duration
param.t180min = 160e-9; % 180deg pulse minimum duration

TripleChirp = chorus(param);

% pulse shape Xepr files
seq_Xepr_shape_files(TripleChirp, 'CHORUS', '4325', path)
close all
seq_Xepr_shape_files(TripleChirp, 'CHORUS', '4337', path, resonator)
close all

toc