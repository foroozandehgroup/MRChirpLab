% Illustration of different way of compressing a chorus pulse sequence
%
% The two sequences give equivalent results but are constructed with
% different parameters

close all
clear all


% -------------------------------------------------------------------------
% CHORUS with B1max/Rmin compresssion, superGaussian pulses
% -------------------------------------------------------------------------

% sequence required parameters
par.w1max = 1.5451e4; % RF power max
par.TBPmin = 150;     % minimum time bandwidth product
par.bw = 300000;      % bandwidth
par.tres = 5e-7;      % time resolution

% sequence optional parameters
par.phase_polynomial_fitting = true;
par.polyfit_degree = 5;
par.display_result = true;

% sequence optional pulses parameters
par.pulse_param.n = 40; % superGaussian index

% chorus pulse sequence
tic
chorus1 = chorus(par);
toc


% -------------------------------------------------------------------------
% CHORUS defined with t90min, t180min, sinsmoothed pulses
% -------------------------------------------------------------------------

% sequence required parameters
par2.t90min = 500e-6;   % 90deg pulse minimum duration
par2.t180min = 1000e-6; % 180deg pulse minimum duration
par2.bw = 300000;
par2.tres = 5e-7;

% sequence optional param2eters
par2.phase_polynomial_fitting = true;
par2.polyfit_degree = 5;
par2.display_result = true;

% sequence optional pulses parameters
par2.pulse_param.amp = "sinsmoothed";
par2.pulse_param.sm = 5; % smoothing percentage

% chorus pulse sequence
tic
chorus2 = chorus(par2);
toc

