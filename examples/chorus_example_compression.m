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
param.w1max = 1.5451e4; % RF power max
param.TBPmin = 150;     % minimum time bandwidth product
param.bw = 300000;      % bandwidth
param.tres = 5e-7;      % time resolution

% sequence optional parameters
param.phase_polynomial_fitting = true;
param.polyfit_degree = 5;
param.display_result = true;

% sequence optional pulses parameters
param.pulse_param.n = 40; % superGaussian index

% chorus pulse sequence
tic
chorus1 = chorus(param);
toc


% -------------------------------------------------------------------------
% CHORUS defined with t90min, t180min, sinsmoothed pulses
% -------------------------------------------------------------------------

% sequence required parameters
param2.t90min = 500e-6;   % 90deg pulse minimum duration
param2.t180min = 1000e-6; % 180deg pulse minimum duration
param2.bw = 300000;
param2.tres = 5e-7;

% sequence optional param2eters
param2.phase_polynomial_fitting = true;
param2.polyfit_degree = 5;
param2.display_result = true;

% sequence optional pulses parameters
param2.pulse_param.type = "sinsmoothed";
param2.pulse_param.sm = 5; % smoothing percentage

% chorus pulse sequence
tic
chorus2 = chorus(param2);
toc






