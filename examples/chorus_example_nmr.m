% recreating results from function CHORUS.m which can be found in
% supplementary information of :
% "Improved ultra-broadband chirp excitation." Foroozandeh M., Nilsson M.
% & Morris G. A. (2019). JMR doi:10.1016/j.jmr.2019.03.007
%
% equivalent to:
% [RF_exc,RF_ref1,RF_ref2,Ix,Iy,Iz,Phase]=CHORUS([500 1000],0.5,300000,5)

clear all; close all;

% sequence required parameters
param.t90min = 500e-6; % 500 90deg pulse minimum duration
param.t180min = 1000e-6; % 1000 180deg pulse minimum duration
param.bw = 300e3; % bandwidth % 
param.tres = 0.5e-6;

% sequence optional parameters
param.phase_polynomial_fitting = true;
param.polyfit_degree = 10;
param.display_result = true;

% sequence optional pulses parameters
param.pulse_param.type = "sinsmoothed";
param.pulse_param.sm = 10;

% chorus pulse sequence
tic
seq = chorus(param);
toc
