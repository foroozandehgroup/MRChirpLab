% Example of pulse generation for EPR

clear all
close all

% sequence required parameters
param.t90min = 64e-9; % 90deg pulse minimum duration
param.t180min = 128e-9; % 180deg pulse minimum duration
param.bw = 400e6; % bandwidth
param.tres = 0.5e-9;

% sequence optional parameters
param.display_result = true;
param.t_delay = 50e-9;

% chorus pulse sequence
chorus = exc_3fs(param);