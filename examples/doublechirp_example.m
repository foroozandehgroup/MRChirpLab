clear all;
close all;

param.bw = 300e3;
param.tres = 5e-7;
param.w1max = 15.451e3;
param.display_result = true;

seq = exc_2fs(param);