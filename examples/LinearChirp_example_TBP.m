% showing compression of chirps with different time bandwidth products

clear all
close all

param.tres = 5e-7;
param.bw = 300e3;
param.Q = (2 / pi) * log(2 / (cosd(90) + 1));

w1_max = 15.45e3;
TBPmin = 150;
param.tp = w1max_TBP_compression(w1_max, TBPmin, param.Q, param.bw);
param.delta_t = param.tp /2;
p1 = LinearChirp(param);

TBPmin = 100;
param.tp = w1max_TBP_compression(w1_max, TBPmin, param.Q, param.bw);
param.delta_t = p1.tp + param.tp /2;
p2 = LinearChirp(param);

TBPmin = 50;
param.tp = w1max_TBP_compression(w1_max, TBPmin, param.Q, param.bw);
param.delta_t = p1.tp + p2.tp + param.tp /2;
p3 = LinearChirp(param);

seq.pulses = {p1 p2 p3};
seq.total_time = p1.tp + p2.tp + p3.tp;
plot_seq(seq)