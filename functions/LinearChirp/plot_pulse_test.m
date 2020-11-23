function plot_pulse_test()
% test of the function plot_pulse_test

param.bw = 150000;
param.tp = 300e-6;
param.w1 = 5e+03;
pulse1 = LinearChirp(param);

plot_pulse(pulse1)
plot_pulse(pulse1, "")
plot_pulse(pulse1, "polar")

plot_pulse(pulse1, "", "empty plot type")
plot_pulse(pulse1, "polar", "polar plot type")
plot_pulse(pulse1, "Xepr", "pulse1 - Xepr")

param.bw = 600000;
param.tp = 600e-6;
param.w1 = 10e+03;
param.n = 4;
pulse2 = LinearChirp(param);

pulses = {pulse1 pulse2};

plot_pulse(pulses, "", ["pulse1" "pulse2"])
plot_pulse(pulses, "polar", ["pulse1" "pulse2"])


end