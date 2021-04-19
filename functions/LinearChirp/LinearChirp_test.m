function LinearChirp_test()

par0.bw = 300000;
par0.tp = 500e-6;

% not enough inputs test
try
    LinearChirp(par0);
catch error
    disp(error.message)
end

par0.w1 = 6.491401671720040e+03;

%3 parameters definition
pulse_3_parameters = LinearChirp(par0);

% Q defined
par0.Q = par0.w1^2 * 2 * pi * par0.tp / par0.bw;

par0 = rmfield(par0, 'bw');
pulse_Q_defined_no_bw = LinearChirp(par0);
par0.bw = 300000;

par0 = rmfield(par0, 'w1');
pulse_Q_defined1_no_w1 = LinearChirp(par0);
par0.w1 = 6.491401671720040e+03;

par0 = rmfield(par0, 'tp');
pulse_Q_defined1_no_tp = LinearChirp(par0);
par0.tp = 500e-6;

par0 = rmfield(par0, 'Q');

% should be the same
pulses = {pulse_3_parameters pulse_Q_defined_no_bw ...
          pulse_Q_defined1_no_w1 pulse_Q_defined1_no_tp};
titles = ["3 parameters" "Q defined no bw" ...
          "Q defined no w1"  "Q defined no tp"];

plot_pulse(pulses, "",titles);


% 6 parameters definitions
par0.delta_t = 500e-6;
pulse_delta_t_500us = LinearChirp(par0);
par0 = rmfield(par0, 'delta_t');

par0.phi0 = pi;
pulse_phi0_180deg = LinearChirp(par0);
par0 = rmfield(par0, 'phi0');

par0.delta_f = 5e4;
pulse_delta_f_50kHz = LinearChirp(par0);
par0 = rmfield(par0, 'delta_f');

pulses = {pulse_3_parameters pulse_delta_t_500us ...
          pulse_phi0_180deg pulse_delta_f_50kHz};
titles = ["3 parameters" "delta_t = 500us" ...
          "phi0 = 180deg"  "delta_f = 50kHz"];
plot_pulse(pulses, "", titles);

% optional parameters
par0.n = 2;
pulse_n_2 = LinearChirp(par0);
par0 = rmfield(par0, 'n');

par0.tres = 2e-6;
pulse_tres_2us = LinearChirp(par0);
par0 = rmfield(par0, 'tres');

pulses = {pulse_3_parameters pulse_n_2 pulse_tres_2us};
titles = ["3 parameters" "n = 2" "tres = 2us"];
plot_pulse(pulses, "", titles);

% WURST type test
par1 = struct('bw', 300e3, 'w1', 6491.402, 'tp', 500e-6, ...
              'tres', 0.5e-6, 'type', 'WURST');

wurst = LinearChirp(par1);

par1.n = 5;
wurst_n5 = LinearChirp(par1);

titles = ["WURST80 (default)", "WURST5"];
plot_pulse({wurst, wurst_n5}, "", titles);

% sinsmoothed type test
par2 = struct('bw', 300e3, 'w1', 6491.402, 'tp', 500e-6, ...
              'tres', 0.5e-6, 'type', 'sinsmoothed');
sinsmoothed_default = LinearChirp(par2);

par2.sm = 40;
sinsmoothed40 = LinearChirp(par2);

pulses = {sinsmoothed_default, sinsmoothed40};
titles = ["sinsmoothed default", "sinsmoothed sm=40"];
plot_pulse(pulses, "polar", titles)

% reverse sweep test
par2.sm = 10;
par2.bw = - par2.bw;
reversed_sweep = LinearChirp(par2);

pulses = {sinsmoothed_default, reversed_sweep};
titles = ["sinsmoothed default", "sinsmoothed reverse sweep"];
plot_pulse(pulses, "polar", titles)

% test warning for high number of points
par2.tres = 5e-9;
pulse_low_tres = LinearChirp(par2);

end