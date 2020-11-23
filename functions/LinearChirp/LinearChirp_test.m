function LinearChirp_test()

param.bw = 300000;
param.tp = 500e-6;

% not enough inputs test
try
    LinearChirp(param);
catch error
    disp(error.message)
end

param.w1 = 6.491401671720040e+03;

%3 parameters definition
pulse_3_parameters = LinearChirp(param);

% Q defined
param.Q = param.w1^2 * 2 * pi * param.tp / param.bw;

param = rmfield(param, 'bw');
pulse_Q_defined_no_bw = LinearChirp(param);
param.bw = 300000;

param = rmfield(param, 'w1');
pulse_Q_defined1_no_w1 = LinearChirp(param);
param.w1 = 6.491401671720040e+03;

param = rmfield(param, 'tp');
pulse_Q_defined1_no_tp = LinearChirp(param);
param.tp = 500e-6;

param = rmfield(param, 'Q');

% should be the same
pulses = {pulse_3_parameters pulse_Q_defined_no_bw ...
          pulse_Q_defined1_no_w1 pulse_Q_defined1_no_tp};
titles = ["3 parameters" "Q defined no bw" ...
          "Q defined no w1"  "Q defined no tp"];

plot_pulse(pulses, "",titles);


% 6 parameters definitions
param.delta_t = 500e-6;
pulse_delta_t_500us = LinearChirp(param);
param = rmfield(param, 'delta_t');

param.phi0 = pi;
pulse_phi0_180deg = LinearChirp(param);
param = rmfield(param, 'phi0');

param.delta_f = 5e4;
pulse_delta_f_50kHz = LinearChirp(param);
param = rmfield(param, 'delta_f');

pulses = {pulse_3_parameters pulse_delta_t_500us ...
          pulse_phi0_180deg pulse_delta_f_50kHz};
titles = ["3 parameters" "delta_t = 500us" ...
          "phi0 = 180deg"  "delta_f = 50kHz"];
plot_pulse(pulses, "", titles);

% optional parameters
param.n = 2;
pulse_n_2 = LinearChirp(param);
param = rmfield(param, 'n');

param.tres = 2e-6;
pulse_tres_2us = LinearChirp(param);
param = rmfield(param, 'tres');

pulses = {pulse_3_parameters pulse_n_2 pulse_tres_2us};
titles = ["3 parameters" "n = 2" "tres = 2us"];
plot_pulse(pulses, "", titles);

% sinsmoothed type test
param.bw = 300000;
param.tp = 500e-6;
param.w1 = 6.491401671720040e+03;
param.type = "sinsmoothed";
sinsmoothed_default = LinearChirp(param);

param.sm = 40;
param.delta_t = 500e-6;
sinsmoothed40 = LinearChirp(param);

pulses = {sinsmoothed_default, sinsmoothed40};
titles = ["sinsmoothed default", "sinsmoothed sm=40"];
plot_pulse(pulses, "polar", titles)

% test warning for high number of points
param.tres = 5e-9;
pulse_low_tres = LinearChirp(param);

end