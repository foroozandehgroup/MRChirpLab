function plot_seq_test()

param.bw = 150000;
param.tp = 300e-6;
param.w1 = 5e+03;
pulse1 = LinearChirp(param);

param.bw = 600000;
param.tp = 600e-6;
param.w1 = 10e+03;
param.n = 10;
param.delta_t = 300e-6 + 600e-6;
pulse2 = LinearChirp(param);

seq.pulses = {pulse1 pulse2};
seq.tres = 1e-7;
seq.total_time = pulse2.t(end);

plot_seq(seq)

ch_param.bw = 300000;
ch_param.t90min = 100e-6;
ch_param.t180min = 200e-6;
ch_param.tres = 5e-7;
ch_param.phase_polynomial_fitting = false;
ch_param.display_result = false;

chorus_seq = chorus(ch_param);

plot_seq(chorus_seq, 'amplitude')
plot_seq(chorus_seq, 'cartesian')

try
    plot_seq(chorus_seq, 'oklm')
catch error
    disp(error.message)
end

end