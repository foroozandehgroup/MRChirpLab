function seq_add_delay_test()

par = struct('t90min', 350e-6, 't180min', 684e-6, 'bw', 200e3, 'tres', 0.5e-6);

seq = chorus(par);

plot_seq(seq)
plot_seq(seq_add_delay(seq, 350e-6, 2))
plot_seq(seq_add_delay(seq, 350e-6, 3))
plot_seq(seq_add_delay(seq, 350e-6, 4))

end