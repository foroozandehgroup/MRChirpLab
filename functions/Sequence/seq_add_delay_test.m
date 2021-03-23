function seq_add_delay_test()

par = struct('t90min', 350e-6, 't180min', 684e-6, 'bw', 200e3, 'tres', 0.5e-6);

seq0 = chorus(par);

seq1 = seq_add_delay(seq0, 350e-6, 1);
seq2 = seq_add_delay(seq0, 350e-6, 2);
seq3 = seq_add_delay(seq0, 350e-6, 3);
seq4 = seq_add_delay(seq0, 350e-6, 4);
seq5 = seq_add_delay(seq0, 350e-6, 5);

disp(seq0.total_time)
disp(seq1.total_time)
disp(seq2.total_time)
disp(seq3.total_time)
disp(seq4.total_time)
disp(seq5.total_time)

plot_seq(seq0)
plot_seq(seq1)
plot_seq(seq2)
plot_seq(seq3)
plot_seq(seq4)
plot_seq(seq5)

end