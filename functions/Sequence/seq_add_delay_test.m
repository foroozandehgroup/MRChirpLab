function seq_add_delay_test()

par = struct('t90min', 350e-6, 't180min', 684e-6, 'bw', 200e3, 'tres', 0.5e-6);

seq1 = chorus(par);
seq2 = seq_add_delay(seq1, 350e-6, 2);
seq3 = seq_add_delay(seq1, 350e-6, 2);
seq4 = seq_add_delay(seq1, 350e-6, 3);
seq5 = seq_add_delay(seq1, 350e-6, 4);

disp(seq1.total_time)
disp(seq2.total_time)
disp(seq3.total_time)
disp(seq4.total_time)
disp(seq5.total_time)

plot_seq(seq1)
plot_seq(seq2)
plot_seq(seq3)
plot_seq(seq4)
plot_seq(seq5)

end