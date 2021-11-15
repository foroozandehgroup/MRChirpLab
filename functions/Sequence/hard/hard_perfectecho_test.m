function hard_perfectecho_test()

seq = hard_perfectecho(12e-6, 24e-6);

plot_seq(seq, "amplitude")

% offsets
offs = linspace(-300e3, 300e3, 151);

magn = magn_calc_rot(seq.pulses, seq.total_time, [0 0 0 0 0]', offs);

plot_magn(magn, offs)


end