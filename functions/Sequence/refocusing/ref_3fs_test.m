function ref_3fs_test()


par = struct('tres', 0.5e-6, 'bw', 400e3, 't180min', 500e-6, ...
             't_delay', 0e-6, 'display_result', true, 'Q180', 5);

seq = ref_3fs(par);

% offsets
n_off = 81;
off = linspace(-seq.bw/2, seq.bw/2, n_off);

opt.magn_init = repmat([0,1,0]', 1, n_off);

opt.pc = seq.pc;

opt.B1 = linspace(0.5, 1.5, 16);

disp('Magnetization computation...')
final_magn = magn_calc_rot(seq.pulses, seq.total_time, off, opt);

plot_magn(final_magn,off,opt.B1)

end