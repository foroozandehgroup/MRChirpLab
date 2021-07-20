function ref_2fs_test()


par = struct('tres', 0.5e-6, 'bw', 300e3, 't180min', 500e-6, ...
             't_delay', 500e-6, 'display_result', true, 'Q180', 2);

seq = ref_2fs(par);

% offsets
n_off = 101;
off = linspace(-seq.bw/2, seq.bw/2, n_off);

opt.magn_init = repmat([0,1,0]', 1, n_off);
% opt.B1 = linspace(0.5, 1.5, 31);
opt.pc = seq.pc;
opt.pc(2,:) = opt.pc(2,:) - pi;

disp('Magnetization computation...')
final_magn = magn_calc_rot(seq.pulses, seq.total_time, off, opt);
    
plot_magn(final_magn,off) %,opt.B1

end