function prochorus_test()
% prochorus test function


%% no display and no optimization
par = struct('bw', 300e3, 'tres', 0.5e-6, 't90min',500e-6, ...
             't180min', 1000e-6, 'comression', false, ...
             'pulse_param', struct('n', 20)); % 

seq1 = prochorus(par);

plot_seq(seq1);

%% adding a delay

% par.t_delay = 335e-6;
seq2 = prochorus(par);
plot_seq(seq2)

%% compressed version
% ~1 min for optimization
par.compression = true;
par.display_result = true;
par.phase_ph_opt = true;
par.ph_opt_alpha = 0;
par.ph_opt_parnb = 4;

% uncomment to avoid optimization
% par.ph_opt_x = [0.00680252 1.17786831 -118.60818 -61.900567];

seq3 = prochorus(par);

off = linspace(-seq3.bw/2, seq3.bw/2, 300);

% opt.pc = seq3.pc; % phase cycling
magn = magn_calc_rot(seq3.pulses, seq3.total_time, off); %, opt
plot_magn(magn,off)

%% uncompressed version
% ~2 min for optimization

par.compression = false;
par.display_result = true;
par.phase_sg_opt = true;
par.ph_opt_alpha = 1;
par.ph_opt_parnb = 4;

% uncomment to avoid optimization
% par.ph_opt_x = [1.334597 5.973826 -121.77302 -36.405419];

seq4 = prochorus(par);

end