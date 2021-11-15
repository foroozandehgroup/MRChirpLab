function chorus_cpmg_test()

% first echo (from chorus)
par = struct('bw', 300e3, 'tres', 0.5e-6, 't90min', 500e-6, ...
             't180min', 1000e-6, 'phase_polynomial_fitting', true, ...
             'tau_echo', 500e-6, 'n', 0, 'display_result', true);

seq = chorus_cpmg(par);

% second echo (from 1st double 180)
par.n = 1;
seq = chorus_cpmg(par);

% third echo (from 2nd double 180)
par.n = 6;
seq = chorus_cpmg(par);

end