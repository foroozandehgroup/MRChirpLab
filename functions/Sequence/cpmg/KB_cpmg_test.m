function KB_cpmg_test()

% first echo (from KB)
par = struct('bw', 300e3, 'tres', 0.5e-6, ...
             't180min', 1000e-6, 'phase_polynomial_fitting', true, ...
             'tau_echo', 500e-6, 'n', 7, 'display_result', true);

seq = KB_cpmg(par);

% second echo (from 1st double 180)
par.n = 2;
seq = KB_cpmg(par);

% third echo (from 2nd double 180)
par.n = 3;
seq = KB_cpmg(par);

end



