function chorus_test()

% sequence parameters
ch_param.bw = 300000; % bandwidth
ch_param.t90min = 1e-6* 500;
ch_param.t180min = 1e-6 * 1000;
ch_param.tres = 5e-7;

disp('no display, no polynomial fitting (no magnetization calculation)')
ch_param.phase_polynomial_fitting = false;
ch_param.display_result = false;
tic
chorus_seq1 = chorus(ch_param);
toc
disp('')

disp('display, no polynomial fitting (1 magnetization computation)')
ch_param.phase_polynomial_fitting = false;
ch_param.display_result = true;
tic
chorus_seq2 = chorus(ch_param);
toc
disp('')

disp('no display and polynomial fitting (1 magnetization computation)')
ch_param.phase_polynomial_fitting = true;
ch_param.display_result = false;
tic
chorus_seq3 = chorus(ch_param);
toc

disp('display and polynomial fitting (2 magnetization computation)')
ch_param.phase_polynomial_fitting = true;
ch_param.display_result = true;
tic
chorus_seq4 = chorus(ch_param);
toc

disp('sinsmoothed chorus')
ch_param.pulse_param.type = "sinsmoothed";
ch_param.pulse_param.sm = 5;

ch_param.phase_polynomial_fitting = true;
ch_param.polyfit_degree = 5;
ch_param.display_result = true;
tic
chorus_seq4 = chorus(ch_param);
toc

end













