function exc_2fs_test()

param.bw = 500000; % bandwidth
param.tres = 5e-7;
param.w1max = 30e3;
param.TBPmin = 150;
param.display_result = true;

seq = exc_2fs(param);

param = rmfield(param, 'TBPmin');
param = rmfield(param, 'w1max');

param.t90min = 886e-6;
seq = exc_2fs(param);

param = rmfield(param, 't90min');
param.t180min = 443e-6;
seq = exc_2fs(param);

end

