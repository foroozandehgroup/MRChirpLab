function doublechirp_test()

param.bw = 500000; % bandwidth
param.tres = 5e-7;
param.w1max = 30e3;
param.Rmin = 150;
param.display_result = true;

seq = doublechirp(param);

param = rmfield(param, 'Rmin');
param = rmfield(param, 'w1max');

param.t90min = 886e-6;
seq = doublechirp(param);

param = rmfield(param, 't90min');
param.t180min = 443e-6;
seq = doublechirp(param);

end

