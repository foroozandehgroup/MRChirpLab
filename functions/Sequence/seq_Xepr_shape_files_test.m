function seq_Xepr_shape_files_test()

% sequence
param.bw = 1e9;
param.t90min = 64e-9;
param.t180min = 128e-9;
param.tres = 0.5e-9;
param.pulse_param.n = 10;
param.display_result = true;

seq = exc_3fs(param);

% pulse shape Xepr files
ID = 'CHORUS';
IDnb = '500';
path = [pwd '\'];

seq_Xepr_shape_files(seq, ID, IDnb, path)
close all

end