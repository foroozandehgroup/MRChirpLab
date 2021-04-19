function pulse_Xepr_shape_file_test()

% pulse parameters
param.bw = 500e6;
param.tp = 128e-9;
param.w1 = 40e3;
param.tres = 0.625e-9;

% pulse
pulse = LinearChirp(param);

% pulse Xepr shape
shape_pulse = pulse_Xepr_format(pulse);

Cx_Xepr = shape_pulse(:,1);
Cy_Xepr = shape_pulse(:,2);
plot(pulse.t, Cx_Xepr, pulse.t, Cy_Xepr)

% pulse shape Xepr file
ID = 'IDnb_nb_phase000';
path = [pwd '\'];
IDnb = '048';

pulse_Xepr_shape_file(pulse, ID, IDnb, path)

end