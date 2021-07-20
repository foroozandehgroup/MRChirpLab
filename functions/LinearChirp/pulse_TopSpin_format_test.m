function pulse_TopSpin_format_test()

param.bw = 300e3;
param.Q = 5;
param.tp = 1000e-6;
param.tres = 0.5e-6;

Inversion_p = MRchirp(param);
plot_pulse(Inversion_p,'polar')

% offsets
n_offs = 100;
offs = linspace(-Inversion_p.bw/2, Inversion_p.bw/2, n_offs);

magn1 = magn_calc_rot({Inversion_p}, Inversion_p.tp, offs);
plot_magn(magn1,offs)

shp = pulse_TopSpin_format(Inversion_p);
Inversion_p.Pr = shp(:,1);
Inversion_p.Pph = shp(:,2);
plot_pulse(Inversion_p,'polar',"TopSpin")

Inversion_p.Pr = Inversion_p.w1*shp(:,1)/100;
Inversion_p.Pph = unwrap(wrapTo2Pi((deg2rad(shp(:,2)))));
plot_pulse(Inversion_p,'polar')

magn2 = magn_calc_rot({Inversion_p}, Inversion_p.tp, offs);
plot_magn(magn2,offs)

end