function pulse_TopSpin_format_test()

param.bw = 300e3;
param.Q = 5;
param.tp = 1000e-6;

p = LinearChirp(param);
plot_pulse(p,'polar')

% offsets
n_offs = 100;
offs = linspace(-p.bw/2, p.bw/2, n_offs);

magn1 = magn_calc_rot({p},p.tp,0,offs);
plot_magn(magn1,offs)

shp = pulse_TopSpin_format(p);
p.Pr = shp(:,1);
p.Pph = shp(:,2);
plot_pulse(p,'polar',"TopSpin")

p.Pr = p.w1*shp(:,1)/100;
p.Pph = unwrap(wrapTo2Pi((deg2rad(shp(:,2)))));
plot_pulse(p,'polar')

magn2 = magn_calc_rot({p},p.tp,0,offs);
plot_magn(magn2,offs)

end