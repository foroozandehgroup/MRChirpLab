function pulse_modif_test()

par.bw = 300e3;
par.tp = 1000e-6;
par.Q = 5;
par.tres = 0.5e-6;
par.type = "superGaussian";

p0 = LinearChirp(par);

p1 = pulse_modif(p0, "w1", 20e3);
p2 = pulse_modif(p0, "bw", 100e3);
p3 = pulse_modif(p0, "tp", 500e-6);
p4 = pulse_modif(p0, "phi0", pi/2);
p5 = pulse_modif(p0, "delta_t", 1000e-6);
p6 = pulse_modif(p0, "delta_f", 150e3);
p7 = pulse_modif(p0, "n", 16);

plot_pulse({p0, p1, p2, p3, p4, p5, p6, p7});
plot_pulse({p0, p1, p2, p3, p4, p5, p6, p7}, 'polar');

p0 = pulse_modif(p0, "type", "sinsmoothed");

p1 = pulse_modif(p0, "w1", 20e3);
p2 = pulse_modif(p0, "bw", -p0.bw); % reverse sweep
p3 = pulse_modif(p0, "bw", 100e3);
p4 = pulse_modif(p0, "tp", 500e-6);
p5 = pulse_modif(p0, "phi0", pi/2);
p6 = pulse_modif(p0, "sm", 25);

plot_pulse({p0, p1, p2, p3, p4, p5, p6});
plot_pulse({p0, p1, p2, p3, p4, p5, p6}, 'polar');

end