function hard90_test()

seq = hard90(12e-6,6e-6);

plot_seq(seq)

% offsets
offs = linspace(-300e3, 300e3, 1024);

magn = magn_calc_rot(seq.pulses, seq.total_time, offs);

plot_magn(magn, offs)


spectrum = magn(1,:) + 1i * magn(2,:);

phi0 = 0 * pi/180;
phi1 = 12e-6*2* pi ; % * 300e3/1000

% 2*pi*off*tau_phase_evo

x = (0:length(spectrum)-512); % *1/length(spectrum)

spectrum2 = imag(spectrum .* exp(-1i * (phi0 + offs * phi1)));

figure();

p=plot(offs,spectrum2);
ylim([-1 1])

set(gcf,'Position',[0 0 560 120])

% p.XDataSource = 'offs';
% p.YDataSource = 'spectrum2';
% spectrum2 = imag(spectrum .* exp(-1i * (180 * pi / 180 + x * 0e-6)));refreshdata

end