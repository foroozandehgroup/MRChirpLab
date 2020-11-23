% testing parameters for an inversion pulse

% pulse parameters for DEER inversion in the range of values from article:
% "Optimising broadband pulses for DEER depends on concentration and 
% distance range of interest" (Scherer et Al., 2020)

clear all; close all; clc;

param.tp = 100e-9; % 36
param.bw = 50e6;
param.w1 = 30e6;
% param.Q = 3;
param.type = 'sinsmoothed';
param.sm = 20;
param.tres = 0.625e-9;

inversion_pulse1 = LinearChirp(param);

ph_cy = 0;
n_offs = 200;
offs = linspace(-200e6, 200e6, n_offs);

final_magnetization_1 = magn_calc_rot({inversion_pulse1}, inversion_pulse1.tp, ph_cy, offs);

disp(['TBP = ' num2str(inversion_pulse1.TBP)])
disp(['w1 = ' num2str(1e-6 * inversion_pulse1.w1)])
disp(['Q = ' num2str(inversion_pulse1.Q)])

param.bw = 100e6;
inversion_pulse2 = LinearChirp(param);

final_magnetization_2 = magn_calc_rot({inversion_pulse2}, inversion_pulse2.tp, ph_cy, offs);

disp(['TBP = ' num2str(inversion_pulse2.TBP)])
disp(['w1 = ' num2str(1e-6 * inversion_pulse2.w1)])
disp(['Q = ' num2str(inversion_pulse2.Q)])

param.bw = 200e6;
inversion_pulse3 = LinearChirp(param);

final_magnetization_3 = magn_calc_rot({inversion_pulse3}, inversion_pulse3.tp, ph_cy, offs);

disp(['TBP = ' num2str(inversion_pulse3.TBP)])
disp(['w1 = ' num2str(1e-6 * inversion_pulse3.w1)])
disp(['Q = ' num2str(inversion_pulse3.Q)])

Iz1 = final_magnetization_1(3,:);
Iz2 = final_magnetization_2(3,:);
Iz3 = final_magnetization_3(3,:);

figure();
plot(1e-6*offs, Iz1, 'Color','#fbda28','Linewidth',1)
hold on
plot(1e-6*offs, Iz2, 'Color','#15b9c3','Linewidth',1)
hold on
plot(1e-6*offs, Iz3, 'Color','#345fe4','Linewidth',1)
l1 = ['TBP = ' num2str(inversion_pulse1.TBP)];
l2 = ['TBP = ' num2str(inversion_pulse2.TBP)];
l3 = ['TBP = ' num2str(inversion_pulse3.TBP)];

legend(l1, l2, l3, 'Interpreter','latex','FontSize',8,'Location','northeast')
hold off
ylim([-1 1]);

h = gca;
set(h,'FontSize',10);
set(h,'TickLabelInterpreter','latex')
set(h, 'xtick', [-180 -120 -60 0 60 120 180]);
set(h, 'ytick', [-1 -0.5 0 0.5 1]);

xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',10);
ylabel('$M_{z}$','Interpreter','latex','FontSize',10)
set(gcf,'units','centimeters','position',[0,0,12,6])

print(gcf,'-depsc','figure_S1')