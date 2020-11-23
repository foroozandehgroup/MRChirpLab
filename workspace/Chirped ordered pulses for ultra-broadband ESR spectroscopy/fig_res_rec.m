% plot receiver profile and resonator profiles (x-LiPc at X-band
% and bisnitroxide at Q-band)

clear all;
close all;

load rec_prof_x
load rec_prof_y
x = x - 50;

load Xband_res_prof_f
load Xband_res_prof_H_f
f = (f - 9.45) * 1e3;

f2 = load('Qband_res_prof_f.mat');
H_f2 = load('Qband_res_prof_H_f.mat');
f2 = f2.f;
H_f2 = H_f2.H_f;
f2 = (f2 - 33.945) * 1e3;

figure()

subplot(3,1,1)
plot(x,data,'Color','#345fe4', 'Color','#345fe4','Linewidth',1)
ylim([0 1])
xlim([-350 350])
h = gca;
set(h,'FontSize',10, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-300 -150 0 150 300]);
set(h, 'ytick', [0 0.25 0.5 0.75 1]);
xlabel('$\Omega$ / $2 \pi$ (MHz)','Interpreter','latex','FontSize',10);
ylabel('$I$/$I_{\mathrm{max}}$','Interpreter','latex','FontSize',10);

subplot(3,1,2)
plot(f, H_f, 'Color','#345fe4','Linewidth',1)
xlim([-300 300])
ylim([0 60])
h = gca;
set(h,'FontSize',10, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-300 -150 0 150 300]);
xlabel('Carrier frequency (MHz)','Interpreter','latex','FontSize',10);
ylabel('$\omega_{1}$/$2\pi$ (MHz)','Interpreter','latex','FontSize',10);

subplot(3,1,3)
plot(f2, H_f2, 'Color','#345fe4','Linewidth',1)
xlim([-150 150])
ylim([0 90])
h = gca;
set(h,'FontSize',10, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-150 -100 -50 0 50 100 150]);
xlabel('Carrier frequency (MHz)','Interpreter','latex','FontSize',10);
ylabel('$\omega_{1}$/$2\pi$ (MHz)','Interpreter','latex','FontSize',10);

annotation('textbox',[0.0, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',10)
annotation('textbox',[0.0, 0.66, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',10)
annotation('textbox',[0.0, 0.36, 0, 0],'interpreter','latex','String','$\bf{c)}$','FontSize',10)

set(gcf,'units','centimeters','position',[0,0,12,24])

print(gcf,'-depsc','figure_S4')
