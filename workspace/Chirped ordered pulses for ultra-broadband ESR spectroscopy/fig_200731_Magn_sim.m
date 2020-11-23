% JBV 28.07.2020
% Sequence pulses generation for article on use of linear chirp for 
% broadband excitation in EPR

clear all
close all

tic

% -------------------------------------------------------------------------
% common parameters for both sequences
% -------------------------------------------------------------------------

% sequence required parameters
param.bw = 400e6;      % bandwidth
param.tres = 0.625e-9; % time resolution

% sequence optional parameters
param.display_result = false;
param.Q180 = 4;

% sequence optional pulses parameters
param.pulse_param.type = "sinsmoothed";
param.pulse_param.sm = 12.5;

param.display_result = false;

% offsets
n_offs = 301;
offs = linspace(-param.bw/2, param.bw/2, n_offs);

% -------------------------------------------------------------------------
% Kunz scheme without/with resonator compensation
% -------------------------------------------------------------------------

param.t180min = 160e-9; % 90deg pulse minimum duration

DoubleChirp = doublechirp(param);

% phase cycling
ph1 = [0 0 0 0];
ph2 = [0 1 2 3];
phrec = [0 2 0 2] - [1 1 1 1];

% CTP = [+1 -2]; % coherence transfer pathway
% phrec = phase_cycle_receiver([ph1; ph2], CTP);

ph_cy = pi/2 * [ph1; ph2; phrec];

% magnetization calculation for display
M_DoubleChirp = magn_calc_rot(DoubleChirp.pulses, DoubleChirp.total_time, ph_cy, offs);

My_DoubleChirp = M_DoubleChirp(2,:);

% -------------------------------------------------------------------------
% CHORUS without/with resonator compensation and polynomial fitting
% -------------------------------------------------------------------------

param.t90min = 80e-9; % 90deg pulse minimum duration
param.t180min = 160e-9; % 180deg pulse minimum duration

% with polynomial fitting
param.phase_polynomial_fitting = true;
param.polyfit_degree = 6;

TripleChirp = chorus(param);

% phase cycling
ph1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ph2 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph3 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];

CTP = [1 -2 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph1; ph2; ph3], CTP);

ph_cy = pi/2 * [ph1; ph2; ph3; phrec];

M_TripleChirp = magn_calc_rot(TripleChirp.pulses, TripleChirp.total_time, ph_cy, offs);

plot_magn_articlesup(offs, M_DoubleChirp, M_TripleChirp)
print(gcf,'-depsc','figure_S2')

My_TripleChirp = M_TripleChirp(2,:);
plot_My_articlesup(offs,My_TripleChirp)
print(gcf,'-depsc','figure_S3')
% plot_My_articlesup(offs,sqrt(M_TripleChirp(1,:).^2 + M_TripleChirp(2,:).^2))

toc

function plot_My_articlesup(offs,My_TripleChirp)

figure()

subplot(1,2,1)
plot(1e-6 * offs, My_TripleChirp,'Color','#345fe4','Linewidth',1)

xlabel('$\Omega$ / $2 \pi$ (MHz)','Interpreter','latex','FontSize',10);
ylabel('$M_{y}$','Interpreter','latex','FontSize',10)

h = gca;
set(h,'FontSize',10, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-150 0  150]);
set(h, 'ytick', [0 1]);


subplot(1,2,2)
plot(1e-6 * offs, My_TripleChirp,'Color','#345fe4','Linewidth',1)
xlim([-200 0])
ylim([0.85 1])

xlabel('$\Omega$ / $2 \pi$ (MHz)','Interpreter','latex','FontSize',10);
ylabel('$M_{y}$','Interpreter','latex','FontSize',10)

h = gca;
set(h,'FontSize',10, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -150 -100 -50 0]);
set(h, 'ytick', [0.9 1]);


disp(['My(-150MHz) = ' num2str(My_TripleChirp(26))])
disp(['My(+150MHz) = ' num2str(My_TripleChirp(176))])

e = 100 * (max(My_TripleChirp(51:151)) - min(My_TripleChirp(51:151)));
disp(['Maxium variation between -100MHz and +100MHz = ' num2str(e) '%'])

set(gcf,'units','centimeters','position',[0,0,16,5])

end

function plot_magn_articlesup(offs, magn1, magn2)

Ix = magn1(1,:);
Iy = magn1(2,:);
Iz = magn1(3,:);
Ixy = sqrt(Ix.^2 + Iy.^2);

Phase = magn_phase(magn1);

Ix2 = magn2(1,:);
Iy2 = magn2(2,:);
Iz2 = magn2(3,:);
Ixy2 = sqrt(Ix2.^2 + Iy2.^2);

Phase2 = magn_phase(magn2);

figure(); % -500 -500 1120 1680

subplot(5, 2, 1)
plot(1e-6*offs, Ix,'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{x}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 2)
plot(1e-6*offs, Ix2, 'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{x}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 3)
plot(1e-6*offs, Iy, 'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{y}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 4)
plot(1e-6*offs, Iy2, 'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{y}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 5)
plot(1e-6*offs, Iz, 'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{z}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 6)
plot(1e-6*offs, Iz2, 'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{z}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 7)
plot(1e-6*offs, Ixy, 'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{xy}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 8)
plot(1e-6*offs, Ixy2, 'Color','#345fe4','Linewidth',1);
ylim([-1 1]);
ylabel('$M_{xy}$','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [-1 0 1]);
xticklabels = get(h, 'XTickLabel');
xticklabels{1} = ''; xticklabels{2} = ''; xticklabels{3} = '';
xticklabels{4} = ''; xticklabels{5} = '';
set(h, 'XTickLabel', xticklabels);

subplot(5, 2, 9)
plot(1e-6*offs, Phase, 'Color','#345fe4','Linewidth',1);
ylim([1 7]);
xlabel('$\Omega$ / $2 \pi$ (MHz)','Interpreter','latex','FontSize',8); 
ylabel('$\phi$ (rad)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [2 4 6]);

subplot(5, 2, 10)
plot(1e-6*offs, Phase2, 'Color','#345fe4','Linewidth',1);
ylim([1 7]);
xlabel('$\Omega$ / $2 \pi$ (MHz)','Interpreter','latex','FontSize',8); 
ylabel('$\phi$ (rad)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [2 4 6]);

set(gcf,'units','centimeters','position',[0,0,16,10])

end







