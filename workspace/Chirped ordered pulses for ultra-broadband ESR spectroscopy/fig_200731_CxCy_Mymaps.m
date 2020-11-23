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

param.display_results = true;

% resonator profile information
resonator.f = load('resonator_profile_f.mat').f;
resonator.H_f = load('resonator_profile_H_f.mat').H_f;
resonator.nu = 9.45;

% Kunz scheme
param.t90min = 320e-9; % 90deg pulse minimum duration

DoubleChirp = doublechirp(param);

% CHORUS
param.t90min = 80e-9; % 90deg pulse minimum duration
param.t180min = 160e-9; % 180deg pulse minimum duration
param.phase_polynomial_fitting = true;
param.polyfit_degree = 6;
% param.t_delay = 90e-9;

TripleChirp = chorus(param);

toc

% -------------------------------------------------------------------------
% plot cartesian components
% -------------------------------------------------------------------------

DoubleChirp.total_time = 700e-9;
TripleChirp.total_time = 700e-9;

% figure 1
subplot(4,1,1)
figure_CxCy(DoubleChirp)
A = get(gca,'position');          
A(4) = 0.75*A(4);
A(2) = A(2)+0.07;
set(gca,'position',A);
subplot(4,1,2)
figure_My_map(DoubleChirp,[0 0 0 0]',1)
A = get(gca,'position');          
A(4) = 1.5*A(4);
A(2) = A(2)+0.03;
set(gca,'position',A);
subplot(4,1,3)
figure_CxCy(TripleChirp)
A = get(gca,'position');          
A(4) = 0.75*A(4);
A(2) = A(2)+0.04;
set(gca,'position',A);
subplot(4,1,4)
figure_My_map(TripleChirp,[0 0 0 0]',2)
A = get(gca,'position');          
A(4) = 1.5*A(4);
A(2) = A(2);
set(gca,'position',A);


set(gcf,'units','centimeters','position',[0,0,8,12])

annotation('textbox',[0, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',8)
annotation('textbox',[0, 0.5, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',8)
print(gcf,'-depsc','figure_1')

% figure 1 - with phase cycle

figure()
subplot(4,1,1)
figure_CxCy(DoubleChirp)
A = get(gca,'position');          
A(4) = 0.75*A(4);
A(2) = A(2)+0.07;
set(gca,'position',A);
subplot(4,1,2)
ph1 = [0 0 0 0];
ph2 = [0 1 2 3];
CTP = [-1 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph1; ph2], CTP);
ph_cy = pi/2 * [ph1; ph2; phrec]; % phase cycling
figure_My_map(DoubleChirp,ph_cy,1)
A = get(gca,'position');          
A(4) = 1.5*A(4);
A(2) = A(2)+0.03;
set(gca,'position',A);
subplot(4,1,3)
figure_CxCy(TripleChirp)
A = get(gca,'position');          
A(4) = 0.75*A(4);
A(2) = A(2)+0.04;
set(gca,'position',A);
subplot(4,1,4)
ph1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ph2 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph3 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];
CTP = [1 -2 +2]; % coherence transfer pathway
phrec = phase_cycle_receiver([ph1; ph2; ph3], CTP);
ph_cy = pi/2 * [ph1; ph2; ph3; phrec]; % phase cycling
figure_My_map(TripleChirp,ph_cy,2)
A = get(gca,'position');          
A(4) = 1.5*A(4);
A(2) = A(2);
set(gca,'position',A);


set(gcf,'units','centimeters','position',[0,0,8,12])

annotation('textbox',[0, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',8)
annotation('textbox',[0, 0.5, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',8)
set(gcf,'units','centimeters','position',[0,0,8,6])

% figure 2

OneChirp.bw = DoubleChirp.bw;
OneChirp.tres = DoubleChirp.tres;
OneChirp.total_time = 400e-9;
OneChirp.pulses{1} = DoubleChirp.pulses{1};

figure()
subplot(2,1,1)
figure_CxCy(OneChirp)
A = get(gca,'position');          
A(4) = 0.75*A(4);
A(2) = A(2)+0.07;
set(gca,'position',A);
subplot(2,1,2)
figure_My_map(OneChirp,[0],1)

A = get(gca,'position');          
A(4) = 1.4*A(4);
A(2) = A(2) + 0.03;
set(gca,'position',A);

set(gcf,'units','centimeters','position',[0,0,8,6])
print(gcf,'-depsc','figure_2')


% figure 3

% resonator profile information
f = load('Xband_res_prof_f.mat').f;
H_f = load('Xband_res_prof_H_f.mat').H_f;
nu = 9.45;

path = [pwd '\'];

DoubleChirp = seq_resonator_easyspin(DoubleChirp, f, H_f, nu);
TripleChirp = seq_resonator_easyspin(TripleChirp, f, H_f, nu);

DoubleChirp.total_time = 640e-9;
TripleChirp.total_time = 640e-9;

figure()
subplot(2,1,1)
article_figure_sup(DoubleChirp)
subplot(2,1,2)
article_figure_sup(TripleChirp)
set(gcf,'units','centimeters','position',[0,0,12,8])
print(gcf,'-depsc','figure_S6')

function figure_My_map(seq,ph_cy,n)

% offsets and time axis
n_offs = 301;
offs = linspace(-seq.bw/2, seq.bw/2, n_offs);
t = 0:seq.tres:seq.total_time;

final_magnetization = magn_calc_rot_full(seq.pulses, t, ph_cy, offs);

y_traj = transpose(squeeze(final_magnetization(n,:,:)));

%  figure('Position', [0 0 1000 400])
imagesc(1e9*t, 1e-6*offs, y_traj)
xlh = xlabel('$t$ (ns)','Interpreter','latex','FontSize',8);
ylabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
h = gca;
set(h,'FontSize',8);set(h,'TickLabelInterpreter','latex')
set(h, 'xtick', [0 200 400 600 800]);
set(h, 'ytick', [-200 0 200]);

end


function figure_CxCy(seq)

%  figure('Position', [0 0 1000 400])

[~,Cx] = seq_sum_component(seq, 'Cx');
[t,Cy] = seq_sum_component(seq, 'Cy');

Cx = Cx * 1e-6;
Cy = Cy * 1e-6;
t = t * 1e9;

plot(t, Cx, 'Color','#345fe4','Linewidth',1) % #090088 #4077be #4126b4
hold on
plot(t, Cy, 'Color','#fbda28','Linewidth',1) % #FFBD39 #f8fb13
% legend('$C_{x}$','$C_{y}$','Interpreter','latex','FontSize',8)
hold off

xlim([0 max(t)])
ylim([-45 45])
h = gca;
set(h,'FontSize',8)
set(h,'XTick',[])
set(h,'YTick',[])

end


function article_figure_sup(seq)

% figure('Position', [0 0 1300 300])

[~,Cx] = seq_sum_component(seq, 'Cx');
[t,Cy] = seq_sum_component(seq, 'Cy');

Cx = Cx * 1e-6;
Cy = Cy * 1e-6;
t = t * 1e9;

plot(t, Cx, 'Color','#345fe4','Linewidth',1) % #090088 #4077be #4126b4
hold on
plot(t, Cy, 'Color','#fbda28','Linewidth',1) % #FFBD39 #f8fb13
legend('$C_{x}$','$C_{y}$','Interpreter','latex','FontSize',10)
hold off

xlabel('$t$ (ns)','Interpreter','latex','FontSize',10);

xlim([0 max(t)])
ylim([-45 45])
h = gca;
set(h,'FontSize',10, 'TickLabelInterpreter','latex')
set(h,'XTick',[0 200 400 600])
set(h,'YTick',[])

end


