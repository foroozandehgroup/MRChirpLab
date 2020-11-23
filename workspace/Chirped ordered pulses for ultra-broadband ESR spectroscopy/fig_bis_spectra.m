% Comparison of bis-nitroxide spectra obtained at Q-band

clear all; clc; close all;

%  --- DATA ---

floc = [pwd '/data/'];

% spectrum from FS and echo

[x0,data0,par0] = eprload(strcat(floc,'014_fidFS'));
data0 = data0 - mean(data0(1:20));
x0 = (x0*2.8e-3) - 33.818; % frequency domain 

% limit: actually recorded at 33.879612GHz (hence manual adjustment...)

% x0 = x0(1:261); % centering
% data0 = data0(1:261);
data0 = flip(data0);
x0 = x0 - (x0(end)+x0(1));

% spectra from chirps

% 2313 - resonator compensation - 12050G
% 33.8/12.050 = 2.8 * 1.005
[x1, data1, par1] = eprload(strcat(floc, '025_bis_2ch_phased.DTA'));

% 3513 family - resonator compensated - 12050G
[x2, data2, par2] = eprload(strcat(floc, '031_chQ_phased.DTA'));

% 3501 family - non resonator compensation - 12050G
[x3, data3, par3] = eprload(strcat(floc, '013_bis_chQ_phased.DTA'));

% Hard pulses - hahn echo - 12064G
% 33.8/12.064 = 2.8 * 1.002
[x4, data4, par4] = eprload(strcat(floc, '001_hard_raw.DTA')); % 

fid = [data4; zeros(2048 - length(data4),1)];
y = fftshift(fft(fid));
% PZphasetool(y); % activate/desactivate
phase_values = [-73;636;1024];
data4 = phase(y, phase_values); % phase correction

fs = 1/0.5e-9;
f = fs * (0:1:length(fid)-1) /length(fid) ; % frequency axis
x4 = f - max(f)/2; % center around zero
x4 = x4 * 1e-9;


disp('Intensity variation from CHORUS to KB:')
% only need to account for video gain (24dB vs 30dB) 
% phase cycling additional scans and AveragesPerScan compensating (1024*16 and 256*64) 
disp("Max:" + num2str(1.995*max(real(data2))/max(real(data1))))
disp("Sum:" + num2str(1.995*sum(abs(real(data2))) / sum(abs(real(data1)))))
               
% resonator vs non resonator
plot2_comp(1.005*x0, data0, 'Field sweep echo', ...
           x3, data3, 'CHORUS - no resonator compensation', ...
           1.005*x0, data0, 'Field sweep echo', ...
           x2, data2, 'CHORUS - resonator compensation');
       
% hardvs2vs3 chirps
plot3_comp(1.005*x0, data0, 'Field sweep echo', x4, data4, 'Hahn echo', ...
           1.005*x0, data0, 'Field sweep echo', x1, data1, 'Kunz scheme', ...
           1.005*x0, data0, 'Field sweep echo', x2, data2, 'CHORUS');

% 2vs3 chirps
plot1_comp(x1,data1, x2,data2);

function plot1_comp(x1, data1,x2,data2)

figure()
plot(1e3*x1, real(data1),'Color','#15b9c3','Linewidth',1);
hold on
% scaled by 1.995 for 6dB difference in video gain
plot(1e3*x2, 1.995*real(data2),'Color','#345fe4','Linewidth',1);
hold off
xlim([-250 250]);
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',10);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',10)
h = gca;
set(h,'FontSize',10, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(gcf,'units','centimeters','position',[0,0,12,8])

print(gcf,'-depsc','figure_S9')

end

function plot2_comp(x1, data1, legend1, x2, data2, legend2, ...
                    x3, data3, legend3, x4, data4, legend4)


% plot
fig = figure();
sbh1 = subplot(1,2,1);
plot(1e3*x1, real(data1)/max(real(data1)),'Color','#15b9c3','Linewidth',1);
hold on
plot(1e3*x2, real(data2)/max(real(data2)),'Color','#345fe4','Linewidth',1);
hold off
xlim([-250 250]);
ylim([0 1]);
% legend({legend1},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0.04, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',8)
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

sbh2 = subplot(1,2,2);
plot(1e3*x3, real(data3)/max(real(data3)),'Color','#15b9c3','Linewidth',1);
hold on
plot(1e3*x4, real(data4)/max(real(data4)),'Color','#345fe4','Linewidth',1);
hold off
xlim([-250 250]);
ylim([0 1]);
% legend({legend2},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0.48, 0.97, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',8)
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

set(gcf,'units','centimeters','position',[0,0,16,4])

print(gcf,'-depsc','figure_S8')

end


function plot3_comp(x1, data1, legend1, x2, data2, legend2, ...
                    x3, data3, legend3, x4, data4, legend4, ...
                    x5, data5, legend5, x6, data6, legend6)


% plot
fig = figure();

sbh1 = subplot(3,1,1);
plot(1e3*x1, real(data1)/max(real(data1)),'Color','#15b9c3','Linewidth',1);
hold on
plot(1e3*x2, real(data2)/max(real(data2)),'Color','#345fe4','Linewidth',1);
hold off
xlim([-250 250]);
ylim([0 1]);
% legend({legend1},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',8)
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);


% scan normalizatio to hard pulse -> 1.28
sbh2 = subplot(3,1,2);
plot(1e3*x3, real(data3)/max(real(data3)),'Color','#15b9c3','Linewidth',1);
hold on
plot(1e3*x4, real(data4)/(max(real(data4))),'Color','#345fe4','Linewidth',1);
hold off
xlim([-250 250]);
ylim([0 1]);
% legend({legend2},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0, 0.97-0.57/2, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',8)
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);


sbh3 = subplot(3,1,3);
plot(1e3*x5, real(data5)/max(real(data5)),'Color','#15b9c3','Linewidth',1);
hold on
plot(1e3*x6, real(data6)/(max(real(data6))),'Color','#345fe4','Linewidth',1);
hold off
xlim([-250 250]);
ylim([0 1]);
% legend({legend2},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0, 0.97-0.57, 0, 0],'Interpreter','latex','String','$\bf{c)}$','FontSize',8)
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);


set(gcf,'units','centimeters','position',[0,0,8,12])

ha=get(gcf,'children');
A = get(ha(1),'position');
B = get(ha(2),'position');
C = get(ha(3),'position');

A(2) = 0.71-0.57;
B(2) = 0.71-0.57/2;
C(2) = 0.71;

set(ha(1),'position',A)
set(ha(2),'position',B)
set(ha(3),'position',C)

print(gcf,'-depsc','figure_5')

end






