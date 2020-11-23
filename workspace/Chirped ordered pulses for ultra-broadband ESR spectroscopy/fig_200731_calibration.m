% figure of 2D amplitude optimization
% need to be in the directory of the script

clear all; close all; clc;

floc = [pwd '/data/'];
fname = '015_2DAO.DTA';
[tpulses,magn,par] = eprload(strcat(floc,fname));
magn = abs(magn');

subplot(2,1,1)
colormap default;
pcolor(tpulses{1},tpulses{2},magn);
shading flat;
c=colorbar;
c.Label.String = 'Signal magnitude (a.u.)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 8;
c.TickLabelInterpreter = 'latex';

annotation('textbox',[0, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',8)
xlabel('$90^{\circ}$ pulse amplitude ($\%$)','Interpreter','latex','FontSize',8);
ylabel('$180^{\circ}$ pulse amplitude ($\%$)','Interpreter','latex','FontSize',8);
set(gca,'FontSize',8);set(gca,'TickLabelInterpreter','latex')

set(gcf,'units','centimeters','position',[0,0,8,12])

floc = [pwd '/Data/'];
fname = '012_2DAO.DTA';
[tpulses,magn,par] = eprload(strcat(floc,fname));
magn = abs(magn');

subplot(2,1,2)
colormap default;
pcolor(tpulses{1},tpulses{2},magn);
shading flat;
c=colorbar;
c.Label.String = 'Signal magnitude (a.u.)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 8;
c.TickLabelInterpreter = 'latex';

annotation('textbox',[0, 0.5, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',8)
xlabel('$1^{\mathrm{st}}$ $180^{\circ}$ pulse amplitude ($\%$)','Interpreter','latex','FontSize',8);
ylabel('$2^{\mathrm{nd}}$ $180^{\circ}$ pulse amplitude ($\%$)','Interpreter','latex','FontSize',8);
set(gca,'FontSize',8);set(gca,'TickLabelInterpreter','latex')

set(gcf,'units','centimeters','position',[0,0,8,12])

print(gcf,'-depsc','figure_3')