% 28 Feb.,2019
% CAESR, ICL

% compare Nut vs. pct at different transmitter levels

clear all; clc;

fs = 10; % font size

tfs = ['data/NutVSpct_25tl.dat';
       'data/NutVSpct_50tl.dat';
       'data/NutVSpct_75tl.dat'];

M1 = textread(tfs(1,:));
pct1 = M1(:,1);
f1 = M1(:,2);
P1 = polyfit(pct1,f1,4);
disp(num2str(P1','%1.6e'));
fit1 = P1(1).*pct1.^4+P1(2).*pct1.^3+P1(3).*pct1.^2+P1(4).*pct1+P1(5);

M2 = textread(tfs(2,:));
pct2 = M2(:,1);
f2 = M2(:,2);
P2 = polyfit(pct2,f2,4);
disp(num2str(P2','%1.6e'));
fit2 = P2(1).*pct2.^4+P2(2).*pct2.^3+P2(3).*pct2.^2+P2(4).*pct2+P2(5);

M3 = textread(tfs(3,:));
pct3 = M3(:,1);
f3 = M3(:,2);
P3 = polyfit(pct3,f3,4);
disp(num2str(P3','%1.6e'));
fit3 = P3(1).*pct3.^4+P3(2).*pct3.^3+P3(3).*pct3.^2+P3(4).*pct3+P3(5);

figure(); hold on;

h(1) = plot(pct1,fit1,'Color','#fbda28','Linewidth',1.5);
h(2) = plot(pct2,fit2,'Color','#15b9c3','Linewidth',1.5);
h(3) = plot(pct3,fit3,'Color','#345fe4','Linewidth',1.5);
h(4) = plot(pct1,f1,'o', 'Color','[0 0 0]','Markersize',3);
h(5) = plot(pct2,f2,'square', 'Color','[0 0 0]','Markersize',3);
h(6) = plot(pct3,f3,'^', 'Color','[0 0 0]','Markersize',3);

legend('TM = 25\%', 'TM = 50\%', 'TM = 75\%', 'Interpreter','latex','FontSize',8,'Location','northwest')

xlabel('AWG amplitude input (\%)','Interpreter','latex','FontSize',fs);
ylabel('$\omega_{1}$/$2\pi$ (MHz)','Interpreter','latex','FontSize',fs);

set(gca,'FontSize',fs, 'TickLabelInterpreter','latex')
set(gca, 'xtick', [0 25 50 75 100]);
set(gca, 'ytick', [0 10 20 30 40 50]);

set(gcf,'units','centimeters','position',[0,0,12,8])

print(gcf,'-depsc','figure_S5')