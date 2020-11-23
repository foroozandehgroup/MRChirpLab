% Comparison of excitation profiles - 2 and 3 chirps

clear all; clc; close all;

profile_plot3_comp('006_OneP_FS_FIDs.DTA', [-18;-3.36;8415], 'Hard pulse', ...
                   '016_2Dfs_sum.DTA', [165.6;15.95;4411], 'Kunz', ...
                   '013_2Dfs_sum.DTA', [182.97;14.69;4411], 'CHORUS')

profile_plot4_comp('016_2Dfs_sum.DTA', [165.6;15.95;4411], '2 chirps', ...
                   '018_2Dfs_sum.DTA', [132.6;15.46;4411], 'res. comp.', ...
                   '013_2Dfs_sum.DTA', [182.97;14.69;4411], '3 chirps', ...
                   '007_2Dfs_sum.DTA', [212.4;15.0484;4410], 'res. comp.')

function profile_plot3_comp(fname1, phase_values1, legend1, ...
                            fname2, phase_values2, legend2, ...
                            fname3, phase_values3, legend3)

% retrieving data
floc = [pwd '\Data\'];
floc(strfind(floc,'\'))='/';

[x1,y1,par1] = eprload(strcat(floc,fname1));
fids = y1;
y_fft = fftshift(fft(fids));
y1 = sum(y_fft,2); % summing into 1D data
y1 = phase(y1, phase_values1);

fs = 1/0.5e-9;
f = fs * (0:1:length(fids)-1) /length(fids) ; % frequency axis
x1 = f - max(f)/2; % center around zero
x1 = x1 * 1e-9;

[x2,y2,par2] = eprload(strcat(floc,fname2));
y2 = phase(y2, phase_values2);

[x3,y3,par3] = eprload(strcat(floc,fname3));
y3 = phase(y3, phase_values3);


% plot
fig = figure();

sbh1 = subplot(3,1,1);
plot(1e3*x1, real(y1)/max(real(y1)),'Color','#345fe4','Linewidth',1);
xlim([-250 250]);
ylim([-0.4 1]);
% legend({legend1},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',8)
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

% (2*4*max(real(y1)) - to normalize to hard pulse
sbh2 = subplot(3,1,2);
plot(1e3*x2, real(y2)/max(real(y2)),'Color','#345fe4','Linewidth',1);
xlim([-250 250]);
ylim([-0.4 1]);
% legend({legend2},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0, 0.97-0.57/2, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',8)
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

% (2*4*4*max(real(y1)) - to normalize to hard pulse
sbh2 = subplot(3,1,3);
plot(1e3*x3, real(y3)/max(real(y3)),'Color','#345fe4','Linewidth',1);
xlim([-250 250]);
ylim([-0.4 1]);
% legend({legend2},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0, 0.97-0.57, 0, 0],'interpreter','latex','String','$\bf{c)}$','FontSize',8)
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

disp('Intensity variation from CHORUS to KB:')
% only need to account for phase cycling additional scans (64 vs 16)
% same number of AveragesPerScan (20) same video gain (30 dB)
disp("Max:     " + num2str(max(real(y3)) / (4*max(real(y2)))))


% comparison with individual peaks
pks2 = islocalmax(abs(real(y2))/max(abs(real(y2))),'MinProminence',0.2);
pknb2 = 0;

% getting rid of artifact
for i = 1:length(pks2)
    if pks2(i) == 1
        pknb2 = pknb2 +1;
    end
    if pknb2 == 4
        pks2(i) = 0;
    end
end

x2_pks = x2(pks2);
y2_pks = y2(pks2);

pks3 = islocalmax(abs(real(y3))/max(abs(real(y3))),'MinProminence',0.2);
x3_pks = x3(pks3);
y3_pks = y3(pks3);

% plot of peak picking
figure();
subplot(2,1,1);
plot(x2,abs(real(y2)),x2_pks,abs(real(y2_pks)),'*')
subplot(2,1,2);
plot(x3,abs(real(y3)),x3_pks,abs(real(y3_pks)),'*')

sum2 = 0;
isum2 = 0;
for i = 1:length(x2_pks)
   if abs(x2_pks(i)) <= 0.15
       %disp(x2_pks(i))
       sum2 = sum2 + abs(real(y2_pks(i)));
       isum2 = isum2 +1;
   end
end

sum3 = 0;
isum3 = 0;
for i = 1:length(x3_pks)
   if abs(x2_pks(i)) <= 0.15
       sum3 = sum3 + abs(real(y3_pks(i)));
       isum3 = isum3 +1;
   end
end
disp("MaxPks over bandwidth of interest: " + num2str(sum3 / (4*sum2)) + ...
     " (" + num2str(isum2) + " and " + num2str(isum2) + " peaks picked)")

print(gcf,'-depsc','figure_4')

end

 
function profile_plot4_comp(fname1, phase_values1, legend1, fname2, phase_values2, legend2, ...
                            fname3, phase_values3, legend3, fname4, phase_values4, legend4)

% retrieving data
floc = [pwd '\Data\'];
floc(strfind(floc,'\'))='/';

[x1,y1,par1] = eprload(strcat(floc,fname1));
[x2,y2,par2] = eprload(strcat(floc,fname2));
[x3,y3,par3] = eprload(strcat(floc,fname3));
[x4,y4,par4] = eprload(strcat(floc,fname4));

y1 = phase(y1, phase_values1);
y2 = phase(y2, phase_values2);
y3 = phase(y3, phase_values3);
y4 = phase(y4, phase_values4);

figure();

subplot(2,2,1);
plot(1e3*x1, real(y1)/max(real(y1)),'Color','#345fe4','Linewidth',1);
xlim([-250 250]);
ylim([-0.4 1]);
% legend({legend1},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0.05, 0.97, 0, 0],'interpreter','latex','String','$\bf{a)}$','FontSize',8)
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

subplot(2,2,3);
plot(1e3*x2, real(y2)/max(real(y2)),'Color','#345fe4','Linewidth',1);
xlim([-250 250]);
ylim([-0.4 1]);
% legend({legend2},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0.05, 0.5, 0, 0],'interpreter','latex','String','$\bf{b)}$','FontSize',8)
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

subplot(2,2,2);
plot(1e3*x3, real(y3)/max(real(y3)),'Color','#345fe4','Linewidth',1);
xlim([-250 250]);
ylim([-0.4 1]);
% legend({legend1},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0.48, 0.97, 0, 0],'interpreter','latex','String','$\bf{c)}$','FontSize',8)
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

subplot(2,2,4);
plot(1e3*x4, real(y4)/max(real(y4)),'Color','#345fe4','Linewidth',1);
xlim([-250 250]);
ylim([-0.4 1]);
% legend({legend1},'Location','northwest','Interpreter','latex','FontSize',8);
annotation('textbox',[0.48, 0.5, 0, 0],'interpreter','latex','String','$\bf{d)}$','FontSize',8)
xlabel('$\frac{\Omega}{2 \pi}$ (MHz)','Interpreter','latex','FontSize',8);
ylabel('Intensity (a.u.)','Interpreter','latex','FontSize',8)
h = gca;
set(h,'FontSize',8, 'TickLabelInterpreter','latex')
set(h, 'xtick', [-200 -100 0 100 200]);
set(h, 'ytick', [0 0.5 1]);

set(gcf,'units','centimeters','position',[0,0,16,8])

print(gcf,'-depsc','figure_S7')

end