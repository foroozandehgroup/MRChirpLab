% Will Myers
% 11 Oct., 2020
% CAESR, ICL

% modified by JBV 29/10/2020

% modified by MF 29/10/2020
% using Lorentzian-Gaussian apodization
% using SVD for denoising
% using normalised data and vertical truncation

% CHEESY data processing

clear all; clc; close all;

tf = dir('data/009_echo3513_CHEESY17hrs.DTA');
[b,a,par] = eprload(tf(1).name);

% Lorenz-Gauss apodization
en = 0.15;
gn = 0.05;

for i = 1: size(a,2)

a(1,i)=a(1,i)/2;
x=linspace(0,12,size(a,1))';
a(:,i)=a(:,i).*exp(en*(x)).*exp(-gn*(x.^2));

end

apod = exp(en*(x)).*exp(-gn*(x.^2));

% aab absolute frequency domain data
Af = fftshift(fft(a,2*size(a,1),1),1); %  time domain data FFT
aab = abs(Af');

% f frequencies of time domain fft
bf = b{1};
f = fdaxis(bf(2)-bf(1),size(Af,1)).*1000;
nt = numel(f);

% fE ELDOR frequencies
fE = b{2};

% 
if max(f) > max(fE)
    mxL = max(fE);
else
    mxL = max(f);
end

mxf = max(fE); mnf = min(fE);
aabi = zeros(nt,nt);   % new array for square pts & freqs.

% interpolate ELDOR frequency dimension to number of points of offsets
for ii = 1:nt
    aabi(:,ii) = interp1(fE,aab(:,ii),f(:),'spline');
end

% increasing size to be able to shear the data
aabit = zeros(3*nt,nt);
aabit(nt+1:2*nt,:) = aabi;
aabis = zeros(size(aabit));

% shear operation
for ii = 1:nt
    for jj = 1:nt
        aabis(nt+ii,jj) = aabit(0.5*nt+ii+jj-1,jj);
    end
end

% zoom
aabist = aabis(nt+1:2*nt,:);

maxData = max(max(aabist));

aabist = aabist/norm(aabist);

% aabist(find(aabist>0.1*maxData)) = 0.1*maxData;

% check if we can reduce the noise with SVD

[u,s,v] = svd(aabist);

sCutoff = 10;

aabist = u(:,1:sCutoff)*s(1:sCutoff,1:sCutoff)*transpose(v(:,1:sCutoff));


fig = figure();
pbaspect([1,1,1]);

subplot(3,3,[2 3]);
A = get(gca,'position');
A(1) = A(1) - 0.1;
A(2) = A(2) - 0.08;
set(gca,'position',A);
plot(f,sum(aabist,1)/max(sum(aabist,1)),'Color','#345fe4','Linewidth',1);
xlim([-120 120]);
set(gca,'TickLength',[0.02, 0.01])
set(gca,'FontSize',8, 'TickLabelInterpreter','latex')
set(gca,'YAxisLocation','right')
set(gca,'XTickLabel',[]);
set(gca,'ytick',[]);
%set(gca, 'ytick', [0 1]);
% set(gca, 'Box','off')
set(gca,'Visible','off')

subplot(3,3,[4 7]);
A = get(gca,'position');
A(1) = A(1) - 0.01;
set(gca,'position',A);
plot(sum(aabist,2)/max(sum(aabist,2)),f,'Color','#345fe4','Linewidth',1)
xlim([0 0.2]);
ylim([-100 100]);
set(gca, 'Xdir', 'reverse')
set(gca,'TickLength',[0.02, 0.01])
set(gca,'FontSize',8, 'TickLabelInterpreter','latex')
% set(gca, 'xtick', [0 0.3]);
set(gca,'xtick',[]);
set(gca,'YTickLabel',[]);
set(gca,'Visible','off')

subplot(3,3,[5 6 8 9]);
A = get(gca,'position');
A(1) = A(1) - 0.1;
set(gca,'position',A);
pcolor(f,f,aabist);
% imagesc(f,f,aabist);

xlim([-120 120]);
ylim([-100 100]);

shading flat;
load('mycmap.mat')
% mycmap = colormap(gca); % to save the current colormap
colormap(mycmap);

set(gca,'layer','top')
set(gca,'TickLength',[0.02, 0.01])
set(gca,'FontSize',8, 'TickLabelInterpreter','latex')
set(gca,'YAxisLocation','right')
set(gca, 'xtick', [-100 -50 0 50 100]);
xlabel('$\frac{\Omega}{2 \pi}$ / MHz','Interpreter','latex','FontSize',8);
set(gca, 'ytick', [-100 -50 0 50 100]);
ylabel('$\nu_{HTA}$ / MHz','Interpreter','latex','FontSize',8);

% annotation('textbox',[0.55, 0.87, 0, 0],'interpreter','latex','String','$\bf{*}$','FontSize',8)

set(gcf,'units','centimeters','position',[0,0,8,8])

print(gcf,'-depsc','figure_6')