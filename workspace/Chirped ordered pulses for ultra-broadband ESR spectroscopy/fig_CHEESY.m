% CHEESY data processing

% MF
% based on a code by WM and JBV

% using phase sensitive data

% 0th and 1st order phase correction (optional)

% filtering detection mirror signal and zero-frequency glitch (optional)
% removing a band anti-diagonal for mirror signal
% removing a band column for zero-frequency glitch

% using SVD for denoising (optional)

% possibility of writing data, eitehr raw or processed data in NMR Topspin

% using Topspin as an interface for processing: FT, phase correction, tilting/shearing, etc.
%% importing data

clear all; clc; close all;

tf = dir('data/009_echo3513_CHEESY17hrs.DTA');
[freq_cell,data,par] = eprload(tf(1).name);

% boolean to select the different post-processing options
phase_correction = false;
antidiag_remove = false;
zerofreq_remove = false;
SVD_denoising = true;


%% phase correction

sw = 500e6;

zf = 1024; % zero filling

frq = linspace(-sw/2,sw/2,zf);
frq = frq(:);

tau = 9500/(360*sw); 
% 9340 degrees 1st order correction
% corresponding to 53 ns timing error (tau)

phi0 = -40; % 0th order correction

Data_col_FT = fftshift(fft(data,zf),1); %  time domain data FFTfdata = sum(Data_col_FT_Interp,2);

Data_col_FT_sum = sum(Data_col_FT,2);

if phase_correction == true
    Data_col_FT_sum_phi = Data_col_FT_sum.*exp(1i*2*pi*frq*tau).*exp(1i*2*pi*phi0/360);
else
    Data_col_FT_sum_phi = abs(Data_col_FT_sum);
end
plot(frq,real(Data_col_FT_sum),'LineWidth',1);hold on; ...
    plot(frq,real(Data_col_FT_sum_phi),'LineWidth',1); hold off;

%% apply phase correction to 2D data and interpolating HTA dimension

[~,r] =size(Data_col_FT);


bf = freq_cell{1};
f = fdaxis(bf(2)-bf(1),size(Data_col_FT,1)).*1000;

Data_col_FT_Interp = zeros(zf,zf);

for k = 1:r
    if phase_correction == true
        Data_col_FT(:,k) = ...
            Data_col_FT(:,k).*exp(1i*2*pi*frq*tau).*exp(1i*2*pi*phi0/360);
    else
        Data_col_FT(:,k) = abs(Data_col_FT(:,k));
    end

end

for ii = 1:zf
    Data_col_FT_Interp(ii,:) = interpft(Data_col_FT(ii,:),zf);
end

imagesc(real(Data_col_FT_Interp))

%% removing anti-diagonal mirror signal
if antidiag_remove == true
    
    

    w = 15; %width of the band
    hw = floor(w/2);

    e = ones(zf,w);
    A = spdiags(e,-hw:hw,zf,zf);
    A = full(A);
    A = flipud(A);
    B = ones(zf,zf);
    C = B - A;
    Data_col_FT_Interp_mirr = C.*Data_col_FT_Interp; % mirror image filtered
    imagesc(real(Data_col_FT_Interp_mirr));
else
    Data_col_FT_Interp_mirr = Data_col_FT_Interp;

end

%% removing zero-frequency glitch

if zerofreq_remove == true
    w = 15;
    e = zeros(zf,w);
    B = ones(zf,zf);

    a1 = floor(zf/2)-floor(w/2);
    a2 = floor(zf/2)+floor(w/2);

    B(:,a1:a2)=e;
    B = transpose(B);

    Data_col_FT_Interp_mirr_DC = B.*Data_col_FT_Interp_mirr;
    imagesc(real(Data_col_FT_Interp_mirr_DC));
else
    Data_col_FT_Interp_mirr_DC = Data_col_FT_Interp_mirr;
end
%% denoising with SVD

if SVD_denoising == true
    [u,s,v] = svd(real(Data_col_FT_Interp_mirr_DC));

    sCutoff = 200;

    Data_col_FT_Interp_mirr_DC = u(:,1:sCutoff)*s(1:sCutoff,1:sCutoff)*transpose(v(:,1:sCutoff));
else
    Data_col_FT_Interp_mirr_DC = real(Data_col_FT_Interp_mirr_DC);
end


imagesc(Data_col_FT_Interp_mirr_DC);
%%
Data_col_FT_Interp_mirr_DC = ...
    transpose(abs(Data_col_FT_Interp_mirr_DC));

% increasing size to be able to shear the data
Data_col_FT_Interp_mirr_DC_extended = zeros(3*zf,zf);
Data_col_FT_Interp_mirr_DC_extended(zf+1:2*zf,:) = ...
    Data_col_FT_Interp_mirr_DC;
Data_col_FT_Interp_mirr_DC_extended_sheared = ...
    zeros(size(Data_col_FT_Interp_mirr_DC_extended));

% shearing data
for ii = 1:zf
    for jj = 1:zf
        Data_col_FT_Interp_mirr_DC_extended_sheared(zf+ii,jj) = ...
            Data_col_FT_Interp_mirr_DC_extended(0.5*zf+ii+jj-1,jj);
    end
end
%% plotting CHEESY with projections

Data_col_FT_Interp_mirr_DC_extended_sheared = ...
    Data_col_FT_Interp_mirr_DC_extended_sheared(zf+1:2*zf,:);
imagesc(Data_col_FT_Interp_mirr_DC_extended_sheared)
%%
fig = figure();
pbaspect([1,1,1]);

xProj = sum(Data_col_FT_Interp_mirr_DC_extended_sheared,1);
yProj = sum(Data_col_FT_Interp_mirr_DC_extended_sheared,2);

subplot(3,3,[2 3]);
A = get(gca,'position');
A(1) = A(1) - 0.1;
A(2) = A(2) - 0.07;
set(gca,'position',A);
plot(f,xProj/max(xProj),'Color','#345fe4','Linewidth',1);
xlim([-150 150]);
set(gca,'Visible','off')

subplot(3,3,[4 7]);
A = get(gca,'position');
A(1) = A(1) + 0.01;
set(gca,'position',A);
plot(yProj/max(yProj),f,'Color','#345fe4','Linewidth',1)
xlim([-0.05 0.2]);
ylim([-100 100]);
set(gca, 'Xdir', 'reverse')
set(gca,'Visible','off')

subplot(3,3,[5 6 8 9]);
A = get(gca,'position');
A(1) = A(1) - 0.1;
set(gca,'position',A);
pcolor(f,f,Data_col_FT_Interp_mirr_DC_extended_sheared);
xlim([-150 150]);
ylim([-100 100]);

shading flat;
load('mycmap.mat')
colormap(mycmap);

set(gca,'layer','top')
set(gca,'TickLength',[0.02, 0.01])
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
set(gca,'YAxisLocation','right')
set(gca, 'xtick', [-100 -50 0 50 100]);
xlabel('$\frac{\Omega}{2 \pi}$ / MHz','Interpreter','latex','FontSize',16);
set(gca, 'ytick', [-100 -50 0 50 100]);
ylabel('$\nu - \nu_{HTA}$ / MHz','Interpreter','latex','FontSize',16);

set(gcf,'units','centimeters','position',[0,0,16,16])

% %% using Topspin as an interface for processing: FT, phase correction, tilting/shearing, etc.
% %% ifft to get filtered time-domain back and write to Topspin as raw data
% 
% % Data_filt = ifft(ifftshift(Data_col_FT_Interp_mirr_DC,1)); %  time domain data FFT
% Data_filt = ifft(ifftshift(Data_col_FT_Interp,1)); %  time domain data FFT
% 
% data=Data_filt; % fid is scaled up with a factor (optional), to avoid very low contours in Topspsin
% 
% infile = '/Users/mforoozandeh/esr/data/esr_data/1/ser';
% 
% [r,c]=size(data);
% redata=zeros(r,c);
% imdata=zeros(r,c);
% inc_fid=cell(1,c);
% 
% for n=1:c
%     
%     redata(:,n)=real(data(:,n));
%     imdata(:,n)=imag(data(:,n));
%     
%     inc_fid{n} = horzcat(redata(:,n),imdata(:,n));
%     
%     inc_fid{n} = (inc_fid{n}).';
%     
%     inc_fid{n} = inc_fid{n}(:);
%     
% end
% 
% ser_mat = cell2mat(inc_fid);
% 
% ser_mat = ser_mat(:);
% 
% fid = fopen(infile,'w','l');
% fwrite(fid, ser_mat, 'int32');
% fclose(fid);
% %% alternatively the frequency domain data can be written to Topspin as processed data
% 
% data=Data_col_FT_Interp_mirr_DC; % fid is scaled up with a factor (optional), to avoid very low contours in Topspsin
% 
% infileReal = '/Users/mforoozandeh/esr/data/esr_data/1/pdata/1/2rr';
% infileImage = '/Users/mforoozandeh/esr/data/esr_data/1/pdata/1/2ii';
% 
% data_2rr = real(data);
% data_2ii = imag(data);
% 
% data_2rr_vec = data_2rr(:);
% data_2ii_vec = data_2ii(:);
% 
% fid = fopen(infileReal,'w','l');
% fwrite(fid, data_2rr_vec, 'int32');
% fclose(fid);
% 
% fid = fopen(infileImage,'w','l');
% fwrite(fid, data_2ii_vec, 'int32');
% fclose(fid);