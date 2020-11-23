function plot_magn(magn, offs)
% Displays the magnetization components of magn at the offsets offs
%
% Input:
%   - magn, magnetization on x, y and z for a certain number of offsets
%   - offs, the offsets vector at which the components of magn are placed
%
% Plot:
%   - Mx, My and Mz, the x,y,z-magnetization components 
%   - Mxy, computed normalized transverse magnetization components
%   - ph, computed magnetization phase


Ix = magn(1,:);
Iy = magn(2,:);
Iz = magn(3,:);
Ixy = sqrt(Ix.^2 + Iy.^2);

ph = magn_phase(magn);

figure('Position', [0 0 560 840]);

subplot(5, 1, 1)
plot(offs, Ix, 'Color',[0 0.25 0.75]);
xlim([offs(1) offs(end)])
ylim([-1 1]);
ylabel('Mx');

subplot(5, 1, 2)
plot(offs, Iy, 'Color',[0 0.25 0.75]);
xlim([offs(1) offs(end)])
ylim([-1 1]);
ylabel('My');

subplot(5, 1, 3)
plot(offs, Iz, 'Color',[0 0.25 0.75]);
xlim([offs(1) offs(end)])
ylim([-1 1]);
ylabel('Mz');

subplot(5, 1, 4)
plot(offs, Ixy, 'Color',[0 0.25 0.75]);
xlim([offs(1) offs(end)])
ylim([-1 1]);
ylabel('Mxy');

subplot(5, 1,5)
plot(offs, ph, 'Color',[0 0.25 0.75]);
xlim([offs(1) offs(end)])
ylim([min(ph) max(ph)]);
ylabel('Phase');

end