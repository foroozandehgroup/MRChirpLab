function plot_magn(magn, off, B1)
% Displays the magnetization components of magn at the offsets offs
%
% Input:
%   - magn, magnetization on x, y and z for a certain number of offsets
%   - offs, the offsets vector at which the components of magn are placed
%   - B1, optional argument for magnetization simulated in offset and B1
%   dimensions (only plots My)
%
% Plot:
%   - Mx, My and Mz, the x,y,z-magnetization components 
%   - Mxy, computed normalized transverse magnetization components
%   - ph, computed magnetization phase


if nargin == 2
    
    Ix = magn(1,:);
    Iy = magn(2,:);
    Iz = magn(3,:);
    Ixy = sqrt(Ix.^2 + Iy.^2);

    ph = magn_phase(magn);

    figure('Position', [0 0 560 840]);

    subplot(5, 1, 1)
    plot(off, Ix, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('Mx');

    subplot(5, 1, 2)
    plot(off, Iy, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('My');

    subplot(5, 1, 3)
    plot(off, Iz, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('Mz');

    subplot(5, 1, 4)
    plot(off, Ixy, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('Mxy');

    subplot(5, 1,5)
    plot(off, ph, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([min(ph) max(ph)]);
    ylabel('Phase');

elseif nargin == 3
        
    Iy_B1 = magn(2, :, :);
    ry = reshape(Iy_B1, length(off), length(B1));
        
    figure('Position', [0 0 1000 400])
    
    subplot(1, 2, 1)
    mesh(off, B1, ry')
    xlabel('\Delta\omega_0');
    ylabel("\omega_1 / \omega_1^0");
    zlabel('My');
    
    subplot(1, 2, 2)
    contour(off, B1, ry.', [-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99])
    xlabel('\Delta\omega_0');
    ylabel("\omega_1 / \omega_1^0");


end