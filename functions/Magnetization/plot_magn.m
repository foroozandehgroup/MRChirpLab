function plot_magn(magn, off, B1, M_B1_comp)
% Displays the magnetization components of magn at the offsets offs
%
% Input:
%   - magn, magnetization on x, y and z for a certain number of offsets
%   - offs, the offsets vector at which the components of magn are placed
%   - B1, optional argument for magnetization simulated in offset and B1
%   dimensions (only plots My)
%   - M_B1_comp, optional argument to choose which component of the 
%   magnetization to plot
%
% Plot:
%   - Mx, My and Mz, the x,y,z-magnetization components 
%   - Mxy, computed normalized transverse magnetization components
%   - ph, computed magnetization phase

if nargin == 2
    
    Mx = magn(1,:);
    My = magn(2,:);
    Mz = magn(3,:);
    Mxy = sqrt(Mx.^2 + My.^2);

    ph = magn_phase(magn);

    figure('Position', [0 0 560 840]);

    subplot(5, 1, 1)
    plot(off, Mx, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('Mx');

    subplot(5, 1, 2)
    plot(off, My, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('My');

    subplot(5, 1, 3)
    plot(off, Mz, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('Mz');

    subplot(5, 1, 4)
    plot(off, Mxy, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([-1 1]);
    ylabel('Mxy');

    subplot(5, 1,5)
    plot(off, ph, 'Color',[0 0.25 0.75]);
    xlim([off(1) off(end)])
    ylim([min(ph) max(ph)]);
    ylabel('Phase (rad)');

elseif nargin > 2
    
    if nargin == 3 % default
        M_B1_comp = "My";
        M_B1 = magn(2, :, :); 
    elseif M_B1_comp == "Mx"
        M_B1 = magn(1, :, :);
    elseif M_B1_comp == "My"
        M_B1 = magn(2, :, :);
    elseif M_B1_comp == "Mz"
        M_B1 = magn(3, :, :);
    elseif M_B1_comp == "Mxy"
        M_B1 = abs(magn(1,:,:)+1i*magn(2,:,:));
    else
        error(['M_B1 should be a string and can only take the ' ...
               'following values: Mx, My, Mz, Mxy'])
    end
        
    r = reshape(M_B1, length(off), length(B1));
        
    figure('Position', [0 0 1000 400])
    
    subplot(1, 2, 1)
    mesh(off, B1, r')
    xlabel('\Delta\omega_0');
    ylabel("\omega_1 / \omega_1^0");
    zlabel(M_B1_comp);
    
    subplot(1, 2, 2)
    contour(off, B1, r.', [-0.99 -0.95 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99])
    xlabel('\Delta\omega_0');
    ylabel("\omega_1 / \omega_1^0");

else
    error('Wrong number of input arguments')
end

end