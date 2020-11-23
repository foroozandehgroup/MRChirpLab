function phi_corr = polyfit_ph(pulse, Phase, Options)
% Computes the phase correction to obtain a uniform phase with a
% polynomial fitting
%     
% Input:
%     - pulse, the pulse onto which to apply the phase correction
%     - Phase, vector containing the phases of each magnetization point
%     - Options, field containing two optional parameters:
%           - deg, degree of the polynomial for the fit
%           - disp, option to set to true to display the polyfit results
%
% Output:
%     - phi_corr, array containing the phase corrections in degrees to be
%     applied to each point of the pulse (in the unsmoothed part).


grumble(pulse, Phase, Options);

if isfield(Options, 'polyfit_degree')
    deg = Options.polyfit_degree;
else
    deg = 10;
end

% unsmoothed part beginning index for the magnetization phase vector
i_sm_Phase = floor(length(Phase) * pulse.sm /100);
% beginning and end of smoothing
begin_sm_Phase = i_sm_Phase;
end_sm_Phase = length(Phase) - i_sm_Phase;

% unsmoothed part selection for the magnetization phase vector
xdata = begin_sm_Phase:end_sm_Phase;
ydata = Phase(begin_sm_Phase:end_sm_Phase);

% S and mu: scaling factor used by matlab
% (to avoid warning, needs to be used with polyval too)
[pol, S, mu] = polyfit(xdata, ydata, deg);

% calculate the phase correction for the pulse phase vector
x_phi_corr = linspace(1, length(Phase), pulse.np);
phi_corr = polyval(pol, x_phi_corr, S, mu);

% optional display
if isfield(Options, 'display_result')
    if Options.display_result == true
        
        yfit = polyval(pol, xdata, S, mu);
        figure();
        plot(xdata, ydata, 'k',xdata,yfit,'r')
        ylabel('phase (rad)')
        hold on
        
        yyaxis right
        plot(x_phi_corr, phi_corr, 'b')

        hold off
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'b';
        legend('phase to fit','fit','phase correction')
        ylabel('phase correction (rad)')
        xlabel('point number')

    end
end

end

function grumble(pulse, Phase, Options)

if isfield(Options, 'display_result')
    if ~islogical(Options.display_result)
        error('display_result must be a boolean')
    end
end

end