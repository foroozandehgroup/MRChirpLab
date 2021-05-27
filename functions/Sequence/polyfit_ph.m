function phi_corr = polyfit_ph(pulse, phase, opt)
% Computes the phase correction to obtain a uniform phase with a
% polynomial fitting
%     
% Input:
%     - pulse, the pulse onto which to apply the phase correction
%     - phase, vector containing the phases of each magnetization point
%     - opt, field containing optional parameters:
%           - polyfit_degree, degree of the polynomial for the fit
%           - display_result, option to set to true to display the polyfit
%             results
%           - start, percentage between 0 and 100 corresponding to the
%           starting point on phase for the polynomial fitting
%           - stop, cf. start but for the stopping point.
%
% Output:
%     - phi_corr, array containing the phase corrections in degrees to be
%     applied to each point of the pulse (in the unsmoothed part).


if nargin < 3
    opt = struct();
end

grumble(pulse, phase, opt);

if isfield(opt, 'polyfit_degree')
    deg = opt.polyfit_degree;
else
    deg = 6;
end

% start and stop point picked at smoothing by default
if isfield(opt, 'start')
    start = opt.start;
else
    start = pulse.sm;
end

if isfield(opt, 'stop')
    stop = opt.stop;
else
    stop = pulse.sm;
end

% start index
i_start = floor(length(phase) * start /100);
% stop index
i_stop = length(phase) - floor(length(phase) * stop /100);

% unsmoothed part selection for the magnetization phase vector
xdata = i_start:i_stop;
ydata = phase(i_start:i_stop);

% S and mu: scaling factor used by matlab
% (to avoid warning, needs to be used with polyval too)
[pol, S, mu] = polyfit(xdata, ydata, deg);

% calculate the phase correction for the pulse phase vector
x_phi_corr = linspace(1, length(phase), pulse.np);
phi_corr = polyval(pol, x_phi_corr, S, mu);

% optional display
if isfield(opt, 'display_result')
    if opt.display_result == true
        
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

function grumble(pulse, phase, opt)

if length(phase) < 20
    error('phase length needs to be superior to 20.') 
end

if isfield(opt, 'display_result')
    if ~islogical(opt.display_result)
        error('display_result must be a boolean')
    end
end

if isfield(opt, 'polyfit_degree')
    if ~isreal(opt.polyfit_degree) || ...
            floor(opt.polyfit_degree) ~= opt.polyfit_degree || ...
            opt.polyfit_degree <= 0
        error('polyfit_degree must be an integer > 0')
    end
end

% checking for unexpected input
input_opt = fieldnames(opt);
for i = 1:length(input_opt)
    
    if ~ismember(input_opt{i}, ...
                 ["polyfit_degree", "display_result", ...
                  "display_result", "start", "stop"])        
        warning(['Careful, ' input_opt{i} ' is not a standard parameter.'])
    end
end

end