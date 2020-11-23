function plot_pulse(pulses, plot_type, titles)
% Plots the cartesian component of one pulse or several pulses contained in
% the structure pulse
%
% Input:
%   Required:
%   - pulses, pulse structure or a cell array of pulse to plot
%   Optional:
%    - plot_type, type of plot, it can take the following values
%   "cartesian", "polar", "Xepr" (default: "" for cartesian)
%   - titles, string or list of string to put titles on the
%   plots/subplots of the pulses. If used, requires an input for plot_type

if ~exist('plot_type','var') || plot_type == ""
    plot_type = "cartesian";
end

if exist('titles','var')
    grumble(pulses, plot_type, titles)
else
    grumble(pulses, plot_type)
end

if size(pulses) == 1
    fig = figure();
else
    fig = figure('Position', [0 0 560 840]);
end

for i = 1:length(pulses)
    
    if size(pulses) == 1
        p = pulses;
    else
        p = pulses{i};
    end
    subplot(length(pulses),1,i)
    
    if plot_type == "cartesian"
        
        plot(p.t, p.Cx, 'Color', [0 0.25 0.75])
        hold on
        plot(p.t, p.Cy, 'Color', [0.9 0.05 0.05])
        hold off
        
        if size(pulses) < 2
            ylabel('C_x and C_y (Hz)')
        end
        
    elseif plot_type == "polar"
        
        plot(p.t, p.Pr, 'Color', [0 0.25 0.75])
        
        if size(pulses) < 2
            ylabel('\omega_1(t) (Hz)')
        end
        
        hold on
        yyaxis right
        plot(p.t, p.Pph, 'Color', [0.9 0.05 0.05])
        hold off
        
        if size(pulses) < 2
            ylabel('\phi (rad)')
        end
        
        plt = gca;
        plt.YAxis(2).Color = 'k';
        
    elseif plot_type == "Xepr"
        
        shape = pulse_Xepr_format(p);
        plot(p.t, shape(:,1), 'Color', [0 0.25 0.75])
        hold on
        plot(p.t, shape(:,2), 'Color', [0.9 0.05 0.05])
        hold off
        
    end
    
    xlim([p.t(1) p.t(end)])
    
    if size(pulses) < 2
        xlabel('t (s)')
    end
    
    if exist('titles','var')
        title(titles(i), 'Interpreter', 'none')
    end
end

end


function grumble(pulses, plot_type, titles)

if ~isstruct(pulses) && ~iscell(pulses)
    error("the input pulses must be a cell array or one pulse structure.")
end

if ~contains(plot_type, ["cartesian", "polar", "Xepr"])
    error("plot_type needs to take one of the following value: '' 'cartesian' 'polar' 'Xepr'")
end

if exist('titles','var')
    if length(pulses) ~= length(titles)
        error("the length of pulses and titles must be the same")
    end
end

end