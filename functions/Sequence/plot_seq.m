function plot_seq(seq, plot_type)
% Plots the components of a sequence in one figure
%
% Input:
%   - seq, sequence to plot
%   - plot_type, type of plot. Can be 'cartesian' for the sequence 
%   cartesian coordinates (by default) or 'amplitude' to plot the enveloppe
%   of the sequence


if nargin < 2
    plot_type = "cartesian";
end

grumble(seq, plot_type)

figure('Position', [0 0 1000 250])

if plot_type == "cartesian"
    [~,Cx] = seq_sum_component(seq, 'Cx');
    [t,Cy] = seq_sum_component(seq, 'Cy');
    plot(t, Cx, 'Color', '#002147', 'Linewidth', 1)
    
    hold on
    plot(t, Cy, 'Color', '#FC3C3C', 'Linewidth', 1)
    hold off
    ylabel('C_x and C_y (Hz)')
    
elseif plot_type == "amplitude"
    [t,Pr] = seq_sum_component(seq, 'Pr');
    plot(t, Pr, 'Color', '#FC3C3C', 'Linewidth', 1)
    ylabel('\omega_1(t) (Hz)')

end

xlabel('t (s)')
xlim([0 max(t)])

end


function grumble(seq, plot_type)

if ~isfield(seq, 'pulses')
    error("The structure seq must contain a structure pulses.")
end

if ~iscell(seq.pulses)
    error("the sequence pulses seq.pulses must be a cell array.")
end

if ~contains(plot_type, ["cartesian", "amplitude"])
    error("plot_type needs to take one of the following value: 'cartesian', 'amplitude'")
end

end