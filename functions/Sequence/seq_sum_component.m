function [t, component_sum] = seq_sum_component(seq, component)
% Returns a the component of a sequence in one vector with its 
% corresponding time vector
%
% Input:
%   - seq, sequence whose component sum is wanted
%   - component, string with the component name (sum is made over 
%   seq.component)
%
% Output:
%   - t, time vector corresponding to the summed component
%   - component_sum, component sum vector

t = [];
component_sum = [];

tres = seq.pulses{1}.tres;

% possible delay at start of the sequence
if seq.pulses{1}.t(1) > 0

    t_delay = tres:tres:seq.pulses{1}.t(1)-tres;
    t = [t t_delay];

    component_sum = [component_sum zeros(1,length(t_delay))];
end

for i = 1:length(seq.pulses)
    
    t = [t seq.pulses{i}.t];
    
    str = "seq.pulses{"+i+"}." + component;
    component_sum = [component_sum eval(str)];
    
    if i+1 <= length(seq.pulses) % check for delay in between pulses
        
        if seq.pulses{i}.t(end) < seq.pulses{i+1}.t(1) - tres
            %seq.pulses{i}.t(1) > t + seq.pulses{i}.tres
            % seq.pulses{i}.t(end) < seq.pulses{i+1}.t(1)
            
            t_delay = seq.pulses{i}.t(end)+tres:tres:seq.pulses{i+1}.t(1) - tres;
            t = [t t_delay];

            component_sum = [component_sum zeros(1,length(t_delay))];
        end
    end
end

% possible delay at the end of the sequence
if seq.pulses{end}.t(end) < seq.total_time

    t_delay = seq.pulses{end}.t(end)+tres:tres:seq.total_time;
    t = [t t_delay];

    component_sum = [component_sum zeros(1,length(t_delay))];
end

end
