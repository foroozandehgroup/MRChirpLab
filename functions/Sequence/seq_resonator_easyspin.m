function compensated_seq = seq_resonator_easyspin(seq, f, H_f, nu, opt_simulate)
% [EPR] Makes the sequence pulses compensate for resonator effects witht
% the function resonator from EasySpin
% 
% Input:
%   - seq, sequence to be compensated
%   - f, frequencies f of the transfer function (GHz)
%   - H_f, transfer function of the resonator
%   - nu, center frequency of the resonator compensation (GHz)
%   - opt_simulate, to indicate 'compensate' or 'simulate' mode
%
% Output:
%   - compensated_seq, compensated sequence
%
% For more details, see pulse_resonator_easyspin

if nargin > 4
    for i = 1:length(seq.pulses)
        seq.pulses{i} = pulse_resonator_easyspin(seq.pulses{i}, f, H_f, nu, opt_simulate);
    end
else
    for i = 1:length(seq.pulses)
        seq.pulses{i} = pulse_resonator_easyspin(seq.pulses{i}, f, H_f, nu);
    end
end

compensated_seq = seq;

end