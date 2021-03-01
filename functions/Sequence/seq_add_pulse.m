function seq = seq_add_pulse(seq, p)
% Adds a pulse at the end of a sequence
%
% Input:
%   - seq, the sequence
%   - p, the pulse to be added at the end of the sequence
%
% Output:
%   - seq, the sequence with th added pulses at the end


% adding pulse duration
seq.tau = [seq.tau p.tp];

% adding pulse to pulses
seq.pulses{end+1} = p;

% position of the pulse
seq.pulses{end}.delta_t = sum(seq.tau(1:end)) - p.tp/2;

% time vector of the pulse
pulse_start = seq.pulses{end}.delta_t - seq.pulses{end}.tp / 2;
seq.pulses{end}.t = pulse_start + linspace(0, p.tres*p.np, p.np);

end