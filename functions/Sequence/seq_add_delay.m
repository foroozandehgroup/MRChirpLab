function seq = seq_add_delay(seq, td, pos)
% Adds a delay in a sequence
%
% Input:
%   - seq, the sequence
%   - td, the delay to be added
%   - pos, the position of the delay in the duration of the sequence 
%   seq.tau
%
% Output:
%   - seq, the sequence with th added pulses at the end

grumble(pos)

% adding delay duration
if pos == length(seq.tau)
    seq.tau = [seq.tau td];
else
    
    seq.tau = [seq.tau(1:pos-1) td seq.tau(pos:end)];
    
    % pulses shift
    for i = 1:length(seq.pulses)

        if seq.pulses{i}.delta_t > sum(seq.tau(1:pos-1))

            % pulse position
            seq.pulses{i}.delta_t = seq.pulses{i}.delta_t + td;

            % time vector of the pulse
            seq.pulses{i}.t = seq.pulses{i}.t + td;

        end
    end
    
end



seq.total_time = sum(seq.tau);

end

function grumble(pos)

if pos < 2
    error("The position pos should be > 1.")
end

end