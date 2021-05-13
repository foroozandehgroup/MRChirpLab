function seq = seq_add_delay(seq, td, pos)
% Adds a delay in a sequence
%
% Input:
%   - seq, the sequence
%   - td, the delay to be added
%   - pos, the position at which the delay is inserted (in seq.tau)
%
% Output:
%   - seq, the sequence with th added pulses at the end


grumble(pos)

if pos == 1 % before the sequence
    
    % adding delay duration
    seq.tau = [td seq.tau];
    
    % pulses shift
    for i = 1:length(seq.pulses)

        % pulse position
        seq.pulses{i}.delta_t = seq.pulses{i}.delta_t + td;

        % time vector of the pulse
        seq.pulses{i}.t = seq.pulses{i}.t + td;
    end
    
elseif pos < length(seq.tau)+1 % in the sequence 
    
    % adding delay duration
    seq.tau = [seq.tau(1:pos-1) td seq.tau(pos:end)];
    
    % pulses shift
    for i = 1:length(seq.pulses)
        
        % select pulses after the added delay
        if seq.pulses{i}.delta_t > sum(seq.tau(1:pos-1)) 
            
            % pulse position
            seq.pulses{i}.delta_t = seq.pulses{i}.delta_t + td;
            
            % time vector of the pulse
            seq.pulses{i}.t = seq.pulses{i}.t + td;
        end
    end
    
elseif pos == length(seq.tau)+1 % after the sequence
    
    % adding delay duration
    seq.tau = [seq.tau td];  
    
end

seq.total_time = sum(seq.tau);

end

function grumble(pos)

if floor(pos)~=pos || pos < 1
    error("The position pos should be a positive integer.")
end

end