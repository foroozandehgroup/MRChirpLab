function ph_rec = phase_cycle_receiver(ph_cy, ctp)
% Computes the receiver phase cycle for a given coherence transfer pathway
% 
% Input:
%   - ph_cy, list with the phase of each pulse for phase cycling. Dimension
%   of number of phase cycle steps*number of pulses
%   - ctp, coherence transfer pathway of each pulse
% 
% Output:
%   - ph_rec: the receiver phase cycling
%
% ph_pulse[i] is the phase cycle associated to the coherence transfer
% pathway ctp[i].
%
% Example: 90  -  180 - 180 - detection
%         
% ctp = [-1, +2, -2]
% +1                 ******  
% 0 **********      *      * 
% -1          ******        *************
%         
% ph_pulses = [ph1, ph2, ph3]
%         
% ph1 =    0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
% ph2 =    0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
% ph3 =    0 0 1 1 2 2 3 3 0 0 1 1 2 2 3 3
% ph_rec = 0 2 2 0 0 2 2 0 2 0 0 2 2 0 0 2


grumble(ph_cy, ctp)

ph_rec = zeros(1, length(ph_cy));
for i = 1:size(ph_cy, 1)
    ph_rec = ph_rec + ph_cy(i,:) * ctp(i);
end

ph_rec = mod(ph_rec + 16, 4); % modulation

end


function grumble(ph_cy, ctp)

if size(ph_cy,1) ~= length(ctp)
    error(['Coherence transfer pathway list ctp has a length '...
           'not equal to the number of pulses'])
end

end




















