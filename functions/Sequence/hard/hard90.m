function seq = hard90(tp, tres)
% Sequence with a single hard 90 degreee pulse
%
% Input:
%   - tp, duration of the pulse
%   - (optional) tres, the time resolution of the pulse 
% Output:
%   - seq, a structure containing the sequence information
%
% Fields contained in seq:
%   - tau, a vector containing the duration of the pulses and delays in
%   order (s)
%   - pulses, a cell array containing the pulse structures (LinearChirp)
%   - total_time, the total time of the pulse sequence (s)

par.tp = tp;
par.alpha = pi/2;

if nargin > 1
    par.tres = tres;
end

seq.pulses{1} = hardpulse(par);
seq.total_time = par.tp;

end