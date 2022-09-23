function seq = hard_spinecho(tp90, tau_echo)
% Creates a spin echo sequence with hard pulses 90x-180y
%
% Input (properties in the structure param):
%   Required:
%   - tp90, the 90 degree pulse duration (s). The 180 degree pulse duration
%   is defined as twice tp90.
%   - tau_echo, the delay between the 90 degree pulses and the 180 degrees
%   pulses
% Output:
%   - seq, a structure containing the sequence information
%
% Fields contained in seq:
%   - tau, a vector containing the duration of the pulses and delays in
%   order (s)
%   - pulses, a cell array containing the pulse structures (LinearChirp)
%   - total_time, the total time of the pulse sequence (s)

tp180 = 2 * tp90;

% sequence durations
seq.tau = [tp90 tau_echo tp180 tau_echo];

% pulse parameters
par90.tp = tp90;
par90.alpha = pi/2;
par180.tp = tp180;
par180.alpha = pi;
par180.phi0 = pi/2; % phase shift

% 90x
par90.delta_t = tp90/2;
seq.pulses{1} = hardpulse(par90);

% 180y
par180.delta_t = sum(seq.tau(1:2)) + tp180/2;
seq.pulses{2} = hardpulse(par180);

seq.tres = seq.pulses{1}.tres;
seq.total_time = sum(seq.tau);

end