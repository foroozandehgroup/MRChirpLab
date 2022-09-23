function tp = w1max_TBPmin_MRChirp(pulse_param, w1max, TBPmin)
% Returns minimum acceptable duration tp of a MRChirp for given pulse
% parameters and TBPmin (pulse defined with w1max)
%
% Input:
%   - pulse_param, a structure containing the parameters of the MRChirp
%   - w1max, desired maximum excitation field value (Hz)
%   - TBPmin, desired minimum time bandwidth product  
%
% Output:
%   - tp, pulse compressed duration (s)
%
% similar to w1max_TBP_compression() but more general

if isfield(pulse_param, 'tp')
    error('tp should not be indicated to use w1max and TBPmin compression')
end

pulse_param.w1 = w1max;
p = MRchirp(pulse_param);

% tp for maximum excitation power
tp = p.tp;

% checking minimum TBP factor condition and ajdusting t_p if necessary
if p.bw * tp < TBPmin
    tp = TBPmin / p.bw;
    disp("TBP limited for Q = " + pulse_param.Q)
else
    disp("B1max limited for Q = " + pulse_param.Q)
end
    
end