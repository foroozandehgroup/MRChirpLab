function tp = w1max_TBP_compression(w1_max, TBPmin, Q, bw)
% Returns minimum acceptable duration tp of a linear chirp for a given
% w1max and TBPmin
%
% Input:
%   - w1max, desired maximum excitation field value (Hz)
%   - TBPmin, desired minimum time bandwidth product  
%   - Q, adiabacitiy factor
%   - bw, bandwidth (Hz)
%
% Output:
%   - tp, pulse duration (s)

    % tp for maximum excitation power
    tp = bw * Q / (w1_max^2 * 2 * pi);
    
    % checking minimum TBP factor condition and ajdusting t_p if necessary
    if bw * tp < TBPmin
        tp = TBPmin / bw;
        disp("TBP limited for Q = " + Q)
    else
        disp("B1max limited for Q = " + Q)
    end
    
end