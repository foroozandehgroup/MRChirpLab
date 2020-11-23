function seq_pulses_disp(seq)
% Displays basic information about the sequence pulses
%
% Input:
%   - seq, sequence whose information are to be displayed
    

data_disp = [];
for p = 1:length(seq.pulses)
    data_disp = [data_disp; 1e3*seq.pulses{p}.tp seq.pulses{p}.TBP 1e-3*seq.pulses{p}.w1];
end

format short g
disp(" ")
header = {'tp (ms)','TBP','w1 (Hz)'};
disp([header; num2cell(data_disp)])
disp("Total time (s): " + seq.total_time)
disp("")
    
end