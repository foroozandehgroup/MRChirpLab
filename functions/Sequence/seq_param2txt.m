function seq_param2txt(seq, ID, IDnb, path, additional_information)
% Creates a file text saving different properties of the sequence
%
% Input:
%   - seq, the sequence whose parameters are to be saved


filename = [path 'param_' IDnb '_' ID '_param.txt'];

header = ['Sequence ID = ' ID 13 10 ...
          'ID number   = ' IDnb 13 10 ...
          'Total time  = ' num2str(seq.total_time * 1e6) ' us' 13 10 ...
          'Bandwidth   = ' num2str(seq.bw * 1e-6) ' MHz' 13 10 ...
          'Date        = ' num2str(datestr(now, 1)) 13 10 ...
          'Time        = ' num2str(datestr(now, 'HH:MM:SS')) 13 10];
      
dlmwrite(filename, header, 'delimiter', '')

pulses_information = '';
for i = 1:length(seq.pulses)
    
    p = seq.pulses{i};
    pulses_information = [pulses_information ...
                          'Pulse ' num2str(i) 13 10 ...
                          'Type                          = ' char(p.type) 13 10 ...
                          'Excitation field amplitude w1 = ' num2str(p.w1 * 1e-3) ' kHz' 13 10 ...
                          'Duration tp                   = ' num2str(p.tp * 1e6) ' us' 13 10 ...
                          'Adiabaticity factor Q         = ' num2str(p.Q) 13 10 ...
                          'Time resolution tres          = ' num2str(p.tres * 1e6) ' us' 13 10 ...
                          'Number of points np           = ' num2str(p.np) 13 10 ...
                          'Smoothing sm                  = ' num2str(p.sm) ' %' 13 10 ...
                          13 10];
                      
end

dlmwrite(filename, pulses_information, '-append', 'delimiter', '')

footer = '';

if isfield(seq, 'phase_polynomial_fitting') && seq.phase_polynomial_fitting == 1
    footer = [footer ... 
              'Polynomial fit used with degree = ' num2str(seq.polyfit_degree) 13 10];
end

if nargin > 4
    dlmwrite(filename, additional_information, '-append', 'delimiter', '')
end

dlmwrite(filename, footer, '-append', 'delimiter', '')

end