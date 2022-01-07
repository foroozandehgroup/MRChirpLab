function seq_sum_2TopSpin(seq, ID_seq)

% directory creation at path
total_time = ['_' num2str(seq.total_time * 1e3) 'ms'];
bw = ['_' num2str(seq.bw * 1e-3) 'kHz'];
ID = [ID_seq bw total_time];
folder_name = [pwd '\' ID];

if exist(folder_name, 'dir') == 7
   delete([folder_name '\*'])
   rmdir(folder_name)
   disp("Already existing directory deleted")
end

mkdir(folder_name)
folder_name = [folder_name '\'];

% saving all opened figures (assumed to be relevant to sequence)
fig_list = findobj(allchild(0), 'flat', 'Type', 'figure');
for i = 1:length(fig_list)
  fig_handle = fig_list(length(fig_list)+1-i);
  fig_path   = [folder_name 'fig' num2str(i) '.png'];
  saveas(fig_handle, fig_path);
end
close all

% saving each pulse phase cycling step
pulse_info = [num2str(seq.p.tp * 1e6) 'us' '_' ...
              num2str(max(seq.p.Pr),'%0.6f') 'Hz']; 

for ipc = 1:length(seq.pc(1,:)) % phase cycle step
            
    p = pulse_phase_correction(seq.pulses{1}, seq.pc(1, ipc));
    if length(seq.pulses) > 1
        for i = 2:length(seq.pulses)
            % phase cycling
            p2 = pulse_phase_correction(seq.pulses{i}, seq.pc(i, ipc));
            % sum
            p = pulse_sum(p, p2);
        end
    end
    
    % information for pulse_TopSpin_file
    p.amp = seq.pulses{1}.amp;
    p.phase = seq.pulses{1}.phase;
    p.bw = seq.pulses{1}.bw;

    % TopSpin pulse file
    ID_pulse = [ID_seq '_' num2str(ipc) '_' pulse_info];
    pulse_TopSpin_file(p, ID_pulse, folder_name)

end

end
