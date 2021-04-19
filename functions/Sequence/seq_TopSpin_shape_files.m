function seq_TopSpin_shape_files(seq, ID, path)
% Saves a sequence pulses as TopSPin shape files in a directory
%
% Input:
%   - seq, sequence to be saved
%   - ID, the identification of the sequence
%   - path, location at which to save the files
%
% Files generated in created directory:
%   - pulses shapes files
%   - text file with the sequence main parameters
%   - figures opened when the function is called (are assumed to be linked
%   to the sequence). Plots of each pulse are previoulsy generated.

if isfield(seq, 'bw') && isfield(seq, 'total_time')
    total_time = ['_' num2str(seq.total_time * 1e3) 'ms'];
    bw = ['_' num2str(seq.bw * 1e-3) 'kHz'];
    ID = [ID bw total_time];
end

% directory creation at path
folder_name = [path ID];

if exist(folder_name, 'dir') == 7
   delete([folder_name '\*'])
   rmdir(folder_name)
   disp("Already existing directory deleted")
end

mkdir(folder_name)
folder_name = [folder_name '\'];

% shape file of each pulse

for i = 1:length(seq.pulses)        

    ID_pulse = [ID '_' num2str(i) '_' ...
                num2str(seq.pulses{i}.tp * 1e6) 'us' '_' ...
                num2str(seq.pulses{i}.w1,'%0.6f') 'Hz']

    pulse_TopSpin_file(seq.pulses{i}, ID_pulse, folder_name)
    % plot_pulse(pulse, "cartesian", convertCharsToStrings(ID_pulse))

end

% all the opened figures assumed to be linked to the sequence
disp('Saving all open figures...')

fig_list = findobj(allchild(0), 'flat', 'Type', 'figure');
for i = 1:length(fig_list)
  fig_handle = fig_list(length(fig_list)+1-i);
  fig_path   = [folder_name 'fig' num2str(i) '.png'];
  saveas(fig_handle, fig_path);
end

end