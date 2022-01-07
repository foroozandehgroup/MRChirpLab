function seq_Xepr_shape_files(seq, ID, IDnb, path, resonator)
% Saves a sequence pulses as Xepr shape files in a directory
%
% Input:
%   - seq, sequence to be saved
%   - ID, the identification of the sequence
%   - IDnb, the identification number of the sequence. The identification
%   number of the pulses are attributed from this one.
%   - path, location at which to save the files
%
% Files generated in created directory:
%   - pulses shapes files
%   - text file with the sequence main parameters
%   - figures opened when the function is called (are assumed to be linked
%   to the sequence). Plots of each pulse are previoulsy generated.


if isfield(seq, 'bw') && isfield(seq, 'total_time')
    total_time = ['_' num2str(seq.total_time * 1e3) 'ms'];
    bw = ['_bw_' num2str(seq.bw * 1e-3) 'kHz'];
    ID = [ID bw total_time];
end

% directory creation at path
folder_name = [path IDnb '_' ID];

if exist(folder_name, 'dir') == 7   
   delete([folder_name '\*'])
   rmdir(folder_name)
   disp("Already existing directory deleted")
end

mkdir(folder_name)
folder_name = [folder_name '\'];

% shape file of each pulse
ph_cy = ["000" "090" "180" "270"];
IDnb_pulse = IDnb;
additional_information = '';
for i = 1:length(seq.pulses)
    for k = 1:length(ph_cy)
        
        % creating a version of each pulse for phase cycling
        % negative phase! (Xepr convention)
        cycle_phase = - deg2rad(str2num(ph_cy(k)));
        pulse2save = pulse_phase_correction(seq.pulses{i}, cycle_phase);
        
        if nargin == 5
            pulse2save = pulse_resonator_easyspin(pulse2save, ...
                                resonator.f, resonator.H_f, resonator.nu);
            pulse2save.type = pulse2save.phase + "_resonator compensation";
            
            additional_information = ['resonator compensation centrered'...
                                      ' at ' num2str(resonator.nu) 'GHz'];   
        end
        
        
        ID_pulse = [num2str(i) '_ph' char(ph_cy(k))];
        
        pulse_Xepr_shape_file(pulse2save, ID_pulse, IDnb_pulse, folder_name)
        
        IDnb_pulse = num2str(str2num(IDnb_pulse) + 1);
        plot_pulse(pulse2save, "Xepr", convertCharsToStrings(ID_pulse))
    end
end

% potential add:
% for loop over created directory to concatenate all the pulses shape files
% into one

% save sequence parameters into a text file in the created directory
seq_param2txt(seq, ID, IDnb, folder_name, additional_information)

% all the opened figures assumed to be linked to the sequence
disp('Saving all open figures...')

fig_list = findobj(allchild(0), 'flat', 'Type', 'figure');
for i = 1:length(fig_list)
  fig_handle = fig_list(length(fig_list)+1-i);
  fig_path   = [folder_name 'fig' num2str(i) '.png'];
  saveas(fig_handle, fig_path);
end

end

