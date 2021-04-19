function pulse_Xepr_shape_file(pulse, ID, IDnb, path)
% Creates and save the pulse shape in a .shp file for Xepr use
%
% Input:
%   - pulse, pulse whose shape is to be saved
%   - ID, identification of the pulse
%   - IDnb, the identification number of the pulse
%   - path, the path at which to save the pulse shape for Xepr
%
% File creation:
%   - IDnb_ID.shp, shape file

filename = [path IDnb '_' ID '.shp'];

header = ['begin shape' IDnb ' "' ID '"'];
dlmwrite(filename, header, 'delimiter', '', 'newline', 'pc')

shape = pulse_Xepr_format(pulse);
dlmwrite(filename, shape, '-append', 'precision', '%.6f', 'newline', 'pc')

footer = ['end shape' IDnb];
dlmwrite(filename, footer,'-append', 'delimiter', '', 'newline', 'pc')

end