function pulse_TopSpin_file(p, ID, path)
% Creates and save the pulse shape in a shape file for use with TopSpin
%
% Input:
%   - p, pulse whose shape is to be saved
%   - ID, identification of the pulse
%   - IDnb, the identification number of the pulse
%   - path, the path at which to save the pulse shape for Xepr
%
% File creation:
%   - IDnb_ID.shp, shape file

filename = [path ID '.txt'];

shape = pulse_TopSpin_format(p);
%     ' ; % to be smoothed ' num2str(p.sm) 13 10 ...
header = [...
    '##TITLE= ' ID 13 10 ...
    '##JCAMP-DX= 5.00 Bruker JCAMP library' 13 10 ...
    '##DATA TYPE= Shape Data' 13 10 ...
    '##ORIGIN= Bruker BioSpin GmbH' 13 10 ...
    '##OWNER= <M.FOROOZANDEH>' 13 10 ...
    '##DATE= ' num2str(datestr(now,1)) 13 10 ...
    '##TIME= ' num2str(datestr(now,'HH:MM:SS')) 13 10 ...
    '##$SHAPE_PARAMETERS= Amplitude/Phase modulation: ' char(p.amp) '/' char(p.phase) ' ; Total Sweep-Width [kHz] ' ...
    num2str(p.bw*1e-3) ' ; Length of Pulse [msec] ' num2str(p.tp*1e3) 13 10 ...
    '##MINX= ' num2str(min(shape(:,1))) 13 10 ...
    '##MAXX= ' num2str(max(shape(:,1))) 13 10 ...
    '##MINY= ' num2str(min(shape(:,2))) 13 10 ...
    '##MAXY= ' num2str(max(shape(:,2))) 13 10 ...
    '##$SHAPE_EXMODE= Adiabatic' 13 10 ...
    '##$SHAPE_TOTROT= 1.800000E02' 13 10 ...
    '##$SHAPE_TYPE= Inversion' 13 10 ...
    '##$SHAPE_USER_DEF= ' 13 10 ...
    '##$SHAPE_REPHFAC= ' 13 10 ...
    '##$SHAPE_BWFAC= 2.735620E02' 13 10 ...
    '##$SHAPE_BWFAC50= ' 13 10 ...
    '##$SHAPE_INTEGFAC= 3.236043E-02' 13 10 ...
    '##$SHAPE_MODE= 0' 13 10 ...
    '##NPOINTS= ' num2str(p.np) 13 10 ...
    '##XYPOINTS= (XY..XY)' 13 10];

footer = '##END=';

dlmwrite(filename, header, 'delimiter', '', 'newline', 'pc')
dlmwrite(filename, shape, '-append', 'delimiter',',','precision', '%.6f', 'newline', 'pc')
dlmwrite(filename, footer,'-append', 'delimiter', '', 'newline', 'pc')

end

