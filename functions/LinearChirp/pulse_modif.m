function p = pulse_modif(p, property, value)
% Modify the property of a linear chirp pulse
%
% Input:
%   - p, pulse
%   - property, the property of the pulse to be modified
%       - list of possible properties: cf. required and optional properties
%       in LinearChirp documentation. For bw/tp/Q, w1 is modified. Q is 
%       modified if w1 is input. TBP and sm are recomputed and can be
%       modified when necessary.
%   - value, the value of the property
%
% Output:
%   - p, the modified pulse
%
% Warning: the pulse is rewritten from scratch, pulse tweaks since it
% previous creation on the coordinates vectors will be lost. The same
% can apply to the properies when changing type.

grumble(p, property);

par = p; % generating new parameters from input pulse
par = rmfield(p,{'TBP' 't' 'Cx' 'Cy' 'Pr' 'Pph'});

% prioritize modification of w1 or Q
if property == "w1"
    par = rmfield(par, 'Q');
else
    par = rmfield(par, 'w1');
end

% type conversion
if property == "type"
   if value == "sinsmoothed"
       
       if par.delta_f ~=0
          warning('delta_f is lost during conversion to sinsmoothed')
       end
       
       par = rmfield(par, 'delta_f');
       
       par = rmfield(par, 'n');
       
   end
end

% modifying the property value
par.(property) = value;

% new pulse
p = LinearChirp(par);

end

function grumble(p, property)


if ~isstruct(p)
    error('p should be a pulse structure')
end
    
if ~isstring(property)
    error('property should be a string')
end
    
end