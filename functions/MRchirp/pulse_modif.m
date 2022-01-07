function p = pulse_modif(p, property, value)
% Modify the property of a phase-modulated pulse
%
% Input:
%   - p, pulse
%   - property, the property of the pulse to be modified
%       - list of possible properties: cf. required and optional properties
%       in MRchirp documentation. For bw/tp/Q/k, w1 is modified. Q/k is 
%       modified if w1 is input. TBP and sm are recomputed and can be
%       modified when necessary.
%   - value, the value of the property
%
% Output:
%   - p, the modified pulse
%
% Warning: the pulse is rewritten from scratch. Tweaks since its previous 
% creation on the coordinates vectors will be lost.

grumble(p, property);

par = p; % generating new parameters from input pulse
par = rmfield(p,{'np' 'TBP' 't' 'Cx' 'Cy' 'Pr' 'Pph'});

% prioritize modification of w1 or Q/k
if property == "w1"
    if p.phase == "chirp" || p.phase == "tanh"
        par = rmfield(par, 'Q');
    end
else
    par = rmfield(par, 'w1');
end

if p.amp == "sech" || p.phase == "tanh"
    par = rmfield(par, 'B');
end

% modifying the property value
par.(property) = value;

% new pulse
p = MRchirp(par);

% type conversion: remove obsolete fields
% TBD

end

function grumble(p, property)


if ~isstruct(p)
    error('p should be a pulse structure')
end
    
if ~isstring(property)
    error('property should be a string')
end
    
end