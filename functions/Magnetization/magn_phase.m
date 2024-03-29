function ph = magn_phase(magn)
% Returns the phase of the magnetization magn as a vector ph
%
% Input:
%   - magn, magnetization on x, y and z for a certain number of offsets
%
% Output:
%   - ph, phase for the offsets of magn


% angle between x and y components
ph = angle(complex(magn(2,:), magn(1,:)));
ph = unwrap(wrapTo2Pi(ph));
   
end