function magn_fin = magn_calc_rot_Rtot(pulses, trec, ph_cy, offsets)
% Calculates the magnetization at time trec after applying the list of
% chirps pulses
%
% Input: 
%   - pulses, cell array containing pulses
%   - trec, time at which to calculate the magnetization
%   - ph_cy, phase cycle  applied to the pulses
%   - offsets, offset vector onto which to calculate the magnetization
%
% Output:
%   - magn_fin, calculated magnetization which is output at x,y, and z for 
%   the different offsets of off in a 2D-matrix.


disp('Magnetization computation...')

np_offs = length(offsets);
magn_fin_pc = zeros(3, np_offs, length(ph_cy(1,:)));

for npc = 1:length(ph_cy(1,:)) % phase cycling loop

    for noffs = 1:np_offs % offest loop
        
        magn = [0;0;1];
        tend = 0;
        
        for i = 1:length(pulses)
            
            % potential delay between pulses 
            if pulses{i}.t(1) > tend + pulses{i}.tres
                delay = pulses{i}.t(1) - tend;
                magn = Rz(2 * pi * offsets(noffs) * delay) * magn;
            end
            			
            % apply pulses
            for m = 1:pulses{i}.np
                magn = Rtot(pulses{i}.Pr(m), offsets(noffs), ...
                            ph_cy(i, npc) + pulses{i}.Pph(m), ...
                            pulses{i}.tres) * magn;
            end
            
            % end time after last pulse
            tend = pulses{i}.t(end);
        end
        
        % potential delay at the end of the pulsesuence
        if tend < trec
            magn = Rz(2 * pi * offsets(noffs) * (trec - tend)) * magn;
        end
        
        % receiver phase
        magn_fin_pc(:, noffs, npc) = Rz(-ph_cy(end, npc)) * magn;
    end
end

% phase cycling sum
magn_fin = sum(magn_fin_pc, 3) / length(ph_cy(1,:)); 

end












