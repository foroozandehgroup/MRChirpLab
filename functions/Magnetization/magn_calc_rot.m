function magn_fin = magn_calc_rot(pulses, trec, ph_cy, offsets)
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


np_offs = length(offsets);
magn_fin_pc = zeros(3, np_offs, length(ph_cy(1,:)));

offsets = 2*pi*offsets;

parfor npc = 1:length(ph_cy(1,:)) % phase cycling loop

    for noffs = 1:np_offs % offest loop
        
        magn = [0;0;1];
        tend = 0;
        
        for i = 1:length(pulses)
            
            % potential delay between pulses 
            if pulses{i}.t(1) > tend + pulses{i}.tres
                delay = pulses{i}.t(1) - tend;
                magn = Rz(offsets(noffs) * delay) * magn;
            end
            
            % pulse phase cycling and apply pulse
            p = pulse_phase_correction(pulses{i}, ph_cy(i, npc));
            
            for m = 1:p.np
                magn = Rrod(2*pi*p.Cx(m), 2*pi*p.Cy(m), offsets(noffs), p.tres) * magn;
            end
            
            % end time after last pulse
            tend = pulses{i}.t(end);
        end
        
        % potential delay at the end of the pulsesuence
        if tend < trec
            magn = Rz(offsets(noffs) * (trec - tend)) * magn;
        end
        
        % receiver phase
        magn_fin_pc(:, noffs, npc) = Rz(-ph_cy(end, npc)) * magn;
    end

end

% phase cycling sum
magn_fin = sum(magn_fin_pc, 3) / length(ph_cy(1,:)); 

end












