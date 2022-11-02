function magn_fin = magn_calc_rot_full(pulses, tcalc, ph_cy, offsets)
% Calculates the magnetization over the times of tcalc after applying a
% list of pulses
%
% Input: 
%   - pulses, cell array containing pulses
%   - tcalc, the time vector onto wich to calculate the magnetization
%   - ph_cy, phase cycle  applied to the pulses. Array of dimension
%   phase cycle steps*pulses length. If no phase cycle is used, use a
%   column vector of zeros, i.e. zeros(1,length(pulses))
%   - offsets, offset vector onto which to calculate the magnetization
%
% Output:
%   - magn_fin, calculated magnetization which is output for the different 
%   offsets of off and time of tcalc in a 3D-matrix.


grumble(tcalc)

np_offs = length(offsets);

magn_fin_pc = zeros(3, length(tcalc), np_offs, length(ph_cy(1,:)));

for npc = 1:length(ph_cy(1,:)) % phase cycling loop

    for noffs = 1:np_offs % offest loop
        
        t = 0;   % time
        i_t = 1; % time iterator
    
        magn_fin_pc(:, i_t, noffs, npc) = [0;0;1]; % initial magnetization on z
        
        for i = 1:length(pulses)
            
            % potential delay between pulses 
            if pulses{i}.t(1) > t + pulses{i}.tres
                while pulses{i}.t(1) > t + pulses{i}.tres
                    t = t + pulses{i}.tres;
                    i_t = i_t + 1;
                    magn_fin_pc(:, i_t, noffs, npc) = ...
                        Rz(2 * pi * offsets(noffs) * pulses{i}.tres) * ...
                        magn_fin_pc(:, i_t-1, noffs, npc);
                    
                end
            end
            			
            % apply pulse
            for m = 1:pulses{i}.np
                i_t = i_t + 1;
                magn_fin_pc(:, i_t, noffs, npc) = ...
                    Rrod(2*pi*pulses{i}.Cx(m), 2*pi*pulses{i}.Cy(m), 2*pi*offsets(noffs), pulses{i}.tres) * ...
                    magn_fin_pc(:, i_t-1, noffs, npc);
                t = t + pulses{i}.tres;
            end
%             
%             for m = 1:p.np
%                 magn =  * magn;
%                                     Rtot(pulses{i}.Pr(m), offsets(noffs), ...
%                          ph_cy(i, npc) + pulses{i}.Pph(m), ...
%                          pulses{i}.tres) * ...
%             end
            
        end
        
        % potential delay at the end of the pulse sequence
        while t < tcalc(end)
            i_t = i_t + 1;
            magn_fin_pc(:, i_t, noffs, npc) = ...
                Rz(2 * pi * offsets(noffs) * pulses{1}.tres) * ...
                magn_fin_pc(:, i_t-1, noffs, npc);
            t = t + pulses{i}.tres;
        end
        
        % receiver phase (applied at each time and offset)
        magn_fin_pc(:, :, noffs, npc) = Rz(-ph_cy(end, npc)) * ...
                                          magn_fin_pc(:, :, noffs, npc);
    end
end

% phase cycling sum
magn_fin = sum(magn_fin_pc, 4) / length(ph_cy(1,:)); 

end

function grumble(tcalc)

if size(tcalc,1) ~= 1 || size(tcalc,2) < 2
    error('tcalc must be an array.');
end
    
end
