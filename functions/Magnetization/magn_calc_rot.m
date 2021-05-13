function magn_end = magn_calc_rot(pulses, tend, off, opt)
% Calculates the magnetization at time tend after applying the list of
% chirps pulses
%
% Input: 
%   - pulses, cell array containing pulses
%   - tend, time at which to calculate the magnetization
%   - off, offset vector onto which to calculate the magnetization
%
% Output:
%   - magn_end, calculated magnetization which is output at x, y, and z for 
%   the different offsets of off in a 2D-matrix.
%
% Optional: 
%   - opt: field for the following options:
%       - magn_init, the initial magnetization, on +z by default
%       - tstart, the starting time of the calculation, 0 by default (s)
%       - pc, phase cycle applied to the pulses 
%       - B1, B1 variations coefficients vector

grumble(pulses, tend, off)

n_off = length(off);
off = 2 * pi * off; % radians conversion

% default values for start time and initial magnetization
if nargin > 3
    grumble_opt(pulses, opt)
    if ~isfield(opt, 'magn_init')
        opt.magn_init = repmat([0,0,1]', 1, n_off);
    end
    if ~isfield(opt, 'tstart')
        opt.tstart = 0;
    end
else
    opt.tstart = 0;
    opt.magn_init = repmat([0,0,1]', 1, n_off);
end

%% no phase cycling
if ~isfield(opt, 'pc') && ~isfield(opt, 'B1')
    
    magn_end = opt.magn_init;

    parfor ioff = 1:n_off % offest loop

        magn = magn_end(:, ioff);
        
        t = opt.tstart;
        
        for i = 1:length(pulses)

            % potential delay between pulses 
            if pulses{i}.delta_t - pulses{i}.tp/2 - 10*eps > t
                delay = (pulses{i}.delta_t - pulses{i}.tp/2) - t;
                magn = Rz(off(ioff) * delay) * magn;
            end

            p = pulses{i};

            for m = 1:p.np
                magn = Rrod(2*pi*p.Cx(m), 2*pi*p.Cy(m), off(ioff), p.tres) * magn;
            end

            % end time after last pulse
            t = pulses{i}.delta_t + pulses{i}.tp/2;
        end

        % potential delay at the end of the pulse sequence
        if t < tend
            magn = Rz(off(ioff) * (tend - t)) * magn;
        end

        % receiver phase
        magn_end(:, ioff) = magn;
    end
    
end

%% phase cycling
if isfield(opt, 'pc') && ~isfield(opt, 'B1')
    
    pc = opt.pc;
    n_pc = length(pc(1,:));
        
    magn_end_pc = repmat(opt.magn_init, 1, 1, n_pc);
    
    parfor ipc = 1:length(pc(1,:)) % phase cycling loop

        for ioff = 1:n_off % offest loop

            magn = magn_end_pc(:, ioff, ipc);
            
            t = opt.tstart;
            
            for i = 1:length(pulses)

                % potential delay between pulses 
                if pulses{i}.delta_t - pulses{i}.tp/2 - 10*eps > t
                    delay = (pulses{i}.delta_t - pulses{i}.tp/2) - t;
                    magn = Rz(off(ioff) * delay) * magn;
                end

                % pulse phase cycling and apply pulse
                p = pulse_phase_correction(pulses{i}, pc(i, ipc));

                for m = 1:p.np
                    magn = Rrod(2*pi*p.Cx(m), 2*pi*p.Cy(m), off(ioff), p.tres) * magn;
                end

                % end time after last pulse
                t = pulses{i}.delta_t + pulses{i}.tp/2;
            end

            % potential delay at the end of the pulse sequence
            if t < tend
                magn = Rz(off(ioff) * (tend - t)) * magn;
            end

            % receiver phase
            magn_end_pc(:, ioff, ipc) = Rz(-pc(end, ipc)) * magn;
        end

    end

    % phase cycling sum
    magn_end = sum(magn_end_pc, 3) / n_pc; 
end

%% no phase cycling and B1 variations
if ~isfield(opt, 'pc') && isfield(opt, 'B1')
    
    B1 = opt.B1;
    n_B1 = length(B1);

    magn_end = zeros(3, n_off, n_B1);
    
    for iB1 = 1:n_B1 % rf amplitudes loop

        for ioff = 1:n_off % offest loop

            magn = opt.magn_init(:, ioff);
            t = opt.tstart;

            for i = 1:length(pulses)

                % potential delay between pulses 
                if pulses{i}.delta_t - pulses{i}.tp/2 - 10*eps > t
                    delay = (pulses{i}.delta_t - pulses{i}.tp/2) - t;
                    magn = Rz(off(ioff) * delay) * magn;
                end

                p = pulses{i};

                % B1 variation
                p.Cx = p.Cx * B1(iB1);
                p.Cy = p.Cy * B1(iB1);

                for m = 1:p.np
                    magn = Rrod(2*pi*p.Cx(m), 2*pi*p.Cy(m), off(ioff), p.tres) * magn;
                end

                % end time after last pulse
                t = pulses{i}.delta_t + pulses{i}.tp/2;
            end

            % potential delay at the end of the pulse sequence
            if t < tend
                magn = Rz(off(ioff) * (tend - t)) * magn;
            end

            % receiver phase
            magn_end(:, ioff, iB1) = magn;
        end
    end    
end

%% phase cycling and B1 variations
if isfield(opt, 'pc') && isfield(opt, 'B1')
    
    pc = opt.pc;
    n_pc = length(pc(1,:));
    
    B1 = opt.B1;
    n_B1 = length(B1);
    
    magn_end_pc = zeros(3, n_off, n_B1, n_pc);

    parfor ipc = 1:n_pc % phase cycling loop

        for iB1 = 1:n_B1 % rf amplitudes loop

            for ioff = 1:n_off % offest loop

                magn = opt.magn_init(:, ioff);
                t = opt.tstart;

                for i = 1:length(pulses)

                    % potential delay between pulses 
                    if pulses{i}.delta_t - pulses{i}.tp/2 - 10*eps > t
                        delay = (pulses{i}.delta_t - pulses{i}.tp/2) - t;
                        magn = Rz(off(ioff) * delay) * magn;
                    end

                    % pulse phase cycling and apply pulse
                    p = pulse_phase_correction(pulses{i}, pc(i, ipc));
                    
                    % B1 variation
                    p.Cx = p.Cx * B1(iB1);
                    p.Cy = p.Cy * B1(iB1);

                    for m = 1:p.np
                        magn = Rrod(2*pi*p.Cx(m), 2*pi*p.Cy(m), off(ioff), p.tres) * magn;
                    end

                    % end time after last pulse
                    t = pulses{i}.delta_t + pulses{i}.tp/2;
                end

                % potential delay at the end of the pulse sequence
                if t < tend
                    magn = Rz(off(ioff) * (tend - t)) * magn;
                end

                % receiver phase
                magn_end_pc(:, ioff, iB1, ipc) = Rz(-pc(end, ipc)) * magn;
            end
        end
    end

    % phase cycling sum
    magn_end = sum(magn_end_pc, 4) / n_pc; 

end

end

function grumble(pulses, tend, off)

if ~iscell(pulses)
    error('pulses should be a cell array')
end

if ~isreal(tend) || tend < 0
    error('tend should be a positive real number')
end

if ~isvector(off) || ~isreal(off)
    error('off should be a vector of real numbers')
end

end

function grumble_opt(pulses, opt)

if isfield(opt, 'pc')
    s = size(opt.pc);
    if s(1) ~= length(pulses)+1
        error(['pc number of lines, length(pc(:,1)), must be equal ' ...
               'to the number of pulses + 1'])
    end
end

end









