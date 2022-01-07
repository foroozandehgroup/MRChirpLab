function seq = chorus_cpmg(param)
% Creates a CPMG chorus pulse sequence
% 
% Input:
%   - param, a structure containing the parameters allowing to defined the
%   sequence
%
% Output:
%   - seq, a structure containing the sequence information
%
% Required fields for param:
%   - requirements which apply for chorus (cf. chorus function 
%   documentation)
%
% Optional fields for param:
%   - tau_echo, a positive number which represent the delay to be added on
%   each side of the refocusing block
%   - n, the number CPMG elements (CHORUS refocusing counted as 1 and (n-1)
%   loop

%%  default values for optional parameters

grumble(param) % additional checks

if ~isfield(param, 'tau_echo')
    param.tau_echo = 0;
end

if ~isfield(param, 'n')
    param.n = 2;
end

if ~isfield(param, 'display_result')
    param.display_result = false;
end

%% sequence structure 
% chorus - [- 180 - tau_echo - 180 - tau_echo ] * (n-1)

% avoid warning messages from chorus funtion
par = rmfield(param, {'n', 'tau_echo'});
par.display_result = false;

% chorus (excitation)
seq = exc_3fs(par);

seq.n = param.n;
seq.tau_echo = param.tau_echo;

% appending inversion n inversion blocks to chorus
for i = 1:param.n-1
        
    % 1st 180 pulse
    seq = seq_add_pulse(seq, seq.pulses{3});
    
    % tau_echo
    seq.tau = [seq.tau param.tau_echo];
    
    % 2nd 180 pulse
    seq = seq_add_pulse(seq, seq.pulses{3});
    
    % tau_echo
    seq.tau = [seq.tau param.tau_echo];
    
end

seq.total_time = sum(seq.tau);

% suggested phase cycling (64 steps - 2 steps on refocusing 180s)

ph1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ph2 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph3 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];

seq.pc = [ph1; ph2; ph3];

for i = 1:param.n-1

    ph4 = [1 1 1 1 3 3 3 3 1 1 1 1 3 3 3 3 1 1 1 1 3 3 3 3 1 1 1 1 3 3 3 3 1 1 1 1 3 3 3 3 1 1 1 1 3 3 3 3 1 1 1 1 3 3 3 3 1 1 1 1 3 3 3 3];
    ph5 = [1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3 1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3];

    if i > 1
       ph4 = 2*(i-1)*ph4;
       ph5 = 2*(i-1)*ph5;
    end
    
    seq.pc = [seq.pc; ph4; ph5];
    
end

phrec = [0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 2 0 2 0 2 0 2 0 2 0 2 0 2 0 2 0];
seq.pc = pi/2 * [seq.pc; phrec];

if param.display_result == true
    
    plot_seq(seq);
    
    % offsets
    n_offs = 100;
    offs = linspace(-seq.bw/2, seq.bw/2, n_offs);
    
    opt.pc = seq.pc;
        
    % magnetization simulation
    final_magn = magn_calc_rot(seq.pulses, seq.total_time, offs, opt);
    plot_magn(final_magn,offs)
    
end

end

function grumble(param)

if isfield(param, 'tau_echo')
    if ~isreal(param.tau_echo) || param.tau_echo < 0
        error('tau_echo must be a positive real number')
    end
end

end