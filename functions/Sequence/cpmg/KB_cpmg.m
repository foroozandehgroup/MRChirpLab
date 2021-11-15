function seq =  KB_cpmg(param)
% Creates a KB CPMG pulse sequence
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
% doublechirp - [ tau_echo/2 - 180 - tau_echo - 180 - tau_echo/2 ] * (n-1)

% avoid warning messages from doublechirp function
par = rmfield(param, {'n', 'tau_echo'});
par.display_result = false;

% double chirp excitation
seq = exc_2fs(par);

seq.n = param.n;
seq.tau_echo = param.tau_echo;

% appending inversion n inversion blocks to chorus
for i = 1:param.n-1

    % 1st 180 pulse
    seq = seq_add_pulse(seq, seq.pulses{2});
    
    % tau_echo
    seq.tau = [seq.tau param.tau_echo];
    
    % 2nd 180 pulse
    seq = seq_add_pulse(seq, seq.pulses{2});
    
    % tau_echo
    seq.tau = [seq.tau param.tau_echo];
    
end

seq.total_time = sum(seq.tau);

% suggested phase cycling (16 steps - 2 steps on refocusing 180s)

ph1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ph2 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];

seq.pc = [ph1; ph2];

for i = 1:param.n-1
    
    ph3 = [1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3];
    ph4 = [1 1 3 3 1 1 3 3 1 1 3 3 1 1 3 3];
    
    if i > 1
       ph3 = 2*(i-1)*ph3;
       ph4 = 2*(i-1)*ph4;
    end
    
    seq.pc = [seq.pc; ph3; ph4];
    
end

phrec = [0 0 0 0 2 2 2 2 0 0 0 0 2 2 2 2];
seq.pc = pi/2 * [seq.pc; phrec];

if param.display_result == true
    
    plot_seq(seq);
    
    % offsets
    n_offs = 500;
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