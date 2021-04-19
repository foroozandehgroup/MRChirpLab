function [value, name] = tp_w1_bw_Q(param)
% Computes the missing parameter of the equation for linear chirp:
% 2 * pi * tp * w1^2 = bw * Q
%
% Input 3 of the following:
%   - tp, duratation of the linear chirp (s)
%   - w1, radiofrequency amplitude of the chirp (Hz)
%   - bw, bandwidth of the chirp (Hz)
%   - Q, adiabaticity factor of the chirp
%
% Output:
%   - value, calculated parameter value
%   - name, calculated parameter name (optional)
%
% The calculated parameter one of the four possible input which was not
% given

grumble(param)

if isfield(param,'tp')
    tp = param.tp;
end

if isfield(param, 'bw')
    bw = param.bw;
end

if isfield(param, 'w1')
    w1 = param.w1;
end

if isfield(param, 'Q')
    Q = param.Q;
end

% last parameter computation
if ~isfield(param,'w1')
    name = 'w1';
    value = sqrt(bw * Q / (2 * pi * tp));
elseif ~isfield(param,'tp')
    name = 'tp';
    value = bw * Q / (2 * pi * w1^2);
elseif ~isfield(param,'bw')
    name = 'bw';
    value = w1^2 * 2 * pi * tp / Q;
elseif ~isfield(param, "Q")
    name = 'Q';
    value = w1^2 * 2 * pi * tp / bw;
end

end

function grumble(param)

if sum([isfield(param,'bw') isfield(param,'w1') ...
    isfield(param,'tp') isfield(param,'Q')]) ~= 3
    error(['3 of the following parameters must be defined: ' ...
           'bw w1 tp Q'])
end

if isfield(param,'tau_p')
    if ~isreal(param.tau_p)
    error('tp must be a real number')
    end
end

if isfield(param,'w1')
    if ~isreal(param.w1) || param.w1 < 0
    error('w1 must be a positive real number')
    end
end

if isfield(param,'bw')
    if ~isreal(param.bw) || param.bw < 0
    error('bw must be a positive real number')
    end
end

if isfield(param,'Q')
    if ~isreal(param.Q) || param.Q < 0
    error('Q must be a positive real number')
    end
end

end























