function p = pulse_sum(p1, p2)
% Creates a pulse from the sum of 2 pulses Cartesian coordinates
%
% Input:
%   - p1 and p2, the pulses to be added. p2 should not be positioned before
%   p1
%
% Output
%   - p, new pulse made from the sum of p1 and p2 Cartesian coordinates


grumble(p1, p2)

%% initialization

p = struct();
p.tres = p1.tres;

% time vector
p.t = p1.t(1):p.tres:p2.t(end);

%% summing cartesian coordinates

% p1 and p2 Cartesian coordinates indexes
j = 1;
k = 1;

p.Cx = [];
p.Cy = [];

for i = 1:length(p.t)
   
   p.Cx(i) = 0;
   p.Cy(i) = 0;
    
   %if abs(p.t(i)-p1.t(j)) < 100*eps(p.t(i)) % float comparison p.t(i)==p1.t(j)
   if (p1.delta_t - p1.tp/2 <= p.t(i)) && (j < length(p1.Cx))
       
       p.Cx(i) = p.Cx(i) + p1.Cx(j);
       p.Cy(i) = p.Cy(i) + p1.Cy(j);
       j = j+1;
   end
   
   %if abs(p.t(i)-p2.t(k)) < 100*eps(p.t(i))
   if (p2.delta_t - p2.tp/2 <= p.t(i)) && (k < length(p2.Cx))   
       
       p.Cx(i) = p.Cx(i) + p2.Cx(k);
       p.Cy(i) = p.Cy(i) + p2.Cy(k);
       k = k+1;
   end
end

%% other pulse properties

% polar coordinates
[p.Pph, p.Pr] = cart2pol(p.Cx, p.Cy);

% number of points
p.np = length(p.t);

% duration (tres correction as p.t contains the center of the pulse)
p.tp = p.t(end)-p.t(1) + p.tres;

% position
p.delta_t = p.t(end) + p.tres/2 - p.tp/2;

% phase
p.phi0 = 0;

end

function grumble(p1, p2)

if p1.tres ~= p2.tres
    error('p1 and p2 should have the same tres')
end

if p1.delta_t > p2.delta_t
   error('p2 should not be positioned after p1') 
end

end