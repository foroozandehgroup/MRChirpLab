function residual_time = chorus_condition_tolerance(C, x)
% display the residual evolution times of a pulse sequence
%
% Input:
%
%   C, condition matrix of the sequence
%   x, vector containing the times (in s)
%   np, number of points in the offset vector (generally 1000)

y = C * x;

residual_time = max(abs(y(:,1)));

end

