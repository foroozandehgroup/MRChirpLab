function C = chorus_condition_matrix(np)
% CHORUS pulse sequence durations equation under matrix form
%
% np (generally 1000) is the number of points in the offset vector alpha
  
alpha = linspace(0,1,np);

for i = 1:np
    % whole chemical shift
    C(i,:)=  [ 1-alpha(i), 2*alpha(i)-1, -1, -(2*alpha(i)-1)]; 

end


end

