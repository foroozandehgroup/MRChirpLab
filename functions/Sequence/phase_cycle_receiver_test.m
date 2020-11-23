function phase_cycle_receiver_test()

ph1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ph2 = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3];
ph3 = [0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3];
ph_rec = [0 2 0 2 2 0 2 0 0 2 0 2 2 0 2 0];
ctp = [-1, +2, -2];
disp(phase_cycle_receiver([ph1; ph2; ph3],ctp))

if ph_rec == phase_cycle_receiver([ph1; ph2; ph3],ctp)
    disp('success')
else
    disp('fail')
end

% test from function documentation
ph1 =    [0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2];
ph2 =    [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
ph3 =    [0 0 1 1 2 2 3 3 0 0 1 1 2 2 3 3];
ph_rec = [0 2 2 0 0 2 2 0 2 0 0 2 2 0 0 2];

ctp = [-1, +2, -2];

if ph_rec == phase_cycle_receiver([ph1; ph2; ph3],ctp)
    disp('success')
else
    disp('fail')
end

end






