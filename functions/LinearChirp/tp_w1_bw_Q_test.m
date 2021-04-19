function tp_w1_bw_Q_test()

p.tp = 500e-6; p.w1 = 6.491401671720040e+03; p.bw = 300e3;

[Q, name] = tp_w1_bw_Q(p)
p.Q = tp_w1_bw_Q(p)

p = rmfield(p, 'bw');
[bw, name] = tp_w1_bw_Q(p)
p.bw = tp_w1_bw_Q(p)

p = rmfield(p, 'w1');
[w1, name] = tp_w1_bw_Q(p)
p.w1 = tp_w1_bw_Q(p)

p = rmfield(p, 'tp');
[tp, name] = tp_w1_bw_Q(p)
p.tp = tp_w1_bw_Q(p)


end























