function MRchirp_test()

par.tres = 0.5e-6;
par.bw = 300000;
par.tp = 500e-6;

% not enough inputs test
disp(' ')
disp('--- Test warning for unsufficient number of input ---')
try
    MRchirp(par);
catch error
    disp(error.message)
end

phase = ["chirp", "tanh"];
amp = ["superGaussian", "sinsmoothed", "linearsmoothed", "WURST", "sech"];

for i = 1:length(phase)
    for j = 1:length(amp)
                
        par.phase = phase(i);
        par.amp = amp(j);

        par.bw = 300000;
        par.tp = 500e-6;
        par.w1 = 6.491401671720040e+03;

        % minimal definition
        p0 = MRchirp(par);
        
        % Q defined - no bw
        if phase(i) == "chirp"
            par.Q = par.w1^2 * 2 * pi * par.tp / par.bw;      
        elseif phase(i) == "tanh"
            par.Q = 4 * pi * par.w1^2 / (par.bw * (10.6/par.tp));
        end
        par = rmfield(par, 'bw');
        p1 = MRchirp(par);

        % Q defined1 - no w1
        par.bw = 300000;
        par = rmfield(par, 'w1');
        p2 = MRchirp(par);
        par.w1 = 6.491401671720040e+03;

        % Q defined - no tp
        par = rmfield(par, 'tp');
        p3 = MRchirp(par);
        par.tp = 500e-6;

        par = rmfield(par, 'Q');

        % should be the same
        pulses = {p0 p1 p2 p3};
        titles = ["3 parameters" "Q defined no bw" ...
                  "Q defined no w1"  "Q defined no tp"];
                
        plot_pulse(pulses, "",titles);
        sgtitle("amp=" + p0.amp + " - phase=" + p0.phase)
        
        % opional parameters
        par.delta_t = 500e-6;
        p_delta_t_500us = MRchirp(par);
        par = rmfield(par, 'delta_t');

        par.phi0 = pi;
        p_phi0_180deg = MRchirp(par);
        par = rmfield(par, 'phi0');

        par.delta_f = par.bw/4;
        p_delta_f = MRchirp(par);
        
        % reverse sweep
        par.bw = - par.bw;
        p_reverse_sweep = MRchirp(par);
        
        par = rmfield(par, 'delta_f');
        
        pulses = {p0 p_delta_t_500us p_phi0_180deg p_delta_f p_reverse_sweep};
        titles = ["No optional parameters" "delta_t = 500us" ...
                  "phi0 = 180deg"  "delta_f = 0.25*bw/4" ...
                  "delta_f = 0.25*bw/4 - reversed sweep"];
        plot_pulse(pulses, "", titles);
        sgtitle("amp=" + p0.amp + " - phase=" + p0.phase)

    end
end

% optional parameters specific to amp/phase values
par.phase = "chirp";
par.amp = "sinsmoothed";

p_sinsmoothed_sm10 = MRchirp(par);
par.sm = 20;
p_sinsmoothed_sm20 = MRchirp(par);

par.amp = "linearsmoothed";
p_linearsmoothed_sm20 = MRchirp(par);

par.amp = "WURST";
p_WURST_n20 = MRchirp(par);
par.n = 80;
p_WURST_n80 = MRchirp(par);

par.amp = "superGaussian";
p_sG_n40 = MRchirp(par);
par.n = 20;
p_sG_n20 = MRchirp(par);

pulses = {p_sinsmoothed_sm10 p_sinsmoothed_sm20 ...
          p_linearsmoothed_sm20 ...
          p_WURST_n20 p_WURST_n80 ...
          p_sG_n40 p_sG_n20};
titles = ["p_sinsmoothed_sm10" "p_sinsmoothed_sm20" ...
          "p_linearsmoothed_sm20" ...
          "p_WURST_n20" "p_WURST_n80" ...
          "p_sG_n40" "p_sG_n20"];
plot_pulse(pulses, "", titles);

% B/tp/k for HS pulses
par.phase = "tanh";
par.amp = "sech";
HS_tp = MRchirp(par);
par.k = 5.3;
HS_k_5 = MRchirp(par);
par = rmfield(par, 'tp');
par.B = 5.3/500e-6;
HS_B_k_5 = MRchirp(par);

pulses = {HS_tp HS_k_5 HS_B_k_5};
titles = ["HS_tp" "HS_k_5" "HS_B_k_5"];
plot_pulse(pulses, "", titles);

disp(' ')
disp('--- Test warning for high number of points ---')
par.tres = 5e-9;
pulse_low_tres = MRchirp(par);
disp(' ')

end