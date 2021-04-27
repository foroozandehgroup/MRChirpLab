function MRchirp_test()

par.tres = 0.5e-6;

par.bw = 300000;
par.tp = 500e-6;

% not enough inputs test
try
    LinearChirp(par);
catch error
    disp(error.message)
end

phase = ["superGaussian", "chirp", "tanh"];
amp = ["superGaussian", "sinsmoothed", "WURST", "sech"];

for i = 1:length(phase)
    for j = 1:length(amp)
                
        par.phase = phase(i);
        par.amp = amp(j);

        par.bw = 300000;
        par.tp = 500e-6;
        par.w1 = 6.491401671720040e+03;

        % minimal definition
        p0 = MRchirp(par);
        
        if phase(i) == "chirp" || phase(i) == "superGaussian"
            % Q defined - no bw
            par.Q = par.w1^2 * 2 * pi * par.tp / par.bw;
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
        end
        
        if phase(i) == "tanh"
            
            % Q defined - no bw
            par.k = 0.3747813;
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

            par = rmfield(par, 'k');

            % should be the same
            pulses = {p0 p1 p2 p3};
            titles = ["3 parameters" "k defined no bw" ...
                      "k defined no w1"  "k defined no tp"];
        end
        
        
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
p.sm = 20;
p_sinsmoothed_sm20 = MRchirp(par);

par.amp = "WURST";
p_WURST_n20 = MRchirp(par);
p.n = 80;
p_WURST_n80 = MRchirp(par);

par.amp = "superGaussian";
p_sG_n40 = MRchirp(par);
p.n = 20;
p_sG_n20 = MRchirp(par);

pulses = {p_sinsmoothed_sm10 p_sinsmoothed_sm20 ...
          p_WURST_n20 p_WURST_n80 ...
          p_sG_n40 p_sG_n20};
titles = ["p_sinsmoothed_sm10" "p_sinsmoothed_sm20" ...
          "p_WURST_n20" "p_WURST_n80" ...
          "p_sG_n40" "p_sG_n20"];
plot_pulse(pulses, "", titles);

% test warning for high number of points
par.tres = 5e-9;
pulse_low_tres = LinearChirp(par);

end