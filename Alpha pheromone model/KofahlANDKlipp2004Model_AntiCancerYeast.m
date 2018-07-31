function [t,alpha_factor_conc] = KofahlANDKlipp2004Model_AntiCancerYeast(reaction_v48_k48,pop_size)
    % Is it the alpha concentration that we want or what did we measure?
    % SHOULD INSTEAD RETURN ANTI-CANCER AGENT CONCENTRATION
    % Initial conditions and constants from article 
        
    nr_species = 38;
    % A: Try different population sizes (but we need to have the populatio size 
    % when we do the measurements or will the alpha system be used in the measurements?7)
    baseline_conc = 100.0; % value from ready-made model, but 1000 nM in article?
    % A: Should find our own value; expression of GFP without alpha pheromone

    % Initial conditions vector
    x0 = zeros(nr_species, 1);
    x0(1) = baseline_conc; % alpha-factor 
    x0(2) = 1666.6666667; % Ste2
    x0(3) = 0.0; % Ste2active
    x0(4) = 1666.6666667; % Gabc
    x0(5) = 0.0; % GaGTP
    x0(6) = 0.0; % Gbc
    x0(7) = 0.0; % GaGDP
    x0(8) = 235.724935791903; % complexC
    x0(9) = 0.0; % complexD
    x0(10) = 158.33176608789; % Ste5
    x0(11) = 158.33176608789; % Ste11
    x0(12) = 105.943298120207; % complexA
    x0(13) = 36.3997016405141; % Ste7
    x0(14) = 686.399701640513; % Fus3
    x0(15) = 77.8753625675829; % complexB
    x0(16) = 1000.0; % Ste20
    x0(17) = 0.0; % complexE
    x0(18) = 0.0; % complexF
    x0(19) = 0.0; % complexG
    x0(20) = 0.0; % complexH
    x0(21) = 0.0; % complexI
    x0(22) = 0.0; % complexL
    x0(23) = 0.0; % Fus3PP
    x0(24) = 0.0; % complexK
    x0(25) = 200.0; % Ste12
    x0(26) = 0.0; % Ste12active
    x0(27) = 0.0; % Bar1
    x0(28) = 0.0; % Bar1active
    x0(29) = 0.0; % Bar1activeEx
    x0(30) = 500.0; % Far1
    x0(31) = 0.0; % Far1PP
    x0(32) = 0.0; % Far1U
    x0(33) = 0.0; % complexM
    x0(34) = 0.0; % complexN
    x0(35) = 300.0; % Cdc28
    x0(36) = 0.0; % Sst2
    x0(37) = 0.0; % p
    % The following species have been added
    x0(38) = baseline_conc; % anti-cancer agent
    

    % A: This is never used? What is it? Steady-state concentrations?
    xss = [2.4703e-32
       -2.696e-29
      -1.0794e-30
           1666.7
      -2.7368e-29
       1.5853e-08
      -1.0457e-25
           239.08
       8.0768e-11
           155.64
           155.64
           105.28
           35.222
           685.22
           75.697
             1000
        2.524e-08
       4.8538e-09
        6.518e-10
       4.0886e-09
       1.9946e-06
       4.0886e-08
       5.5849e-06
       1.1162e-06
              200
        0.0020105
      -1.6914e-16
      -3.0921e-19
              200
      -1.1707e-63
      -7.0238e-63
              500
      -1.1505e-69
      -2.1773e-59
              300
       3.8989e-10
                0];

    %Creating linespace
    % t=linspace(0,9000);
    tI = [0 250];
    %Solving equations
    %lsode replaced by ODE45
    %ODE45 is equivalent of lsode in matlab

    % A: Add threshold
    alpha_threshold = 0.0; % Find in literature?
    % A: Assume the threshold is for the Ste2 receptor in the sense that
    % the MAPK cascade is not activated if the concentration of alpha
    % pheromone is too low
    if (pop_size * baseline_conc) > alpha_threshold
        % Add extra input value
        disp('Running model')
        [t,x] = ode23s(@(t, x) f(t, x, reaction_v48_k48, nr_species, pop_size), tI, x0);
    end 

    plot(t,x(:,1)); % Alpha pheromone
    

    % plot(t,x(:,26));
    % title('Analysis of the cell cycle arrest knockout on the transcription factor Ste12');
    % ylabel('[nM]');xlabel('Time(min)');
    % legend('Ste12active concentration with knockout')
    % A: DID THEY REMOVE SOMETHING???? SHOULD REALLY COMPARE THE MODEL WITH
    % THE ARTICLE
    

    alpha_factor_conc = x(:,1); 
end


