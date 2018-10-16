function [t, protein_conc] = alpha_pheromone_system(start_conc, tEnd, Q)
    % For p28: Q = (120/28) / ((120 / 28) + 1)
    % For myrosinase: Q = (120 / 527) / ((120 / 527) + 1)
    
    % Initial conditions 
    x0 = zeros(38, 1);
    x0(1) = start_conc; % alpha-factor 
    x0(2) = 1666.6666667; % Ste2
    x0(3) = 0.0; % Ste2active
    x0(4) = 1666.6666667; % Gabc
    x0(5) = 0.0; % GaGTP
    x0(6) = 0.0; % Gbc
    x0(7) = 0.0; % GaGDP
    x0(8) = 235.724935791903; % complex C
    x0(9) = 0.0; % complex D
    x0(10) = 158.33176608789; % Ste5
    x0(11) = 158.33176608789; % Ste11
    x0(12) = 105.943298120207; % complex A
    x0(13) = 36.3997016405141; % Ste7
    x0(14) = 686.399701640513; % Fus3
    x0(15) = 77.8753625675829; % complex B
    x0(16) = 1000.0; % Ste20
    x0(17) = 0.0; % complex E
    x0(18) = 0.0; % complex F
    x0(19) = 0.0; % complex G
    x0(20) = 0.0; % complex H
    x0(21) = 0.0; % complex I
    x0(22) = 0.0; % complex L
    x0(23) = 0.0; % Fus3PP
    x0(24) = 0.0; % complex K
    x0(25) = 200.0; % Ste12
    x0(26) = 0.0; % Ste12active
    x0(27) = 0.0; % Bar1
    x0(28) = 0.0; % Bar1active
    x0(29) = 0.0; % Bar1activeEx
    x0(30) = 0.0; 500.0; % Far1
    x0(31) = 0.0; % Far1PP
    x0(32) = 0.0; % Far1U
    x0(33) = 0.0; % complex M
    x0(34) = 0.0; % complex N
    x0(35) = 300.0; % Cdc28
    x0(36) = 0.0; % Sst2
    % The following species have been added
    x0(37) = 0.0; % mRNA of anti-cancer agent
    x0(38) = 0.0; % anti-cancer agent
    
    tI = [0 tEnd];
    [t, x] = ode23s(@(t, x) f_alpha_pheromone_system(t, x, Q), tI, x0);

    protein_conc = x(:, 38);

end


