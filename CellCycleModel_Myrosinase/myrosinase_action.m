function proliferation = myrosinase_action(myrosinase_0, mpf_peak, tf)

    Casp3_0 = 1.0e-5;

    % Myrosinase
    x0_m = zeros(3, 1);
    x0_m(1) = myrosinase_0; % Myrosinase
    x0_m(2) = 0.0; % Sulforaphane extracellular
    x0_m(3) = 0.0; % Sulforapane intracellular

    % p53 system
    x0_p = zeros(4, 1);
    x0_p(1) = 0.1403; % p53
    x0_p(2) = 3.336e-5; % Mdm2
    x0_p(3) = 0.5219; % Cop1
    x0_p(4) = 3.414e-3; % p53/Cop1

    % Cell cycle arrest
    x0_c = zeros(15, 1);
    x0_c(1) = 0.25e-6; % Cdc25 active
    x0_c(2) = 0.0; % Cdc25Ps216 active
    x0_c(3) = 0.25e-6; % Cdc25 inactive
    x0_c(4) = 0.5e-5; % Cdc25Ps216 inactive
    x0_c(5) = 2.5e-7; % Chk1P
    x0_c(6) = 0.5; % protein 14-3-3
    x0_c(7) = 0.0075; % Cdc25Ps216/14-3-3 inactive
    x0_c(8) = 0.0; % Wee1 phosphorylated
    x0_c(9) = 0.0; % p21
    x0_c(10) = 0.25e-8; % MPF
    x0_c(11) = 0.0; % p21/MPF
    x0_c(12) = 0.25e-6; % preMPF
    x0_c(13) = 0.00025; % Wee1
    x0_c(14) = 0.24999975; % Chk1
    x0_c(15) = 0.05; % Cad3/ATM (transducer)

    % Apoptosis
    x0_a = zeros(32, 1);
    x0_a(1) = 1.0e-4; % Caspase8
    x0_a(2) = 0.0; % Caspase8/Bid
    x0_a(3) = 0.0040; % Apaf-1
    x0_a(4) = 0.0; % Cytochrome C/Apaf-1
    x0_a(5) = 0.0; % Cytochrome C
    x0_a(6) = 0.0040; % Cytochrome C mitochondria
    x0_a(7) = 0.0; % Bax_2
    x0_a(8) = 0.0; % tBid mitochondria
    x0_a(9) = 0.0; % tBid
    x0_a(10) = 0.0; % tBid/Bax
    x0_a(11) = 0.0040; % Bax
    x0_a(12) = 0.0040; % Bcl_2
    x0_a(13) = 0.0; % Caspase3/Bcl_2
    x0_a(14) = 0.0; % Caspase3/Bid
    x0_a(15) = 0.0040; % Bid
    x0_a(16) = 0.0; % Caspase3/IAP
    x0_a(17) = Casp3_0; % Caspase3
    x0_a(18) = 0.0; % Apoptosome
    x0_a(19) = 0.0040; % Procaspase9
    x0_a(20) = 0.0; % Apoptosome/Procaspase9
    x0_a(21) = 0.0; % Apoptosome/Procaspase9_2
    x0_a(22) = 0.0; % Apoptosome/Caspase9_2/Procaspase3
    x0_a(23) = 0.0; % Caspase9/Procaspase3
    x0_a(24) = 0.0040; % Procaspase3, with added reactions
    x0_a(25) = 0.0; % Caspase9/IAP
    x0_a(26) = 0.0; % Apoptosome/Caspase9/IAP
    x0_a(27) = 0.0; % Apoptosome/Caspase9_2/IAP
    x0_a(28) = 0.0040; % IAP
    x0_a(29) = 0.0; % Apoptosome/Caspase9
    x0_a(30) = 0.0; % Caspase9
    x0_a(31) = 0.0; % Apoptosome/Caspase9_2
    x0_a(32) = 0.0; % p21/Procaspase3

    % Fluctiations in the experienced sulforaphane concentration
    uptake_ratio = normrnd(1.0, 0.5); 
    if uptake_ratio < 0 % Probability for this should be small (less than 3%)
        uptake_ratio = 1.0; 
    end
    x0_m(2) = uptake_ratio * x0_m(2);

    k_uptake_mean = 1.207e-7;
    k_uptake = normrnd(k_uptake_mean, k_uptake_mean / 2);
    if k_uptake < 0 
        k_uptake = k_uptake_mean; 
    end

    % Solve system
    proliferation = 1;
    tI = [0 tf];
    x0 = [x0_m; x0_p; x0_c; x0_a];
    [t, x] = ode23s(@(t, x) f_myrosinase_action(t, x, k_uptake), tI, x0);

    [~, I] = max(x(:, 17));
    peak = t(I);
    
    % Apoptosis or cell cycle arrest
    if (max(x(:, 39)) > Casp3_0) || (peak > 1.25 * mpf_peak)
        proliferation = 0;
    end
    
end