function [proliferation, peak, err] = p28_action(p28_0, tf, mpf_peak, source_p53, d_p53, f1, d_Mdm2, f2, d_cop1, k_deg_cop1, mean_uptake_rate, deg_basal, d_p28, k_67_a)
% Initial concentrations from Hamada et al. unless otherwise stated.
% Concentrations in microM
    Casp3_0 = 1.0e-5;
    err = 0;

    x0_p = zeros(7, 1);

    k_cop1_1 = 0.1; % Moscetti et al. 
    k_cop1_2 = 1.1e-3; % Moscetti et al. 

    x0_p(1) = (d_cop1 * d_Mdm2 * (k_cop1_2 + k_deg_cop1) * source_p53)/...
    (d_Mdm2 * f2 * k_cop1_1 * k_deg_cop1 + d_cop1 * (d_Mdm2 * d_p53 + deg_basal * f1) *  (k_cop1_2 + k_deg_cop1));
    x0_p(2) = f1 / d_Mdm2; % Mdm2
    x0_p(3) = f2 / d_cop1; % COP1
    x0_p(4) = (d_Mdm2 * f2 * k_cop1_1 * source_p53)/...
    (d_Mdm2 * f2 * k_cop1_1 * k_deg_cop1 + d_cop1 *(d_Mdm2 * d_p53 + deg_basal * f1) * (k_cop1_2 + k_deg_cop1));
    x0_p(5) = 0.0; % p28 intracellular

    if (max(x0_p) > 5.0) || (min(x0_p) < 0) 
        x0_p(:) = 0;
        err = 1;
    end
    x0_p(6) = 0.0; % p53/p28
    x0_p(7) = p28_0; % p28 extracellular

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

    % Solve system
    x0 = [x0_p; x0_c; x0_a];

    % Fluctiations in the experienced p28 concentration
    uptake_ratio = normrnd(1.0, 0.5); 
    if uptake_ratio < 0 % Probability for this should be small (less than 3%)
        uptake_ratio = 1.0; 
    end
    x0(7) = uptake_ratio * x0(7);

    % Fluctuations in uptake rate
    uptake_rate = normrnd(mean_uptake_rate, mean_uptake_rate / 2);
    if uptake_rate < 0 
        uptake_rate = mean_uptake_rate; 
    end
    

    tI = [0 tf];
    proliferation = 1;
    [t, x] = ode23s(@(t, x) f_p28_action(t, x, source_p53, d_p53, f1, d_Mdm2, f2, d_cop1, k_deg_cop1, uptake_rate, deg_basal, d_p28, k_67_a), tI, x0);
    
    [msg, ID] = lastwarn;
    
    if strcmp(ID, 'MATLAB:illConditionedMatrix')
        err = 1;
        lastwarn('reset');
    end
    
    [~, I] = max(x(:, 17));
    peak = t(I);
    
    % Apoptosis or cell cycle arrest
    if (max(x(:, 39)) > Casp3_0) || (peak > 1.25 * mpf_peak)
        proliferation = 0;
    end
    
end