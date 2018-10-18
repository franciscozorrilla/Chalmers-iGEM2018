function xdot = f_myrosinase_action(t, x, uptake_rate, kbax_m, kcytc_m, kp21_m)
% Values of constants from Hamuda et al. unless otherwise stated.
% Fluxes assumed to be in microM/s

    x_m = x(1:3);
    x_p = x(4:7);
    x_c = x(8:22);
    x_a = x(23:54);

    % Myrosinase
    km_myr = 160; % Hochkoeppler et al. [nM * s-1]
    kcat_myr = 1.19e2; % Hochkoeppler et al. [s-1]

    % p53 system
    source_p53 = 7.319e-3; 
    d_p53 = 3.419e-6;
    f1 = 2.447e-8;
    d_Mdm2 = 7.336e-4;
    deg_basal = 0.1144; 
    k_cop1_1 = 1.0e-1; 
    k_cop1_2 = 1.1e-3; 
    k_deg_cop1 = 2.1436; 
    f2 = 6.819e-4;
    d_cop1 = 1.306e-3;

    % Cell cycle arrest
    k1_c = 1.0;
    k1_0_c = 1.0;
    k2_c = 0.001;
    k2_0_c = 0.01;
    k2_1_c = 0.01;
    k_ctak1 = 0.0;
    k_ctak1_0 = 0.0;
    k3_c = 100.0;
    k_3_c = 1.0;
    k4_c = 0.0; 
    k5_c = 1.0; 
    k6_c = 0.01;
    k7_c = 1.0;
    k_Plk1 = 0.0;
    k7_0_c = 0.01;
    k7_1_c = 1.0e-4;
    k8_c = 0.1;
    k8_0_c = 0.0;
    k8_1_c = 1.0;
    k9_c = 1.0;
    k9_0_c = 0.0;
    k9_1_c = 1.0;
    k10_c = 0.0;
    k11_c = 0.0;
    k12_c = 1.0;
    k13_c = 1.0;
    k14_c = 5.0e-3;
    k14_0_c = 1.0;
    k15_c = 0.01;
    k16_c = 2.0e-4;
    k17_c = 1.0;
    k18_c = 1.0;
    k_Plk1_0 = 0.0;
    k18_0_c = 0.01;
    k20_c = 0.1;
    k21_c = 0.01;
    k22_c = 1.0;
    k23_c = 0.1;
    k_23_c = 1.0;
    kex = 0.0;
    kin = 1.0e-5;

    % Apoptosis
    k_10_a = 10.0;
    k_00_a = 0.5;
    k_f0_a = 0.1;
    k_11_a = 5.0;
    k_01_a = 0.5;
    k_11b_a = 5.0e4;
    k_01b_a = 0.5;
    k_12_a = 10.0;
    k_02_a = 0.5;
    k_13_a = 10.0;
    k_03_a = 0.5;
    k_f3_a = 0.1;
    k_14_a = 5.0;
    k_04_a = 0.5;
    k_14b_a = 5.0;
    k_04b_a = 0.5;
    k_15_a = 5.0;
    k_05_a = 0.0035;
    k_15b_a = 5.0;
    k_05b_a = 0.0035;
    k_15c_a = 5.0;
    k_05c_a = 0.0035;
    k_16_a = 10.0;
    k_06_a = 0.5;
    k_16b_a = 10.0;
    k_06b_a = 0.5;
    k_f6_a = 0.0010;
    k_f6b_a = 0.1;
    k_17_a = 5.0;
    k_07_a = 0.0035;
    k_18_a = 10.0;
    k_08_a = 0.5;
    k_f8_a = 0.1;
    k_19_a = 10.0;
    k_09_a = 0.5;
    k_f9_a = 0.1;
    k11_a = 10.0;
    k12_a = 10.0;
    k12b_a = 10.0;
    k13_a = 10.0;
    k14_a = 10.0;
    p53_thresh = 0.40; 
    u = 0.0060;
    u_Bax = 0.0060;
    u_Bcl_2 = 0.0060;
    P_Apaf_1 = 3.0e-4;
    P_IAP = 3.0e-5;
    P_Pro3 = 3.0e-4;
    P_Pro9 = 3.0e-4;
    P_Bid = 3.0e-5;
    P_oBcl_2 = 8.0e-5;
    P_oBax = 3.0e-5;
    P_Cytc_mito = 3.0e-4;
    p = 4.0;

    % Added constants
    k_64_a = 1.0;
    k_65_a = 0.1;
    k_66_a = 1.0;
    k_67_a = 1.344e-6;

    % Myrosinase
    gluc_conc = 650; % Approximated glucosinolate concentration
    M1 = kcat_myr * x_m(1) * gluc_conc / (km_myr + gluc_conc);
    M2 = uptake_rate * x_m(2); 
    M_Bax = kbax_m * x_m(3);
    M_Cytc = kcytc_m * x_m(3);
    M_p21 = kp21_m * x_m(3);

    xdot_m = zeros(3, 1);
    % Myrosinase
    xdot_m(1) = 0.0; % Assumption: constant enzyme concentration
    % Sulforaphane extracellular
    xdot_m(2) = M1 - M2; 
    % Sulforaphane intracellular
    xdot_m(3) = M2;

    % p53 system
    P1 = source_p53 - d_p53 * x_p(1) - deg_basal * x_p(1) * x_p(2); 
    P2 = f1 - d_Mdm2 *  x_p(2);
    P3 = k_cop1_1 * x_p(1) * x_p(3) - k_cop1_2 * x_p(4);
    P4 = f2 - d_cop1 * x_p(3);
    P5 = k_deg_cop1 * x_p(4);

    xdot_p = zeros(4, 1);
    % p53
    xdot_p(1) = P1 - P3;
    % Mdm2
    xdot_p(2) = P2;
    % COP1
    xdot_p(3) = P4 - P3 + P5; 
    % p53/COP1
    xdot_p(4) = P3 - P5;

    % Cell cycle arrest
    C1 = k1_c * x_c(14)* x_c(15);
    C1_ = k1_0_c * x_c(5);
    C2 = k2_c * x_c(5) * x_c(1) + k_ctak1 * x_c(1);
    C2_ = k2_0_c * x_c(2);
    C2_0 = k2_1_c * x_c(5) * x_c(3) + k_ctak1_0 * x_c(3);
    C3 = k_3_c * x_c(7);
    C3_ = k3_c * x_c(6) * x_c(4);
    C4 = k4_c;
    C5 = k5_c * x_c(15);
    C6 = k6_c * x_p(1);
    C7 = k7_c * x_c(10) * x_c(3) + k_Plk1 * x_c(3);
    C7_= k7_0_c * x_c(1);    
    C7_1 = k7_1_c * x_c(1);
    C8 = k8_c * x_c(10) * x_c(13) + k8_0_c * x_c(13);
    C8_ = k8_1_c * x_c(8);
    C9 = k9_c * (x_c(1) + x_c(2)) * x_c(12) + k9_0_c * x_c(12);
    C9_ = k9_1_c * x_c(10) * x_c(13);
    C12 = k12_c * x_c(6);
    C13 = k13_c;
    C14 = k14_c / (1 + k14_0_c * x_p(1));
    C15 = k15_c * x_c(10)^2;
    C16 = k16_c;
    C17 = k17_c * x_c(8);
    C18 = k18_c * x_c(10) * x_c(4) + k_Plk1_0 * x_c(4);
    C18_ = k18_0_c *x_c(2);
    C20 = k20_c * x_p(1);
    C21 = k21_c;
    C22 = k22_c * x_c(9);
    C23 = k23_c * x_c(9) * x_c(10);
    C23_ = k_23_c * x_c(11);
    Cex = kex * x_c(7);
    Cin = kin;

    % Added reactions
    R_A1 = k_64_a * x_c(9) * x_a(17);
    R_A2 = k_65_a * x_c(9) * x_a(24);
    R_A3 = k_66_a * x_a(32);
    R_A4 = k_67_a * ((x_p(1)^4)/((x_p(1)^4)  +  p53_thresh^4));

    xdot_c = zeros(15, 1);
    % Cdc25 active
    xdot_c(1) = - C2 +  C2_ - C7_ +  C7 - C7_1;
    % Cdc25Ps216 active
    xdot_c(2) = C2 - C2_ +  C18 - C18_;
    % Cdc25 inactive
    xdot_c(3) = C7_ - C7 - C2_0 +  C3 +  Cin;
    % Cdc25Ps216 inactive
    xdot_c(4) = - C18 +  C18_ +  C2_0 - C3_;
    % Chk1P
    xdot_c(5) = C1 - C1_;
    % protein 14-3-3
    xdot_c(6) = C3 - C12 +  C13 +  C6 - C3_;
    % Cdc25Ps216_14-3-3 inactive
    xdot_c(7) = - C3 +  C3_ - Cex;
    % Wee1 phosphorylated
    xdot_c(8) = C8 - C8_ - C17;	
    % p21, with added reactions
    xdot_c(9) = C20 +  C21 - C22 - C23 +  C23_ - R_A2 - R_A1 + R_A3 + M_p21;	
    % MPF
    xdot_c(10) = - C23 +  C23_ +  C9 - C9_ - C15;	
    % p21/MPF
    xdot_c(11) = C23 - C23_;	
    % preMPF
    xdot_c(12) = - C9 +  C9_ +  C14;	
    % Wee1
    xdot_c(13) = C16 - C8 +  C8_;	
    % Chk1
    xdot_c(14) = - C1 +  C1_;	
    % Cad3/ATM
    xdot_c(15) = C4 - C5;	

    % Apoptosis
    A_f0 = k_f0_a * x_a(2);
    A_Casp8 = - u * x_a(1);
    A_0 = k_10_a * x_a(1) * x_a(15) - k_00_a * x_a(2);
    A_Apaf_1 = P_Apaf_1 - u * x_a(3);
    A_1 = k_11_a * x_a(5) * x_a(3) - k_01_a * x_a(4);
    A_1b = k_11b_a * x_a(4)^p - k_01b_a * x_a(18);
    A_Cytc = - u * x_a(5);
    A_Cytc_mito = P_Cytc_mito - u * x_a(6);
    A_14 = k14_a * x_a(7) * x_a(6);
    A_Bax_2 = - u * x_a(7);
    A_tBid_mito = - u * x_a(8);
    A_12a = k12_a * x_a(8) * x_a(11);
    A_11 = k11_a * x_a(9);
    A_tBid = - u * x_a(9);
    A_tBidBax = - u * x_a(10);
    A_12b = k12b_a * x_a(10) * x_a(11);
    P_Bax = P_oBax * (1  +  x_p(1)^4/(x_p(1)^4  +  p53_thresh^4)); 
    A_Bax = P_Bax - u_Bax * x_a(11);
    A_13 = k13_a * x_a(12) * x_a(11);
    P_Bcl_2 = P_oBcl_2 * p53_thresh^4/(x_p(1)^4  +  p53_thresh^4);
    A_9 = k_19_a * x_a(17) * x_a(12) - k_09_a * x_a(13);
    A_Bcl_2 = P_Bcl_2 - u_Bcl_2 * x_a(12);
    A_f9 = k_f9_a * x_a(13);
    A_8 = k_18_a * x_a(17) * x_a(15) - k_08_a * x_a(14);
    A_f8 = k_f8_a * x_a(14);
    A_Bid = P_Bid - u * x_a(15);
    A_Casp3 = - u * x_a(17);
    A_7 = k_17_a * x_a(17) * x_a(28) - k_07_a * x_a(16);
    A_Pro9 = P_Pro9 - u * x_a(19);
    A_2 = k_12_a * x_a(18) * x_a(19) - k_02_a * x_a(20);
    A_3 = k_13_a * x_a(20) * x_a(19) - k_03_a * x_a(21);
    A_f3 = k_f3_a * x_a(21);
    A_6b = k_16b_a * x_a(31) * x_a(24) - k_06b_a * x_a(22);
    A_f6b = k_f6b_a * x_a(22);
    A_f6 = k_f6_a * x_a(23);
    A_6 = k_16_a * x_a(30) * x_a(24) - k_06_a * x_a(23);
    A_Pro3 = P_Pro3 - u * x_a(24);
    A_5 = k_15_a * x_a(30) * x_a(28) - k_05_a * x_a(25);
    A_5b = k_15b_a * x_a(29) * x_a(28) - k_05b_a * x_a(26);
    A_5c = k_15c_a * x_a(31) * x_a(28) - k_05c_a * x_a(27);
    A_IAP = P_IAP - u * x_a(28);
    A_4 = k_14_a * x_a(31) - k_04_a * x_a(29) * x_a(30);
    A_4b = k_14b_a * x_a(29) - k_04b_a * x_a(18) * x_a(30);
    A_Casp9 = - u * x_a(30);

    xdot_a = zeros(32, 1);
    % Caspase8
    xdot_a(1) =  - A_0 + A_f0 + A_Casp8 + R_A4;
    % Caspase8/Bid
    xdot_a(2) = A_0 - A_f0;
    % Apaf-1
    xdot_a(3) = - A_1 + A_Apaf_1;
    % CytochromeC/Apaf-1
    xdot_a(4) = A_1 - 7 * A_1b;
    % CytochromeC
    xdot_a(5) = A_14 - A_1 + A_Cytc + M_Cytc;
    % CytochromeC mitochondria
    xdot_a(6) = - A_14 + A_Cytc_mito - M_Cytc;
    % Bax_2
    xdot_a(7) = A_12b + A_Bax_2;
    % tBid mitochondria
    xdot_a(8) = A_11 - A_12a + A_tBid_mito;
    % tBid
    xdot_a(9) = A_f0 + A_f8 - A_11 + A_12b + A_tBid;
    % tBid/Bax
    xdot_a(10) = A_12a - A_12b + A_tBidBax;
    % Bax
    xdot_a(11) = - A_12a - A_12b - A_13 + A_Bax + M_Bax;
    % Bcl-2
    xdot_a(12) = - A_9 - A_13 + A_Bcl_2;
    % Caspase3/Bcl-2
    xdot_a(13) = A_9 - A_f9;
    % Caspase3/Bid
    xdot_a(14) = A_8 - A_f8;
    % Bid
    xdot_a(15) =  - A_0 - A_8 + A_Bid;
    % Caspase3/IAP
    xdot_a(16) = A_7;
    % Caspase3
    xdot_a(17) = A_f6 + A_f6b - A_7 - A_8 + A_f8 - A_9 + A_f9 + A_Casp3;
    % Apoptosome
    xdot_a(18) = A_1b - A_2 + A_4b;
    % Procaspase9
    xdot_a(19) =  - A_2 - A_3 + A_Pro9;
    % Apoptosome/Procaspase9
    xdot_a(20) = A_2 - A_3;
    % Apoptosome/Procaspase9_2
    xdot_a(21) = A_3 - A_f3;
    % Apoptosome/Caspase9_2/Procaspase3
    xdot_a(22) = A_6b - A_f6b;
    % Caspase9/Procaspase3
    xdot_a(23) = A_6 - A_f6;
    % Procaspase3, with added reactions
    xdot_a(24) =  - A_6 - A_6b + A_Pro3 - R_A2 + R_A3;
    % Caspase9/IAP
    xdot_a(25) = A_5;
    % Apoptosome/Caspase9/IAP
    xdot_a(26) = A_5b;
    % Apoptosome/Caspase9_2/IAP
    xdot_a(27) = A_5c;
    % IAP
    xdot_a(28) =  - A_5 - A_5b - A_5c - A_7 + A_IAP;
    % Apoptosome/Caspase9
    xdot_a(29) = A_4 - A_4b - A_5b;
    % Caspase9
    xdot_a(30) = A_4 + A_4b - A_5 - A_6 + A_f6 + A_Casp9;
    % Apoptosome/Capsase9_2
    xdot_a(31) = A_f3 - A_4 - A_5c - A_6b + A_f6b;
    % Species added:
    xdot_a(32) = R_A2 - R_A3;

    xdot = [xdot_m; xdot_p; xdot_c; xdot_a];

end