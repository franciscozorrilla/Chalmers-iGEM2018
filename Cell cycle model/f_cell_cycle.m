function xdot = f_cell_cycle(t, x, DNAdamage)

k_1p = 0.000025;
k_2p = 1.5;
k_3p = 0.001;
k_4p = 1.0e-8;
k_5p = 0.004;
k_6p = 6.0;
k_7p = 0.005;
k_8p = 10.0;
k_9p = 0.02;
k_10p = 0.00094;
k_ex = 1.0; % Används denna?
Km = 9.5;
k_damp = 0.02;
k_deg = 0.772;

k_1c = 0.2;
k_2c = 1.0;
k_3c = 1.0;
k_4c = 1.0;
k_5c = 0.0001;
k_6c = 0.001;
k_7c = 0.01;
k_8c = 1.0e-5;
k_9c = 0.01;
k_10c = 1.0;
k_11c = 0.01;
k_12c = 1.0;
k_13c = 0.01;
k_14c = 0.01;
k_15c = 1.0;
k_16c = 1.0;
k_17c = 100.0;
k_18c = 1.0;
k_19c = 0.005;
k_20c = 1.0;
k_21c = 1.0;
k_22c = 0.01;
k_23c = 1.0;
k_24c = 0.1;
k_25c = 0.1;
k_26c = 0.01; 
k_27c = 1.0;
k_28c = 0.1;
k_29c = 1.0;
k_30c = 1.0;
k_31c = 2.0e-4; %e-44??

k_1a = 0.006;
k_2a = 0.5;
k_3a = 10.0;
k_4a = 0.1;
k_5a = 0.006;
k_6a = 0.00003;
k_7a = 0.006;
k_8a = 10.0;
k_9a = 0.006;
k_10a = 8.0e-5;
k_11a = 0.006;
k_12a = 3.0e-5;
k_13a = 10.0;
k_14a = 0.006;
k_15a = 10.0;
k_17a = 0.006;
k_18a = 10.0;
k_19a = 0.006;
k_20a = 0.0003;
k_21a = 0.006;
k_22a = 10.0;
k_23a = 0.006;
k_24a = 0.0003;
k_25a = 0.006;
k_26a = 5.0;
k_27a = 0.5;
k_28a = 5.0e4;
k_29a = 0.5;
k_30a = 0.006;
k_31a = 0.5;
k_32a = 5.0;
k_33a = 0.5;
k_34a = 10.0;
k_35a = 10.0;
k_36a = 0.5;
k_37a = 0.5;
k_38a = 5.0;
k_39a = 0.1;
k_40a = 0.0003;
k_41a = 0.1;
k_42a = 0.0003;
k_43a = 0.006;
k_44a = 0.001;
k_45a = 10.0;
k_46a = 0.5;
k_47a = 0.0035;
k_48a = 5.0;
k_49a = 5.0;
k_50a = 0.0035;
k_51a = 5.0;
k_52a = 0.0035;
k_53a = 3.0e-5;
k_54a = 0.006;
k_55a = 0.0035;
k_56a = 5.0;
k_57a = 10.0;
k_58a = 0.5;
k_59a = 0.1;
k_60a = 10.0;
k_61a = 0.5;
k_62a = 0.1;
k_63a = 0.006;
k_64a = 1.0;
k_65a = 0.1;
k_66a = 1.0;
k_67a = 0.5;
k_68a = 10.0;
k_69a = 0.006;

p53_threshold = 0.4;
p = 4.0;

nr_species = 51;
xdot = zeros(nr_species, 1);

% 1 = p53, 2 = Mdm2, 3 = II, 4 = Chk1, 5 = Chk1p, 6 = transducer, 7 = iMPF,
% 8 = aMPF, 9 = Wee1, 10 = aCdc25, 11 = aCdc25p, 12 = p21, 13 = p21:aMPF,
% 14 = procaspase3, 15 = caspase3, 16 = p21:proscapase3, 17 = iCdc25, 18 =
% iCdc25p, 19 = 14_3_3, 20 = aCdc25p:14_3_3, 21 = iCdc25:14-3-3, 22 =
% Wee1p, 23 = Apaf_1, 24 = Cyt_c, 25 = Cyt_c:Apaf_1, 26 = caspase8,
% 27 = apoptosome, 28 = proscapase9, 29 = apoptosome:procaspase9, 30 =
% apoptosome:caspase9, 31 = caspase9, 32 = apoptosome:proscapase9_2, 33 =
% apoptosome:caspase9_2, 34 = IAP, 35 = Cyt_cmito, 36 =
% apoptosome:proscapase9_2:procaspase3, 37 = apoptosome:caspase9:IAP, 38 =
% caspase9:IAP, 39 = caspase9:procaspase3, 40 = apoptosome:caspase9_2:IAP,
% 41 = caspase3:IAP, 42 = Bid, 43 = caspase3:Bid, 44 = Bcl_2, 45 =
% caspase3:Bcl_2, 46 = Bax_2, 47 = caspase8:Bid, 48 = Bax, 49 = tBid, 50
% = tBid_mito, 51 = tBid:Bax

% Vad hände med 21:an? Finns den ens? Är det samma som 20?

% p53 oscillation system
Deg0 = 0.0556;
Degt = Deg0 - k_deg * (DNAdamage * exp(-k_4p * t) - DNAdamage * exp(-k_damp * DNAdamage * t));
xdot(1) = k_1p + k_2p * (DNAdamage * exp(-k_4p * t)) - k_3p * x(1) - Degt * x(1) * x(2);
xdot(2) = k_10p - k_9p * x(2) + k_8p * ((x(3)^9.0)/(Km^9.0 + x(3)^9.0));
xdot(3) = k_6p * x(1) * ((DNAdamage * exp(-k_4p * t))/(1.0 + k_5p * x(1) * x(2))) - k_7p * x(3);

% G2/M phase cell cycle arrest
xdot(4) = k_4c * x(5) - k_3c * x(4) * x(6);
xdot(5) = k_3c * x(4) * x(6) - k_4c * x(5);
xdot(6) = k_1c * (DNAdamage * exp(-k_4p * t)) - k_2c * x(6);
xdot(7) = (k_19c/(1.0 + x(1))) + k_21c * x(8) * x(9) - k_20c * (x(10) + x(11)) * x(7);
xdot(8) = k_20c * (x(10) + x(11)) * x(7) + k_23c * x(13) - k_21c * x(8) * x(9) - k_24c * x(8) * x(12) - k_22c * x(8) * x(8);
xdot(9) = k_31c + k_29c * x(22) - k_28c * x(8) * x(9);
xdot(10) = k_10c * x(8) * x(17) + k_7c * x(11) - k_5c * x(10) - k_6c * x(5) * x(10) - k_9c * x(10); % Note, table said iCdc25 but there is already one for that
xdot(11) = k_6c * x(5) * x(19) + k_12c * x(8) * x(18) - k_7c * x(11) - k_13c * x(11);
xdot(12) = -k_65a * x(12) * x(14) - k_64a * x(12) * x(15) + k_66a * x(16) + k_25c * x(1) + k_26c + k_23c * x(13) - k_27c * x(12) - k_24c * x(8) * x(12);
xdot(13) = k_24c * x(8) * x(12) - k_23c * x(13);
xdot(17) = k_9c * x(10) + k_8c - k_10c * x(8) * x(17) - k_11c * x(5) * x(17);
xdot(18) = k_13c * x(11) + k_11c * x(5) * x(17) - k_12c * x(8) * x(18) - k_17c * x(18) * x(19);
xdot(19) = k_14c * x(1) + k_15c - k_17c * x(18) * x(19) - k_16c * x(19);
xdot(20) = k_17c * x(18) * x(19) - k_18c * x(20); % x(21) exhanged for x(20)
xdot(21) = 0.0; % Does not exist?
xdot(22) = k_28c * x(8) * x(9) - (k_29c + k_30c) * x(22);

% Apoptosis induction system
xdot(14) = -k_65a * x(12) * x(14) + k_66a * x(16) - k_45a * x(31) * x(14) + k_46a * x(39) - k_57a * x(36) + k_58a * x(36) + k_42a - k_43a * x(14);
xdot(15) = k_44a * x(39) + k_41a * x(36) - k_56a * x(15) * x(34) + k_55a * x(41) - k_68a * x(15) * x(42) + (k_59a + k_67a) * x(43) - k_61a * x(15) * x(44) + (k_60a + k_62a) * x(45) - k_63a * x(15); % Note: used k_55a*caspase3:IAP
xdot(16) = k_65a * x(12) * x(14) - k_66a * x(16);
xdot(23) = -k_26a * x(24) * x(23) + k_27a * x(25) + k_24a - k_25a * x(23);
xdot(24) = k_22a * x(46) * x(35) - k_26a * x(24) * x(23) + k_27a * x(25) - k_23a * x(24);
xdot(25) = k_26a * x(24) * x(23) - k_27a * x(25) - 7.0 * (k_28a * x(25)^p + k_29a * x(27));
xdot(26) = -k_3a * x(26) * x(42) + k_2a * x(47) + k_4a * x(47) - k_1a * x(26);
xdot(27) = k_28a * x(25)^p - k_29a * x(27) - k_34a * x(27) * x(28) + k_33a * x(29) + k_32a * x(30) - k_31a * x(27) * x(31);
xdot(28) = -k_34a * x(27) * x(28) + k_33a * x(29) - k_35a * x(29) * x(28) + k_36a * x(32) + k_40a - k_30a * x(28);
xdot(29) = k_34a * x(27) * x(28) - k_33a * x(29) - k_35a * x(29) * x(28) + k_36a * x(32);
xdot(30) = k_38a * x(33) - k_37a * x(30) * x(31) - k_32a * x(30) + k_31a * x(27) * x(31) - k_48a * x(30) * x(34) + k_47a * x(37);
xdot(31) = k_38a * x(33) - k_37a * x(30) * x(31) + k_32a * x(30) - k_31a * x(27) * x(31) - k_51a * x(31) * x(34) + k_52a * x(38) - k_45a * x(31) * x(14) + k_46a * x(39) + k_44a * x(39) - k_69a * x(31);
xdot(32) = k_35a * x(29) * x(28) - k_36a * x(32) - k_39a * x(32); 
xdot(33) = k_39a * x(32) - k_38a * x(33) + k_37a * x(30) * x(31) - k_49a * x(33) * x(34) + k_50a * x(40) - k_57a * x(33) * x(14) + k_58a * x(36) + k_41a * x(36); 
xdot(34) = -k_51a * x(31) * x(34) + k_52a * x(38) - k_48a * x(30) * x(34) + k_47a * x(37) - k_49a * x(33) * x(34) + k_50a * x(40) - k_56a * x(15) * x(34) + k_55a * x(41) + k_53a - k_54a * x(34);
xdot(35) = -k_22a * x(46) * x(35) + k_20a - k_21a * x(35);
xdot(36) = k_57a * x(33) * x(14) - k_58a * x(36) - k_41a * x(36);
xdot(37) = k_48a * x(30) * x(34) - k_47a * x(37);
xdot(38) = k_51a * x(31) * x(34) - k_52a * x(38);
xdot(39) = k_45a * x(31) * x(14) - k_46a * x(39) - k_44a * x(39);
xdot(40) = k_49a * x(33) * x(34) - k_50a * x(40);
xdot(41) = k_56a * x(15) * x(34) - k_55a * x(41);
xdot(42) = -k_3a * x(26) * x(42) + k_2a * x(47) - k_68a * x(15) * x(42) + k_67a * x(43) + k_6a - k_5a * x(42);
xdot(43) = k_68a * x(15) * x(42) - k_67a * x(43) - k_59a * x(43);
xdot(44) = -k_61a * x(15) * x(44) + k_60a * x(45) - k_13a * x(44) * x(48) + ((k_10a * p53_threshold^4.0)/(x(1)^4.0 + p53_threshold^4.0)) - k_11a * x(44);
xdot(45) = k_61a * x(15) * x(44) - k_60a * x(45) - k_62a * x(45);
xdot(46) = k_18a * x(51) * x(48) - k_19a * x(46);
xdot(47) = k_3a * x(26) * x(42) - k_2a * x(47) - k_4a * x(47);
xdot(48) = -k_15a * x(50) * x(48) - k_18a * x(51) * x(48) - k_13a * x(44) * x(48) + ((k_12a * (1.0 + x(1)^4.0))/(x(1)^4.0 + p53_threshold^4.0)) - k_14a * x(48);
xdot(49) = k_4a * x(47) + k_59a * x(43) - k_8a * x(49) + k_18a * x(51) * x(48) - k_7a * x(49);
xdot(50) = k_8a * x(49) - k_15a * x(50) * x(48) - k_9a * x(50);
xdot(51) = k_15a * x(50) * x(48) - k_18a * x(51) * x(48) - k_17a * x(51);

end