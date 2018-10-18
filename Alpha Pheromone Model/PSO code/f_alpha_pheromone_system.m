function xdot = f_alpha_pheromone_system(t, x, v48_k48, v49_k49, v50_k50, v51_k51, pop_size)
% Alpha pheromone pathway with production of anti-cancer protein
% Constants and reactions from Kofahl and Klipp (2004) unless otherwise
% stated

    v1_k1 = 0.03;
    v2_k2 = 0.0012;
    v3_k3 = 0.6;
    v4_k4 = 0.24;
    v5_k5 = 0.024;
    v6_k6 = 0.0036;
    v7_k7 = 0.24;
    v8_k8 = 0.33; 
    v9_k9 = 2000.0;
    v10_k10 = 0.1;
    v11_k11 = 5.0;
    v12_k12 = 1.0;
    v13_k13 = 3.0;
    v14_k14 = 1.0;
    v15_k15 = 3.0;
    v16_k16 = 3.0;
    v17_k17 = 100.0;
    v18_k18 = 5.0;
    v19_k19 = 1.0;
    v20_k20 = 10.0;
    v21_k21 = 5.0;
    v22_k22 = 47.0;
    v23_k23 = 5.0;
    v24_k24 = 345.0;
    v25_k25 = 5.0;
    v26_k26 = 50.0;
    v27_k27 = 5.0;
    v28_k28 = 140.0;
    v29_k29 = 10.0;
    v30_k30 = 1.0;
    v31_k31 = 250.0;
    v32_k32 = 5.0;
    v33_k33 = 50.0;
    v34_k34 = 18.0;
    v35_k35 = 10.0;
    v36_k36 = 0.1;
    v37_k37 = 0.1;
    v38_k38 = 0.01;
    v39_k39 = 18.0;
    v40_k40 = 1.0;
    v41_k41 = 0.002; 
    v42_k42 = 0.1;
    v43_k43 = 0.01;
    v44_k44 = 0.01;
    v45_k45 = 0.1;
    v46_k46 = 200.0;    
    v47_k47 = 1.0;
    
    v1 = x(1) * x(29) * v1_k1;
    v2 = x(2) * x(1) * v2_k2;
    v3 = x(3) * v3_k3;
    v4 = x(3) * v4_k4;
    v5 = x(2) * v5_k5;
    v6 = x(3) * x(4) * v6_k6;
    v7 = x(5) * v7_k7;
    v8 = x(5) * x(36) * v8_k8;
    v9 = x(7) * x(6) * v9_k9;
    v10 = x(6) * x(8) * v10_k10;
    v11 = x(9) * v11_k11;
    v12 = x(10) * x(11) * v12_k12;
    v13 = x(12) * v13_k13;
    v14 = x(13) * x(14) * v14_k14;
    v15 = x(15) * v15_k15;
    v16 = x(12) * x(15) * v16_k16;
    v17 = x(8) * v17_k17;
    v18 = x(9) * x(16) * v18_k18;
    v19 = x(17) * v19_k19;
    v20 = x(17) * v20_k20;
    v21 = x(17) * v21_k21;
    v22 = x(18) * v22_k22;
    v23 = x(18) * v23_k23;
    v24 = x(19) * v24_k24;
    v25 = x(19) * v25_k25;
    v26 = x(20) * v26_k26;
    v27 = x(20) * v27_k27;
    v28 = x(21) * v28_k28;
    v29 = x(22) * x(14) * v29_k29;
    v30 = x(24) * v30_k30;
    v31 = x(24) * v31_k31;
    v32 = x(22) * v32_k32;
    v33 = x(23) * v33_k33;
    v34 = x(25) * x(23) * v34_k34;
    v35 = x(26) * v35_k35;
    v36 = x(26) * x(27) * v36_k36;
    v37 = x(28) * v37_k37;
    v38 = x(28) * v38_k38;
    v39 = x(30) * x(23)^2/(100^2 + x(23)^2) * v39_k39; 
    v40 = x(31) * v40_k40;
    v41 = x(30) * x(35) * v41_k41;
    v42 = x(6) * x(31) * v42_k42;
    v43 = x(33) * v43_k43;
    v44 = x(34) * v44_k44;
    v45 = x(31) * x(35) * v45_k45;
    v46 = x(23)^2 / (4^2+x(23)^2) * v46_k46;
    v47 = x(36) * v47_k47;
    
    % The following reactions have been added
    % Production of mRNA
    v48 = x(26) * v48_k48; 
    % Degradation of mRNA
    v49 = x(37) * v49_k49; 
    % Production of protein
    v50 = x(37) * v50_k50;
    % Degradation of protein
    v51 = x(38) * v51_k51;

    % ODEs
    xdot = zeros(38, 1);
    % Alpha-factor
    xdot(1) = - pop_size * v1;
    % Ste2
    xdot(2) = - v2 + v3 - v5;
    % Ste2a
    xdot(3) = v2 - v3 - v4;
    % Gabc
    xdot(4) = - v6 + v9;
    % GaGTP 
    xdot(5) = v6 - v7 - v8;
    % Gbc 
    xdot(6) = v6 - v9 - v10 + v11 + v21 + v23 + v25 + v27 + v32 - v42 + v43;
    % GaGDP 
    xdot(7) = v7 + v8 - v9;
    % Complex C 
    xdot(8) = - v10 + v11 + v16 - v17;
    % Complex D 
    xdot(9) = v10 - v11 - v18 + v19;
    % Ste5
    xdot(10) = - v12 + v13 + v17 + v21 + v23 + v25 + v27 + v32;
    % Ste11
    xdot(11) = - v12 + v13 + v17 + v21 + v23 + v25 + v27 + v32;
    % Complex A 
    xdot(12) = v12 - v13 - v16;
    % Ste7 
    xdot(13) = - v14 + v15 + v17 + v21 + v23 + v25 + v27 + v32;
    % Fus3
    xdot(14) = - v14 + v15 + v17 + v21 + v23 + v25 + v27 - v29 + v30 + v33;
    % Complex B
    xdot(15) = v14 - v15 - v16;
    % Ste20 
    xdot(16) = - v18 + v19 + v21 + v23 + v25 + v27 + v32;
    % Complex E
    xdot(17) = v18 - v19 - v20 - v21;
    % Complex F
    xdot(18) = v20 - v22 - v23;
    % Complex G
    xdot(19) = v22 - v24 - v25;
    % Complex H
    xdot(20) = v24 - v26 - v27;
    % Complex I
    xdot(21) = v26 - v28 + v31;
    % Complex L
    xdot(22) = v28 - v29 + v30 - v32;
    % Fus3PP 
    xdot(23) = v28 - v33 - v34 + v35;
    % Complex K 
    xdot(24) = v29 - v30 - v31;
    % Ste12 
    xdot(25) = - v34 + v35;
    % Ste12 active 
    xdot(26) = v34 - v35;
    % Bar1
    xdot(27) = - v36 + v37;
    % Bar1 active, 
    xdot(28) = v36 - v37 - v38;
    % Bar1 active ex 
    xdot(29) = v38;
    % Far1, temporarily removed
    xdot(30) = 0.0;
    % Far1PP 
    xdot(31) = v39 - v40 - v42 + v43 + v44 - v45;
    % Far1 ubiquitin 
    xdot(32) = v41;
    % Complex M 
    xdot(33) = v42 - v43;
    % Complex N 
    xdot(34) = - v44 + v45;
    % Cdc28 
    xdot(35) = v44 - v45;
    % Sst2
    xdot(36) = v46 - v47;
    % The following species has been added
    % mRNA
    xdot(37) = v48 - v49;
    % Protein
    xdot(38) = v50 - v51;
    
end