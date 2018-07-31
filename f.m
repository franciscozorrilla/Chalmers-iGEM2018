function xdot = f(t,x,reaction_v48_k48,nr_species,pop_size)
    % TO DO: COMPARE ALL REACTIONS TO ARTICLE   
    % Compartment: id = Extracellular, name = Extracellular, constant
    compartment_extracellular = 1.0; % Change?

    % Compartment: id = compartment, name = cell, constant
    compartment_intracellular = 1.0;
    
    % The compartment-things are only 1.0 so won't do any difference

    % Reaction: id = v1, name = alpha degradation
    % Local Parameter:   id =  k1, name = k1
    reaction_v1_k1 = 0.03;
    reaction_v1 = compartment_extracellular*x(1)*x(29)*reaction_v1_k1;

    % Reaction: id = v2
    % Local Parameter:   id =  k2, name = k2
    reaction_v2_k2 = 0.0012;
    reaction_v2 = compartment_intracellular*x(2)*x(1)*reaction_v2_k2;

    % Reaction: id = v3
    % Local Parameter:   id =  k3, name = k3
    reaction_v3_k3 = 0.6;
    reaction_v3 = compartment_intracellular*x(3)*reaction_v3_k3;

    % Reaction: id = v4
    % Local Parameter:   id =  k4, name = k4
    reaction_v4_k4 = 0.24;
    reaction_v4 = compartment_intracellular*x(3)*reaction_v4_k4;

    % Reaction: id = v5
    % Local Parameter:   id =  k5, name = k5
    reaction_v5_k5 = 0.024;
    reaction_v5 = compartment_intracellular*x(2)*reaction_v5_k5;

    % Reaction: id = v6
    % Local Parameter:   id =  k6, name = k6
    reaction_v6_k6 = 0.0036;
    reaction_v6 = compartment_intracellular*x(3)*x(4)*reaction_v6_k6;

    % Reaction: id = v7
    % Local Parameter:   id =  k7, name = k7
    reaction_v7_k7 = 0.24;
    reaction_v7 = compartment_intracellular*x(5)*reaction_v7_k7;

    % Reaction: id = v8
    % Local Parameter:   id =  k8, name = k8
    reaction_v8_k8 = 0.033;
    reaction_v8 = compartment_intracellular*x(5)*x(36)*reaction_v8_k8;

    % Reaction: id = v9
    % Local Parameter:   id =  k9, name = k9
    reaction_v9_k9 = 2000.0;
    reaction_v9 = compartment_intracellular*x(7)*x(6)*reaction_v9_k9;

    % Reaction: id = v10
    % Local Parameter:   id =  k10, name = k10
    reaction_v10_k10 = 0.1;
    reaction_v10 = compartment_intracellular*x(6)*x(8)*reaction_v10_k10;

    % Reaction: id = v11
    % Local Parameter:   id =  k11, name = k11
    reaction_v11_k11 = 5.0;
    reaction_v11 = compartment_intracellular*x(9)*reaction_v11_k11;

    % Reaction: id = v12
    % Local Parameter:   id =  k12, name = k12
    reaction_v12_k12 = 1.0;
    reaction_v12 = compartment_intracellular*x(10)*x(11)*reaction_v12_k12;

    % Reaction: id = v13
    % Local Parameter:   id =  k13, name = k13
    reaction_v13_k13 = 3.0;
    reaction_v13 = compartment_intracellular*x(12)*reaction_v13_k13;

    % Reaction: id = v14
    % Local Parameter:   id =  k14, name = k14
    reaction_v14_k14 = 1.0;
    reaction_v14 = compartment_intracellular*x(13)*x(14)*reaction_v14_k14;

    % Reaction: id = v15
    % Local Parameter:   id =  k15, name = k15
    reaction_v15_k15 = 3.0;
    reaction_v15 = compartment_intracellular*x(15)*reaction_v15_k15;

    % Reaction: id = v16
    % Local Parameter:   id =  k16, name = k16
    reaction_v16_k16 = 3.0;
    reaction_v16 = compartment_intracellular*x(12)*x(15)*reaction_v16_k16;

    % Reaction: id = v17
    % Local Parameter:   id =  k17, name = k17
    reaction_v17_k17 = 100.0;
    reaction_v17 = compartment_intracellular*x(8)*reaction_v17_k17;

    % Reaction: id = v18
    % Local Parameter:   id =  k18, name = k18
    reaction_v18_k18 = 5.0;
    reaction_v18 = compartment_intracellular*x(9)*x(16)*reaction_v18_k18;

    % Reaction: id = v19
    % Local Parameter:   id =  k19, name = k19
    reaction_v19_k19 = 1.0;
    reaction_v19 = compartment_intracellular*x(17)*reaction_v19_k19;

    % Reaction: id = v20
    % Local Parameter:   id =  k20, name = k20
    reaction_v20_k20 = 10.0;
    reaction_v20 = compartment_intracellular*x(17)*reaction_v20_k20;

    % Reaction: id = v21
    % Local Parameter:   id =  k21, name = k21
    reaction_v21_k21 = 5.0;
    reaction_v21 = compartment_intracellular*x(17)*reaction_v21_k21;

    % Reaction: id = v22
    % Local Parameter:   id =  k22, name = k22
    reaction_v22_k22 = 47.0;
    reaction_v22 = compartment_intracellular*x(18)*reaction_v22_k22;

    % Reaction: id = v23
    % Local Parameter:   id =  k23, name = k23
    reaction_v23_k23 = 5.0;
    reaction_v23 = compartment_intracellular*x(18)*reaction_v23_k23;

    % Reaction: id = v24
    % Local Parameter:   id =  k24, name = k24
    reaction_v24_k24 = 345.0;
    reaction_v24 = compartment_intracellular*x(19)*reaction_v24_k24;

    % Reaction: id = v25
    % Local Parameter:   id =  k25, name = k25
    reaction_v25_k25 = 5.0;
    reaction_v25 = compartment_intracellular*x(19)*reaction_v25_k25;

    % Reaction: id = v26
    % Local Parameter:   id =  k26, name = k26
    reaction_v26_k26 = 50.0;
    reaction_v26 = compartment_intracellular*x(20)*reaction_v26_k26;

    % Reaction: id = v27
    % Local Parameter:   id =  k27, name = k27
    reaction_v27_k27 = 5.0;
    reaction_v27 = compartment_intracellular*x(20)*reaction_v27_k27;

    % Reaction: id = v28
    % Local Parameter:   id =  k28, name = k28
    reaction_v28_k28 = 140.0;
    reaction_v28 = compartment_intracellular*x(21)*reaction_v28_k28;

    % Reaction: id = v29
    % Local Parameter:   id =  k29, name = k29
    reaction_v29_k29 = 10.0;
    reaction_v29 = compartment_intracellular*x(22)*x(14)*reaction_v29_k29;

    % Reaction: id = v30
    % Local Parameter:   id =  k30, name = k30
    reaction_v30_k30 = 1.0;
    reaction_v30 = compartment_intracellular*x(24)*reaction_v30_k30;

    % Reaction: id = v31
    % Local Parameter:   id =  k31, name = k31
    reaction_v31_k31 = 250.0;
    reaction_v31 = compartment_intracellular*x(24)*reaction_v31_k31;

    % Reaction: id = v32
    % Local Parameter:   id =  k32, name = k32
    reaction_v32_k32 = 5.0;
    reaction_v32 = compartment_intracellular*x(22)*reaction_v32_k32;

    % Reaction: id = v33
    % Local Parameter:   id =  k33, name = k33
    reaction_v33_k33 = 50.0;
    reaction_v33 = compartment_intracellular*x(23)*reaction_v33_k33;

    % Reaction: id = v34
    % Local Parameter:   id =  k34, name = k34
    reaction_v34_k34 = 18.0;
    reaction_v34 = compartment_intracellular*x(25)*x(23)*reaction_v34_k34;

    % Reaction: id = v35
    % Local Parameter:   id =  k35, name = k35
    reaction_v35_k35 = 10.0;
    reaction_v35 = compartment_intracellular*x(26)*reaction_v35_k35;

    % Reaction: id = v36
    % Local Parameter:   id =  k36, name = k36
    reaction_v36_k36 = 0.1;
    reaction_v36 = compartment_intracellular*x(26)*x(27)*reaction_v36_k36;

    % Reaction: id = v37
    % Local Parameter:   id =  k37, name = k37
    reaction_v37_k37 = 0.1;
    reaction_v37 = compartment_intracellular*x(28)*reaction_v37_k37;

    % Reaction: id = v38
    % Local Parameter:   id =  k38, name = k38
    reaction_v38_k38 = 0.01;
    reaction_v38 = compartment_intracellular*x(28)*reaction_v38_k38;

    % Reaction: id = v39
    % Local Parameter:   id =  k39, name = k39
    reaction_v39_k39 = 18.0;
    reaction_v39 = 0*compartment_intracellular*x(30)*x(23)*x(23)/(100*100+x(23)*x(23))*reaction_v39_k39; % 0*?

    % Reaction: id = v40
    % Local Parameter:   id =  k40, name = k40
    reaction_v40_k40 = 1.0;
    reaction_v40 = compartment_intracellular*x(31)*reaction_v40_k40;

    % Reaction: id = v41
    % Local Parameter:   id =  k41, name = k41
    reaction_v41_k41 = 0.02;
    reaction_v41 = compartment_intracellular*x(30)*x(35)*reaction_v41_k41;

    % Reaction: id = v42
    % Local Parameter:   id =  k42, name = k42
    reaction_v42_k42 = 0.1;
    reaction_v42 = compartment_intracellular*x(6)*x(31)*reaction_v42_k42;

    % Reaction: id = v43
    % Local Parameter:   id =  k43, name = k43
    reaction_v43_k43 = 0.01;
    reaction_v43 = compartment_intracellular*x(33)*reaction_v43_k43;

    % Reaction: id = v44
    % Local Parameter:   id =  k44, name = k44
    reaction_v44_k44 = 0.01;
    reaction_v44 = compartment_intracellular*x(34)*reaction_v44_k44;

    % Reaction: id = v45
    % Local Parameter:   id =  k45, name = k45
    reaction_v45_k45 = 0.1;
    reaction_v45 = compartment_intracellular*x(31)*x(35)*reaction_v45_k45;

    % Reaction: id = v46
    % Local Parameter:   id =  k46, name = k46
    reaction_v46_k46 = 200.0;
    reaction_v46 = compartment_intracellular*x(23)^(2)/(4^(2)+x(23)^(2))*reaction_v46_k46;

    % Reaction: id = v47
    % Local Parameter:   id =  k47, name = k47
    reaction_v47_k47 = 1.0;
    reaction_v47 = compartment_intracellular*x(36)*reaction_v47_k47;
    
    % A: the following reactions have been added and constants need to be
    % found
    
    % Reaction: id = v48
    % Local Parameter:   id =  k48, name = k48
    % reaction_v48_k48 = 0.0; 
    % Protein production proportional to Ste12 concentration (should we also do it in two steps?)
    reaction_v48 = compartment_intracellular*x(25)*reaction_v48_k48; 
  
    % Reaction: id = 49
    % Local Parameter:   id =  k49, name = k49
    reaction_v49_k49 = reaction_v48_k48; 
    % Production of alpha-factor
    reaction_v49 = compartment_intracellular*x(25)*reaction_v49_k49; 
    % A: Don't really understand the thing with the compartments
    
    % A: NOTE: Our measurements will not include removal due to colon movement, 
    % so should optimize without it?
    % Reaction: id = 50 (removal of aplha pheromone)
    % Local Parameter: id = k50, name = k50
    reaction_v50_k50 = 0.0; % Should find in literature
    reaction_v50 = compartment_extracellular*x(1)*reaction_v50_k50;
    
    
    % Reaction: id = 50 (removal of aplha pheromone)
    % Local Parameter: id = k50, name = k50
    reaction_v51_k51 = reaction_v50_k50; 
    reaction_v51 = compartment_extracellular*x(38)*reaction_v51_k51;
    
    % ODEs
    xdot = zeros(nr_species,1);
    
    % Species:   id = alpha, name = alpha-factor, affected by kineticLaw
    % A: added production
    % NOTE: Assume all alpha-factor is secreted
    xdot(1) = (pop_size/(compartment_extracellular))*((-1.0 * reaction_v1) + (1.0 * reaction_v49) + (-1.0 * reaction_v50));

    % Species:   id = Ste2, name = Ste2, affected by kineticLaw
    xdot(2) = (1/(compartment_intracellular))*((-1.0 * reaction_v2) + ( 1.0 * reaction_v3) + (-1.0 * reaction_v5));

    % Species:   id = Ste2a, name = Ste2active, affected by kineticLaw
    xdot(3) = (1/(compartment_intracellular))*(( 1.0 * reaction_v2) + (-1.0 * reaction_v3) + (-1.0 * reaction_v4));

    % Species:   id = Gabc, name = G???, affected by kineticLaw
    xdot(4) = (1/(compartment_intracellular))*((-1.0 * reaction_v6) + ( 1.0 * reaction_v9));

    % Species:   id = GaGTP, name = G?GTP, affected by kineticLaw
    xdot(5) = (1/(compartment_intracellular))*(( 1.0 * reaction_v6) + (-1.0 * reaction_v7) + (-1.0 * reaction_v8));

    % Species:   id = Gbc, name = G??, affected by kineticLaw
    xdot(6) = (1/(compartment_intracellular))*(( 1.0 * reaction_v6) + (-1.0 * reaction_v9) + (-1.0 * reaction_v10) + ( 1.0 * reaction_v11) + ( 1.0 * reaction_v21) + ( 1.0 * reaction_v23) + ( 1.0 * reaction_v25) + ( 1.0 * reaction_v27) + ( 1.0 * reaction_v32) + (-1.0 * reaction_v42) + ( 1.0 * reaction_v43));

    % Species:   id = GaGDP, name = G?GDP, affected by kineticLaw
    xdot(7) = (1/(compartment_intracellular))*(( 1.0 * reaction_v7) + ( 1.0 * reaction_v8) + (-1.0 * reaction_v9));

    % Species:   id = complexC, name = complexC, affected by kineticLaw
    xdot(8) = (1/(compartment_intracellular))*((-1.0 * reaction_v10) + ( 1.0 * reaction_v11) + ( 1.0 * reaction_v16) + (-1.0 * reaction_v17));

    % Species:   id = complexD, name = complexD, affected by kineticLaw
    xdot(9) = (1/(compartment_intracellular))*(( 1.0 * reaction_v10) + (-1.0 * reaction_v11) + (-1.0 * reaction_v18) + ( 1.0 * reaction_v19));

    % Species:   id = Ste5, name = Ste5, affected by kineticLaw
    xdot(10) = (1/(compartment_intracellular))*((-1.0 * reaction_v12) + ( 1.0 * reaction_v13) + ( 1.0 * reaction_v17) + ( 1.0 * reaction_v21) + ( 1.0 * reaction_v23) + ( 1.0 * reaction_v25) + ( 1.0 * reaction_v27) + ( 1.0 * reaction_v32));

    % Species:   id = Ste11, name = Ste11, affected by kineticLaw
    xdot(11) = (1/(compartment_intracellular))*((-1.0 * reaction_v12) + ( 1.0 * reaction_v13) + ( 1.0 * reaction_v17) + ( 1.0 * reaction_v21) + ( 1.0 * reaction_v23) + ( 1.0 * reaction_v25) + ( 1.0 * reaction_v27) + ( 1.0 * reaction_v32));

    % Species:   id = complexA, name = complexA, affected by kineticLaw
    xdot(12) = (1/(compartment_intracellular))*(( 1.0 * reaction_v12) + (-1.0 * reaction_v13) + (-1.0 * reaction_v16));

    % Species:   id = Ste7, name = Ste7, affected by kineticLaw
    xdot(13) = (1/(compartment_intracellular))*((-1.0 * reaction_v14) + ( 1.0 * reaction_v15) + ( 1.0 * reaction_v17) + ( 1.0 * reaction_v21) + ( 1.0 * reaction_v23) + ( 1.0 * reaction_v25) + ( 1.0 * reaction_v27) + ( 1.0 * reaction_v32));

    % Species:   id = Fus3, name = Fus3, affected by kineticLaw
    xdot(14) = (1/(compartment_intracellular))*((-1.0 * reaction_v14) + ( 1.0 * reaction_v15) + ( 1.0 * reaction_v17) + ( 1.0 * reaction_v21) + ( 1.0 * reaction_v23) + ( 1.0 * reaction_v25) + ( 1.0 * reaction_v27) + (-1.0 * reaction_v29) + ( 1.0 * reaction_v30) + ( 1.0 * reaction_v33));

    % Species:   id = complexB, name = complexB, affected by kineticLaw
    xdot(15) = (1/(compartment_intracellular))*(( 1.0 * reaction_v14) + (-1.0 * reaction_v15) + (-1.0 * reaction_v16));

    % Species:   id = Ste20, name = Ste20, affected by kineticLaw
    xdot(16) = (1/(compartment_intracellular))*((-1.0 * reaction_v18) + ( 1.0 * reaction_v19) + ( 1.0 * reaction_v21) + ( 1.0 * reaction_v23) + ( 1.0 * reaction_v25) + ( 1.0 * reaction_v27) + ( 1.0 * reaction_v32));

    % Species:   id = complexE, name = complexE, affected by kineticLaw
    xdot(17) = (1/(compartment_intracellular))*(( 1.0 * reaction_v18) + (-1.0 * reaction_v19) + (-1.0 * reaction_v20) + (-1.0 * reaction_v21));

    % Species:   id = complexF, name = complexF, affected by kineticLaw
    xdot(18) = (1/(compartment_intracellular))*(( 1.0 * reaction_v20) + (-1.0 * reaction_v22) + (-1.0 * reaction_v23));

    % Species:   id = complexG, name = complexG, affected by kineticLaw
    xdot(19) = (1/(compartment_intracellular))*(( 1.0 * reaction_v22) + (-1.0 * reaction_v24) + (-1.0 * reaction_v25));

    % Species:   id = complexH, name = complexH, affected by kineticLaw
    xdot(20) = (1/(compartment_intracellular))*(( 1.0 * reaction_v24) + (-1.0 * reaction_v26) + (-1.0 * reaction_v27));

    % Species:   id = complexI, name = complexI, affected by kineticLaw
    xdot(21) = (1/(compartment_intracellular))*(( 1.0 * reaction_v26) + (-1.0 * reaction_v28) + ( 1.0 * reaction_v31));

    % Species:   id = complexL, name = complexL, affected by kineticLaw
    xdot(22) = (1/(compartment_intracellular))*(( 1.0 * reaction_v28) + (-1.0 * reaction_v29) + ( 1.0 * reaction_v30) + (-1.0 * reaction_v32));

    % Species:   id = Fus3PP, name = Fus3PP, affected by kineticLaw
    xdot(23) = (1/(compartment_intracellular))*(( 1.0 * reaction_v28) + (-1.0 * reaction_v33) + (-1.0 * reaction_v34) + ( 1.0 * reaction_v35));

    % Species:   id = complexK, name = complexK, affected by kineticLaw
    xdot(24) = (1/(compartment_intracellular))*(( 1.0 * reaction_v29) + (-1.0 * reaction_v30) + (-1.0 * reaction_v31));

    % Species:   id = Ste12, name = Ste12, affected by kineticLaw
    xdot(25) = (1/(compartment_intracellular))*((-1.0 * reaction_v34) + ( 1.0 * reaction_v35));

    % Species:   id = Ste12a, name = Ste12active, affected by kineticLaw
    xdot(26) = (1/(compartment_intracellular))*(( 1.0 * reaction_v34) + (-1.0 * reaction_v35));

    % A: REMOVE BAR1 FROM GENOME (set change of species 27, 28, 29 to 0.0)
    % Species:   id = Bar1, name = Bar1, affected by kineticLaw
    % xdot(27) = (1/(compartment_intracellular))*((-1.0 * reaction_v36) + ( 1.0 * reaction_v37));
    xdot(27) = 0.0; % Or just remove?
    
    % Species:   id = Bar1a, name = Bar1active, affected by kineticLaw
    % xdot(28) = (1/(compartment_intracellular))*(( 1.0 * reaction_v36) + (-1.0 * reaction_v37) + (-1.0 * reaction_v38));
    xdot(28) = 0.0;
    
    % Species:   id = Bar1aex, name = Bar1activeEx, affected by kineticLaw
    % xdot(29) = (1/(compartment_extracellular))*(( 1.0 * reaction_v38));
    xdot(29) = 0.0;
    
    % Species:   id = Far1, name = Far1, affected by kineticLaw
    xdot(30) = (1/(compartment_intracellular))*((-1.0 * reaction_v39) + ( 1.0 * reaction_v40) + (-1.0 * reaction_v41));

    % Species:   id = Far1PP, name = Far1PP, affected by kineticLaw
    xdot(31) = (1/(compartment_intracellular))*(( 1.0 * reaction_v39) + (-1.0 * reaction_v40) + (-1.0 * reaction_v42) + ( 1.0 * reaction_v43) + ( 1.0 * reaction_v44) + (-1.0 * reaction_v45));

    % Species:   id = Far1U, name = Far1ubiquitin, affected by kineticLaw
    xdot(32) = (1/(compartment_intracellular))*(( 1.0 * reaction_v41));

    % Species:   id = complexM, name = complexM, affected by kineticLaw
    xdot(33) = (1/(compartment_intracellular))*(( 1.0 * reaction_v42) + (-1.0 * reaction_v43));

    % Species:   id = complexN, name = complexN, affected by kineticLaw
    xdot(34) = (1/(compartment_intracellular))*((-1.0 * reaction_v44) + ( 1.0 * reaction_v45));

    % Species:   id = Cdc28, name = Cdc28, affected by kineticLaw
    xdot(35) = (1/(compartment_intracellular))*(( 1.0 * reaction_v44) + (-1.0 * reaction_v45));

    % Species:   id = Sst2, name = Sst2, affected by kineticLaw
    xdot(36) = (1/(compartment_intracellular))*(( 1.0 * reaction_v46) + (-1.0 * reaction_v47) + (1.0 * reaction_v51));

    % Species:   id = p, name = p
    % WARNING speciesID: p, constant= false, boundaryCondition = true but is not involved in assignmentRule, rateRule or events !
    xdot(37) = 0.0; 
    
    % A: The following species have been added
    % Species:   id = aca, name = anti-cancer agent
    xdot(38) = (1/(compartment_intracellular))*(( 1.0 * reaction_v48));
end