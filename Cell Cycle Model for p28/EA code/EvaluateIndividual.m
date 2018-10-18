function f = EvaluateIndividual(x)
    
    tf = 86400;
    p28_0 = [0, 5, 50, 100, 200];
    dataPts = [1.0, 0.87, 0.91, 0.89, 0.81]; % Yamada et al. (2009)  
    modelPts = zeros(size(dataPts));
    
    source_p53 = x(1);
    d_p53 = x(2);
    f1 = x(3); 
    d_Mdm2 = x(4);
    f2 = x(5);
    d_cop1 = x(6);
    k_deg_cop1 = x(7);
    uptake_rate = x(8);
    deg_basal = x(9);
    d_p28 = x(10);
    k_67_a = x(11);

    flag = 0;
    nrReps = 10;
    for i = 1:size(dataPts, 2)
        
        if i == 1
            mpf_peak = tf;
        end
        
        for k = 1:nrReps
            [prol, peak, err] = p28_action(p28_0(i), tf, mpf_peak, source_p53,...
                d_p53, f1, d_Mdm2, f2, d_cop1, k_deg_cop1, uptake_rate, deg_basal, d_p28, k_67_a);
            modelPts(:, i) = modelPts(:, i) + prol;
            if err == 1
               flag = 1;
               break;
            end
        end
        if flag == 1
            break;
        end
        
        if i == 1
            mpf_peak = peak; % As baseline for cell cycle arrest
        end
    end
    
    if flag == 0
        modelPts = modelPts / nrReps;
        f = 1 / sum(sum((dataPts - modelPts).^2));
        % To favor variations in proliferation
        f = f + length(unique(modelPts)) - 1;
    else
        f = 0;
    end
    
end