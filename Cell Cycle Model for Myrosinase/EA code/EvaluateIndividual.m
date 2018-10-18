function f = EvaluateIndividual(x)
    
    tf = 86400;
    sulphoraphane_0 = [0, 10, 20, 40, 80];
	dataPts = [1.0, 0.88, 0.76, 0.59, 0.29]; % Liu et al. (2017)
	modelPts = zeros(size(dataPts));
    
    kbax_m = x(1);
    kcytc_m = x(2);
    kp21_m = x(3); 
    uptake_rate = x(4);
    
    flag = 0;
    nrReps = 100;
    mpf_peak = tf;
    for i = 1:size(dataPts, 2)
        for k = 1:nrReps
            [prol, peak, err] = myrosinase_action(sulphoraphane_0(i), tf, mpf_peak, kbax_m, kcytc_m, kp21_m, uptake_rate);
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
            mpf_peak = peak;
        end
    end
    
    if flag == 0
        modelPts = modelPts / nrReps;
        f = 1 / sum(sum((dataPts - modelPts).^2));
        % To favor variations in cell proliferation
        f = f + length(unique(modelPts)) - 1;       
    else
        f = 0;
    end
    
end