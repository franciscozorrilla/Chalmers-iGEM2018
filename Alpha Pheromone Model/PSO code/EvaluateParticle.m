function objectiveFunctionValue = EvaluateParticle(particlePosition, data)
% Finds the objective function of a particle

    tPts = data(2:end, 1);
    exprData = data(2:end, 2); 
    alphaConc = data(1, 2);
    
    v48_k48 = particlePosition(1);
    v49_k49 = particlePosition(2); 
    v50_k50 = particlePosition(3);
    v51_k51 = particlePosition(4);
    
    exprModel = zeros(size(exprData));
    for i = 1:length(alphaConc)
        [t, x] = alpha_pheromone_system(alphaConc(i), v48_k48, v49_k49, v50_k50, v51_k51, 1, tPts(end));
        exprModel(:, i) = interp1(t, x, tPts, 'neighbour');
    end
    
    exprModel = exprModel / max(max(exprModel));
    exprModel = exprModel - exprModel(1);
    
    objectiveFunctionValue = sum((exprData - exprModel).^2, 'all');
    
end