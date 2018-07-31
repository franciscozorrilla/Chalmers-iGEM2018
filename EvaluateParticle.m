function objectiveFunctionValue = EvaluateParticle(particlePosition)
% Finds the objective function of a particle
    alpha_factor_conc_measured = KofahlANDKlipp2004Model_AntiCancerYeast(0.1);
    % This is just for testing, we should add our own measurements here
    % Should we do interpolation to be able to compare the obtained values
    % with our observed values?
    % https://se.mathworks.com/help/matlab/interpolation.html
    
    reaction_v48_k48 = particlePosition(1);
    alpha_factor_conc = KofahlANDKlipp2004Model_AntiCancerYeast(reaction_v48_k48);
    
    objectiveFunctionValue = sum((alpha_factor_conc_measured - alpha_factor_conc).^2);
    % Or other objective function?
end