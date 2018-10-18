function population = InitializePopulation(populationSize, nGenes)

    population = zeros(populationSize, nGenes);
    for i = 1:populationSize
        for j = 1:nGenes
            population(i, j) = rand;  
        end
    end
    
end