function mutatedChromosome = Mutate(chromosome, mutationProbability)

    nGenes = size(chromosome, 2);
    mutatedChromosome = chromosome;
    
    Cr = 0.1;
    for j = 1:nGenes
        r = rand;
        if (r < (mutationProbability / 2))
            mutatedChromosome(j) = rand;
        elseif ((mutationProbability / 2) < r < mutationProbability)
            mutatedChromosome(j) = mutatedChromosome(j) - Cr / 2 + Cr * r;
        end
    end
    
end