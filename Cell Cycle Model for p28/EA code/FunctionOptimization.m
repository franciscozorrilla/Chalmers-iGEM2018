% Code implemented based on the course Stochastic Optimization Algorithms, Chalmers University of Technology (Whade, M., 2017)

clear all;
clc;

% Define parameters
populationSize = 50;
nrVars = 11;
numberOfGenes = nrVars * 5;
crossoverProbability = 0.75;
mutationProbability = 0.1; 
minimumMutationProbability = 0.025;
mutationWeight = 0.99; 
tournamentSelectionParameter = 0.75;
numberOfGenerations = 1000;
variableRange = 1.0;
fitness = zeros(populationSize, 1);

% Initialize population
decodedPopulation = zeros(populationSize, nrVars);
population = InitializePopulation(populationSize, numberOfGenes);
%population(:, 1) = dlmread('bestie');

globalMaximumFitness = 0.0;
xBestGlobal = zeros(1, nrVars);
bestChromosomeGlobal = zeros(1, numberOfGenes);

for iGeneration = 1:numberOfGenerations

    maximumFitness = 0.0; 
    xBest = zeros(1, nrVars);
    bestIndividualIndex = 0;
    
    for i = 1:populationSize
        
        chromosome = population(i, :);
        x = DecodeChromosome(chromosome, nrVars);
        decodedPopulation(i, :) = x;
        fitness(i) = EvaluateIndividual(x);
        
        if (fitness(i) > maximumFitness)
            maximumFitness = fitness(i);
            bestIndividualIndex = i;
            xBest = x;
        end
        
        if (fitness(i) > globalMaximumFitness)
            globalMaximumFitness = fitness(i);
            xBestGlobal = x;
            bestChromosomeGlobal = population(i, :);
        end
    end   
    
    disp('xBest')
    disp(xBest)
    disp('maximumFitness')
    disp(maximumFitness)
    % To check that the initial concentrations are ok
    disp(InitConc(xBest))
    
    tempPopulation = population;
    
    for i = 1:2:populationSize
        i1 = TournamentSelect(fitness, tournamentSelectionParameter);
        i2 = TournamentSelect(fitness, tournamentSelectionParameter);
        
        chromosome1 = population(i1,:);
        chromosome2 = population(i2,:);

        r = rand;
        if (r < crossoverProbability)
            newChromosomePair = Cross(chromosome1, chromosome2);
            tempPopulation(i,:) = newChromosomePair(1,:);
            tempPopulation(i+1,:) = newChromosomePair(2,:);
        else
            tempPopulation(i,:) = chromosome1;
            tempPopulation(i+1,:) = chromosome2;
        end
    end 
    
    for i = 1:populationSize
        originalChromosome = tempPopulation(i,:);
        mutatedChromosome = Mutate(originalChromosome, mutationProbability);
        tempPopulation(i,:) = mutatedChromosome;
    end
    
    tempPopulation(1,:) = population(bestIndividualIndex,:);
    tempPopulation(2,:) = bestChromosomeGlobal;
    population = tempPopulation;
    
    if mutationProbability > minimumMutationProbability
        mutationProbability = mutationWeight * mutationProbability;
    end
    
end
 
disp('xBest')
disp(xBest)
disp('maximumFitness')
disp(maximumFitness)
    
