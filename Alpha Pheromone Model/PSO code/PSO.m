clear all;
clc;

% Read data from file
exprData = xlsread('ExprData_iGEMDelft2012.xlsx'); % iGEM Delft (2012)
% Convert days to minutes
exprData(:, 1) = exprData(:, 1) * 24 * 60;
% Remove irrelevant data
exprData = exprData(1:210, [1, 6]);
% Scale data
exprData(2:end, 2) = exprData(2:end, 2) - exprData(2, 2);
exprData(2:end, 2) = exprData(2:end, 2)/max(max(exprData(2:end, 2)));

% Define parameters
numberOfParticles = 40;
dimensions = 4;
timeStepLength = 1;
runningTime = 1000; 
inertiaWeight = 1.4;
inertiaWeightMinimum = 0.3;
inertiaWeightBeta = 0.99;
positionMin = 1e-12; 
positionMax = 10000; 
maximumVelocity = 0.3 * (positionMax - positionMin) / timeStepLength; 

bestPositions = zeros(numberOfParticles, dimensions);
particleMinima = Inf * ones(numberOfParticles, 1);
globalBestPosition = zeros(dimensions, 1);
globalMinimum = Inf;

% Initialize particles
particles = InitializeParticles(numberOfParticles, dimensions, timeStepLength, positionMin, positionMax);
functionValues = zeros(numberOfParticles, 1);

% Optimize
for t = 1:timeStepLength:runningTime
    disp(t)
    
    for i = 1:numberOfParticles
        functionValues(i) = EvaluateParticle(particles(i).Position, exprData);

        if (functionValues(i) < particleMinima(i))
            particleMinima(i) = functionValues(i);
            bestPositions(i,:) = particles(i).Position;
        end

        if (functionValues(i) < globalMinimum)
            globalMinimum = functionValues(i);
            globalBestPosition = particles(i).Position;
        end
    end

    particles = UpdateParticles(particles, bestPositions, globalBestPosition, timeStepLength, inertiaWeight, maximumVelocity, positionMin, positionMax);

    if (inertiaWeight > inertiaWeightMinimum)
        inertiaWeight = inertiaWeight * inertiaWeightBeta;
    end
end

disp(globalBestPosition)
disp(globalMinimum)



