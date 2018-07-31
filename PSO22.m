clear all;
clc;

% OBSTACLE TO SOLVE: ODE23 DOES NOT GIVE A FIXED STEP SIZE, SO HOW DO ONE
% OPTIMIZE??? COULD USE NLINFIT AS DELFT 2012 DID BUT THAT WOULD NOT BE AS
% FUN :'( (https://se.mathworks.com/help/stats/nlinfit.html)

numberOfParticles = 40;
dimensions = 1;
timeStepLength = 1;
runningTime = 1000; 
inertiaWeight = 1.4;
inertiaWeightMinimum = 0.3;
inertiaWeightBeta = 0.99;
positionMin = -10; % Will start with this
positionMax = 10;
maximumVelocity = (positionMax - positionMin) / timeStepLength;

bestPositions = zeros(numberOfParticles, dimensions);
particleMinima = Inf * ones(numberOfParticles, 1);
globalBestPosition = zeros(dimensions, 1);
globalMinimum = Inf;

particles = InitializeParticles(numberOfParticles, dimensions, timeStepLength, positionMin, positionMax);
functionValues = zeros(numberOfParticles, 1);

for t = 1:timeStepLength:runningTime
    for i = 1:numberOfParticles

        functionValues(i) = EvaluateParticle(particles(i).Position);

        if (functionValues(i) < particleMinima(i))
            particleMinima(i) = functionValues(i);
            bestPositions(i,:) = particles(i).Position;
        end

        if (functionValues(i) < globalMinimum)
            globalMinimum = functionValues(i);
            globalBestPosition = particles(i).Position;
        end
    end

    particles = UpdateParticles(particles, bestPositions, globalBestPosition, timeStepLength, inertiaWeight, maximumVelocity);

    if (inertiaWeight > inertiaWeightMinimum)
        inertiaWeight = inertiaWeight * inertiaWeightBeta;
    end
end

disp(sprintf('x: %.6f',globalBestPosition(1)));
disp(sprintf('y: %.6f',globalBestPosition(2)));
disp(sprintf('Function value: %.6f', globalMinimum));



