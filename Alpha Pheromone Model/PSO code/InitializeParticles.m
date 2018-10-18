function particles = InitializeParticles(numberOfParticles, dimensions, timeStepLength, positionMin, positionMax)
% Initializes positions and velocities of the particles in a PSO
    
    alpha = 1;  
    positionRange = positionMax - positionMin;
  
    particles = [];
    for i = 1:numberOfParticles
        position = zeros(dimensions, 1);
        velocity = zeros(dimensions, 1);
        maxFactor = -10;

        for j = 1:dimensions
            r1 = rand;
            r2 = rand * maxFactor;
            position(j) = positionMin + r1 * positionRange * 10^r2;
            
            r = rand; 
            velocity(j) = (alpha / timeStepLength) * ((-positionRange / 2) + r * positionRange);
        end

        tmpParticle = struct('Position', position, 'Velocity', velocity);
        particles = [particles; tmpParticle];
    end
    
end
