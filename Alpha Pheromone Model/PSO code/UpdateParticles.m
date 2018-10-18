function particles = UpdateParticles(particles, bestPositions, globalBestPosition, timeStepLength, inertiaWeight, maximumVelocity, positionMin, positionMax)
% Updates the positions and velocities of the particles
    
    cognitiveComponentConstant = 2;
    socialComponentConstant = 2;
    
    numberOfParticles = size(particles, 1); 
    dimensions = size(particles(1).Position, 1);
    pCr = 0.01;
    
    for i = 1:numberOfParticles
        q = rand;
        r = rand;
        for j = 1:dimensions
            cognitiveComponent = cognitiveComponentConstant * q * (bestPositions(i,j) - particles(i).Position(j)) / timeStepLength;
            socialComponent = socialComponentConstant * r * (globalBestPosition(j) - particles(i).Position(j)) / timeStepLength;
            
            particles(i).Velocity(j) = inertiaWeight * particles(i).Velocity(j) + cognitiveComponent + socialComponent;
		
            r = rand;          
            if r > pCr
                particles(i).Velocity(j) = - maximumVelocity + 2 * rand * maximumVelocity;
            end
                
            if (abs(particles(i).Velocity(j)) > maximumVelocity)
                particles(i).Velocity(j) = sign(particles(i).Velocity(j)) * maximumVelocity; 
            end
		
            particles(i).Position(j) = particles(i).Position(j) + particles(i).Velocity(j) * timeStepLength;
            
            if (particles(i).Position(j) > positionMax) 
                particles(i).Position(j) = positionMax;
            elseif (particles(i).Position(j) < positionMin)
                particles(i).Position(j) = positionMin;
            end
            
        end
    end
    
end