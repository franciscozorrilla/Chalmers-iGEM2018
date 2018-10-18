function x = DecodeChromosome(chromosome, nrVars)

    maxFactor = -8;

    nrGenes = size(chromosome, 2);
    nrParts = floor(nrGenes / nrVars);
    
    x = zeros(nrVars, 1);
    if nrParts == 1
        x = abs(chromosome);
    else 
        for i = 1:nrVars
            factor = 10;
            for j = 1:(nrParts - 1)
                x(i) = x(i) + factor * chromosome((i - 1) * nrParts + j);
                factor = factor / 10;
            end
            x(i) = x(i) * 10^(maxFactor * chromosome(i * nrParts));
        end
    end
    
end
        
    