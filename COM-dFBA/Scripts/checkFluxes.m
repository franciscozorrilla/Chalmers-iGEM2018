function [metEquations, metIdx, bounds, blocked] = checkFluxes(superModel,metNames)
%checkFluxes
%   superModel
%      Sub Model ID       ID of organism
%      Rection Name       Vector containing exchange reaction names
%      Reaction ID        Vector containing exchange reaction IDs          
%      Reaction Index     Vector containing exchange reaction indexes
%      Exchange Fluxes    Vector containing exchange flux vlaues  
%      Bounds             Vector containing lower and upper bounds
%      FBA Solution       Structure containing FBA solution of model
%      Sub Models         GEMs initialized by initModels()
%   
%   metNames              Cell vector containing string name of each
%                         metabolite to be matched between models in
%                         superModel
%                         
%   metEquations          Array of cells containing exchnage reaction
%                         equations of metabolites present in metNames and
%                         models in superModel
%
%   metIdx                Array of exchange reaction equation indexes of
%                         metabolites present in metNames and models in 
%                         superModel
%
%   bounds                Array of exchange reaction lower/upper bounds for
%                         each metabolite metNames and present in superModel      
%
%   blocked               Arbitrary reaction with lower bound and upper
%                         bound set to zero found in each model in
%                         superModel. Used as an index when model is not
%                         producing/consuming metabolite from metNames
%  
%   Usage: [metEquations, metIdx, blocked] = checkFluxes(superModel,metNames)
%
%   Francisco Zorrilla, 12-09-2018


metIdx = [];
%Find exchange reaction index for each defined metabolite in each model
for organism = 1:length(superModel.organismID)
    tempModel = superModel.subModels{organism};
    [exchangeRxns,exchangeRxnsIndexes]=getExchangeRxns(tempModel,'both');
    equationStrings=constructEquations(tempModel,exchangeRxnsIndexes);
    for upMet = 1:length(metNames)
        for orgFlux = 1:length(exchangeRxnsIndexes)
            if contains(equationStrings{orgFlux},metNames{upMet},'IgnoreCase',true) & (length(equationStrings{orgFlux}) <= (length(metNames{upMet})+10))
                metIdx(upMet,organism) = exchangeRxnsIndexes(orgFlux);
            end
        end
    end
end

%find a reactions with lower and upper bound equal to zero in each model
blocked = zeros(1,6);
for organism = 1:length(superModel.organismID)
    while blocked(1,organism) ==0
        for rxn = 1:length(superModel.subModels{organism}.rxns)
            if superModel.subModels{organism}.lb(rxn) == 0 & superModel.subModels{organism}.ub(rxn) == 0
                blocked(1,organism) = rxn;
            end
        end
    end
end

%set metIdx equal to blocked reaction for metabolites not uptaken or
%secreted by model
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(metIdx)
        if metIdx(upMet,organism) == 0
            metIdx(upMet,organism) = blocked(1,organism);
        end
    end
end

equationStrings = {};
%check exchange reactions
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(metIdx)
        if metIdx(upMet,organism) ~= blocked(1,organism)
            equationStrings{upMet,organism}=constructEquations(superModel.subModels{organism},metIdx(upMet,organism));
        end
    end
end
metEquations = equationStrings;

bounds = [];
%check exchange reactions
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(metIdx)
        if metIdx(upMet,organism) ~= blocked(1,organism)
            bounds{upMet,organism} = [superModel.subModels{organism}.lb(metIdx(upMet,organism)) superModel.subModels{organism}.ub(metIdx(upMet,organism))];
        end
    end
end