function [superModel] = initModels(models,whichFluxes)
% initModels
%   Runs initial FBA on each model in models structure. Creates superModel 
%   structure.
%
%   models          Cell structure containing models
%   whichFluxes         Use 1 to print exchange fluxes
%                       Use 2 to print all fluxes
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
%   Usage: [superModel] = initModels(models,whichFluxes)
%
%   Francisco Zorrilla, 12-09-2018

ids = {};
exchangeRxns = {};
exchangeRxnsNames = {};
exchangeRxnsIndexes = {};
FBAsol = {};
exchangeFluxes = {};
bounds = {};

if whichFluxes == 1
    %populate cells for superModel structure RUN THIS TO GET ONLY EXCHANGE FLUXES
    for q = 1:length(models) 
        [exchangeRxns{q},exchangeRxnsIndexes{q}]=getExchangeRxns(models{q}); %get exchange reaction IDs and their indexes
        exchangeRxnsNames{1,q} = transpose({models{q}.rxnNames{exchangeRxnsIndexes{q}}}); % get exchange reaction names
        FBAsol{q} = solveLP(models{q},1); %solve initial FBA with default parameters in models
        ids{q} = models{q}.id; %get model ids
        exchangeFluxes{q} = FBAsol{q}.x(exchangeRxnsIndexes{q}); %get exchange flux values
        bounds{q} = {models{q}.lb models{q}.ub}; %get bounds for each exchange reaction
    end
    %populate superModel structure
    superModel = struct;
    superModel.organismID = ids;
    superModel.rxnName = exchangeRxnsNames;
    superModel.rxnID = exchangeRxns;
    superModel.rxnIndex = exchangeRxnsIndexes;
    superModel.exchangeFluxes = exchangeFluxes;
    superModel.bounds = bounds;
    superModel.FBAsol = FBAsol;
    superModel.subModels = models;
end

if whichFluxes ==2
    %populate cells for superModel structure RUN THIS TO GET ALL FLUXES
    for q = 1:length(models) 
        FBAsol{q} = solveLP(models{q},1); %solve initial FBA with default parameters in models
        ids{q} = models{q}.id; %get model ids
        exchangeFluxes{q} = FBAsol{q}.x; %get exchange flux values
        bounds{q} = {models{q}.lb models{q}.ub}; %get bounds for each exchange reaction
        rxnsNames{q} = models{q}.rxnNames;
        rxns{q} = models{q}.rxns;
        exchangeRxnsIndexes{q} = transpose(1:length(models{q}.rxns)); 
    end
    %populate superModel structure
    superModel = struct;
    superModel.organismID = ids;
    superModel.rxnName = rxnsNames;
    superModel.rxnID = rxns;
    superModel.rxnIndex = exchangeRxnsIndexes;
    superModel.exchangeFluxes = exchangeFluxes;
    superModel.bounds = bounds;
    superModel.FBAsol = FBAsol;
    superModel.subModels = models;
end

end