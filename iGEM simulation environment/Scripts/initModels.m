function [superModel] = initModels(models)
% initModels
%   Runs initial FBA on each model in models structure. Creates superModel 
%   structure.
%
%   models          Cell structure containing models
%
%   superModel
%      Sub Model ID       ID of organism
%      Rection Name       Vector containing exchange reaction names
%      Reaction ID        Vector containing exchange reaction IDs          
%      Reaction Index     Vector containing exchange reaction indexes
%      Exchange Fluxes    Vector containing exchange flux vlaues  
%      Bounds             Vector containing lower and upper bounds
%      Conversion Factor  Based on biomass of each GEM
%      FBA Solution       Structure containing FBA solution of model
%      Sub Models         GEMs initialized by initModels()
%
%   Usage: [superModel, exEnvironment] = initModels(models)
%
%   Francisco Zorrilla, 05-08-2018


exComp = [];
for z =1:length(models) % find extracellular compartment in each model
   for o=1:length(models{z}.compNames)
        if strcmpi(models{z}.compNames(o),'extracellular')
             exComp(z) = o;
        elseif strcmpi(models{z}.compNames(o),'e')
             exComp(z) = o;
        end
   end
   
   %before mergComp, get IDs/idxs of all metabolites present in extracellular comp, then add transport rxns for these from [s] to [e] after compMerge
   exMetId = [];
   for r=1:length(models{z}.metComps)
       if models{z}.metComps(r) == exComp(z)
           exMetId(r) = true;
       else
           exMetId(r) = false;
       end
   end
   exMetId = transpose(logical(exMetId));
   
   %metsToAdd = {};  %INSTEAD OF ADDING METS AGAIN, ADD EXCH RXNS FOR 'e' METS 
   %metsToAdd = struct;
   %metsToAdd.mets = models{z}.mets(exMetId);
   %metsToAdd.metNames = models{z}.metNames(exMetId);
   %metsToAdd.compartments = 'e'; %this will be fed into addMets(), look at help file to add more parameters 
   
   if isempty(exComp(z)) % IF ERROR IN THIS LINE, READ EM BELLOW
        EM=['The extracellular compartement in model' models{z}.id ' could not be identified'];
        dispEM(EM,false);
   else %rename extracellular compartments for consistency
        models{z}.comps(exComp(z))= {'e'};
        models{z}.compNames(exComp(z))= {'extracellular'};
   end
   
   [models{z}, addedRxns]=addExchangeRxns(models{z},'both',models{z}.mets(exMetId));
   
   %merge compartments in each GEM for simpler models
   models{z} = mergeCompartments(models{z}); 
   models{z}.comps(1,1) = {['c_' models{z}.id]};
   models{z}.comps(2,1) = {['e']};
   models{z}.compNames(1,1) = {[models{z}.id ' Cytoplasm']};
   models{z}.compNames(2,1) = {['Extracellular']};
   
   %I = ismember(metsToAdd.mets,models{z}.mets); %make sure mets to add are not present in model
   %metsToAdd.mets = metsToAdd.mets(~I);
   %metsToAdd.metNames = metsToAdd.metNames(~I);
   %models{z}=addMets(models{z},metsToAdd); % works but had to comment out 2 different checks (metName and metId) from addMets fucntion
   

   
   if isfield(models{z},'compOutside') % delete this field if present
       models{z} = rmfield(models{z},'compOutside');
   end
end

ids = {};
exchangeRxns = {};
exchangeRxnsNames = {};
exchangeRxnsIndexes = {};
FBAsol = {};
exchangeFluxes = {};
bounds = {};
conversionFactor = {};

%populate cells for superModel structure
for q = 1:length(models) 
    [exchangeRxns{q},exchangeRxnsIndexes{q}]=getExchangeRxns(models{q}); %get exchange reaction IDs and their indexes
    exchangeRxnsNames{1,q} = transpose({models{q}.rxnNames{exchangeRxnsIndexes{q}}}); % get exchange reaction names
    FBAsol{q} = solveLP(models{q},1); %solve initial FBA with default parameters in models
    ids{q} = models{q}.id; %get model ids
    exchangeFluxes{q} = FBAsol{q}.x(exchangeRxnsIndexes{q}); %get exchange flux values
    bounds{q} = {models{q}.lb models{q}.ub}; %get bounds for each exchange reaction
    conversionFactor{q} = abs(FBAsol{q}.f); %get biomass flux for conversion factor calculation
end

%populate superModel structure
superModel = struct;
superModel.OrganismID = ids;
superModel.RxnName = exchangeRxnsNames;
superModel.RxnID = exchangeRxns;
superModel.RxnIndex = exchangeRxnsIndexes;
superModel.ExchangeFluxes = exchangeFluxes;
superModel.Bounds = bounds;
superModel.ConversionFactor = conversionFactor;
superModel.FBAsol = FBAsol;
superModel.SubModels = models;

% check FBA solutions and constrain fluxes
for w=1:length(models) 
    fixLB = [];
    fixUB = [];
   for r=1:length(superModel.ExchangeFluxes{w})
       if superModel.ExchangeFluxes{w}(r)>= 5   %randomly chose 5 for constraint, check with Ed and literature values
           fixUB(r) = true;
       else
           fixUB(r) = false;
       end
       if superModel.ExchangeFluxes{w}(r) <= -5
           fixLB(r) = true;
       else
           fixLB(r) = false;
       end
       fixLB = logical(transpose(fixLB));
       fixUB = logical(transpose(fixUB));
   end
   superModel.SubModels{w}.lb(superModel.RxnIndex{w}(fixLB)) = -5;
   superModel.SubModels{w}.ub(superModel.RxnIndex{w}(fixUB)) = 5;
end

%update superModel structure with new FBA solutions based on adjusted
%bounds
FBAsol = {}; %clear these cells to avoid using old information
bounds = {};
exchangeFluxes = {};

for y=1:length(models)
    FBAsol{y} = solveLP(superModel.SubModels{y},1); %solve FBA with new bounds
    bounds{y} = {superModel.SubModels{y}.lb superModel.SubModels{y}.ub}; %populate with new bounds
    exchangeFluxes{y} = FBAsol{y}.x(exchangeRxnsIndexes{y}); %get exchange flux values
end

superModel.FBAsol = FBAsol;
superModel.Bounds = bounds;
superModel.ExchangeFluxes = exchangeFluxes;

%identify blocked reactions in models to send warning
%blockedRxns = {};
%
%for w=1:length(models) 
%   for r=1:length(superModel.ExchangeFluxes{w})
%       if superModel.SubModels{w}.lb(r) == 0 && superModel.SubModels{w}.ub(r) == 0
%           blockedRxns{r}= superModel.RxnName{w}(r);
%       end
%   end
%   superModel.BlockedRxns{w} = [blockedRxns];
%end


end