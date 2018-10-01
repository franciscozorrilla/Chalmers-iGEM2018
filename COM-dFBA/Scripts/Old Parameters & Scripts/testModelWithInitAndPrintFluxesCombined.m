modelSbo=setParam(modelSbo,'lb','r_4066',0); %set MFalpha2 production
modelSbo=setParam(modelSbo,'lb','r_4067',0.002); %set Myrosinase production

models = {modelSbo}

%populate cells for superModel structure RUN THIS TO GET ALL FLUXES
for q = 1:length(models) 
    FBAsol{q} = solveLP(models{q},1); %solve initial FBA with default parameters in models
    ids{q} = models{q}.id; %get model ids
    exchangeFluxes{q} = FBAsol{q}.x; %get exchange flux values
    bounds{q} = {models{q}.lb models{q}.ub}; %get bounds for each exchange reaction
    conversionFactor{q} = abs(FBAsol{q}.f); %get biomass flux for conversion factor calculation
    rxnsNames{q} = models{q}.rxnNames;
    rxns{q} = models{q}.rxns;
    exchangeRxnsIndexes{q} = transpose(1:length(models{q}.rxns)); 
end
%populate superModel structure
superModel = struct;
superModel.OrganismID = ids;
superModel.RxnID = rxns;
superModel.ExchangeFluxes = exchangeFluxes;
superModel.Bounds = bounds;
superModel.ConversionFactor = conversionFactor;
superModel.FBAsol = FBAsol;
superModel.SubModels = models;
superModel.RxnName = rxnsNames;
superModel.RxnIndex = exchangeRxnsIndexes;

exch = {};
out = {};
eqns = {}; 
lb = {};
ub = {};
idx = {};
subSys = {};
%RUN THIS TO GET ALL FLUXES
for t=1:length(superModel.OrganismID)
    eqns{t} = constructEquations(superModel.SubModels{t},superModel.RxnID{t},true);
    lb{t} = superModel.SubModels{t}.lb;
    ub{t} = superModel.SubModels{t}.ub;
    idx{t} = superModel.RxnIndex{t};
    subSys{t} = superModel.SubModels{t}.subSystems;
    exch{t} = [ num2cell(idx{t}) superModel.RxnID{t} superModel.RxnName{t} num2cell(superModel.ExchangeFluxes{t}) eqns{t} num2cell(lb{t}) num2cell(ub{t}) subSys{t}];
    out{t} = exch{t}(find(superModel.ExchangeFluxes{t}),:);
end

out{1}