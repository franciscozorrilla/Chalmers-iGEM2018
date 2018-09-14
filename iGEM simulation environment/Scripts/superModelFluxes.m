function fluxes = superModelFluxes(superModel,whichFluxes)
%SuperModelFluxes 
%   superModel          Output from initModels() function
%   whichFluxes         Use 1 to print exchange fluxes
%                       Use 2 to print all fluxes
%
%   Usage: fluxes = superModelFluxes(superModel,whichFluxes)
%
%   Francisco Zorrilla, 12-09-2018

exch = {};
fluxes = {};
eqns = {}; 
lb = {};
ub = {};
idx = {};
subSys = {};

if whichFluxes == 1
%RUN THIS TO GET EXCHANGE FLUXES
    for t=1:length(superModel.organismID)
        eqns{t} = constructEquations(superModel.subModels{t},superModel.rxnID{t},true);
        lb{t} = superModel.subModels{t}.lb(superModel.rxnIndex{t});
        ub{t} = superModel.subModels{t}.ub(superModel.rxnIndex{t});
        idx{t} = superModel.rxnIndex{t};
        exch{t} = [ num2cell(idx{t}) superModel.rxnName{t} num2cell(superModel.exchangeFluxes{t}) superModel.rxnID{t} eqns{t} num2cell(lb{t}) num2cell(ub{t})];
        fluxes{t} = exch{t}(find(superModel.exchangeFluxes{t}),:);
    end
end

if whichFluxes == 2
%RUN THIS TO GET ALL FLUXES
    for t=1:length(superModel.organismID)
        eqns{t} = constructEquations(superModel.subModels{t},superModel.rxnID{t},true);
        lb{t} = superModel.subModels{t}.lb;
        ub{t} = superModel.subModels{t}.ub;
        idx{t} = superModel.rxnIndex{t};
        subSys{t} = superModel.subModels{t}.subSystems;
        exch{t} = [ num2cell(idx{t}) superModel.rxnID{t} superModel.rxnName{t} num2cell(superModel.exchangeFluxes{t}) eqns{t} num2cell(lb{t}) num2cell(ub{t}) subSys{t}];
        fluxes{t} = exch{t}(find(superModel.exchangeFluxes{t}),:);
    end
end

end

