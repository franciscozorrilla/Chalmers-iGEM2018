function out = superModelFluxes(superModel)
%SuperModelFluxes 
%   superModel          Output from initModels() function
%

exch = {};
out = {};
eqns = {}; 
lb = {};
ub = {};

for t=1:length(superModel.OrganismID)
    eqns{t} = constructEquations(superModel.SubModels{t},superModel.RxnID{t},true);
    lb{t} = superModel.SubModels{t}.lb(superModel.RxnIndex{t});
    ub{t} = superModel.SubModels{t}.ub(superModel.RxnIndex{t});
    exch{t} = [ superModel.RxnName{t} num2cell(superModel.ExchangeFluxes{t}) superModel.RxnID{t} eqns{t} num2cell(lb{t}) num2cell(ub{t})];
    out{t} = exch{t}(find(superModel.ExchangeFluxes{t}),:);
end


end

