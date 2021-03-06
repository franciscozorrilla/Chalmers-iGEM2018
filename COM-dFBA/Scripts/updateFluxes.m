function [superModel,fluxes,bound,lbub] = updateFluxes(superModel,params)
%updateFluxes
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
%   params
%      Vmax               Matrix with Vmax information for uptake of
%                         metabolite m by organism n
%      Ks                 Matrix with Ks information for uptake of
%                         metabolite m by organism n
%      mets               Matrix with initial concentration and flow
%                         concentration of metabolites
%      other              Vector with miscellaneous terms not associated
%                         with one particular organism or metabolite
%      metIdx             Array of exchange reaction equation indexes of
%                         metabolites present in metNames and models in 
%                         superModel
%      blocked            Arbitrary reaction with lower bound and upper
%                         bound set to zero found in each model in
%                         superModel. Used as an index when model is not
%                         producing/consuming metabolite from metNames
%
%   superModelConstrained Updated superModel structure
%
%   fluxesConstrained     Array with one cell per organism. Each cell
%                         contains a matrix with 7 columns and one row per
%                         exchange reaction in the FBA solution for that 
%                         origanism. Columns contain information about 
%                         rxnIndex rxnName fluxValue rxnID rxnEqn lb ub
% 
%   boundsConstrained     Contains information about updated bound
%  
%   lbubConstrained       Contains ifnormation about which bound was
%                         updated. 'u' means the upper bound was updated, 
%                         'l' means the lower bound was updated, and 'x' 
%                         means that no bound was updated.
%
%   Usage: [superModelConstrained,fluxesConstrained,boundsConstrained,lbubConstrained] = updateFluxes(superModel,params);
%
%   Francisco Zorrilla, 18-09-2018

Vmax = params.Vmax;
Ks = params.Ks;
mets = params.mets;
metIdx = params.metIdx;
blocked = params.blocked;

initialConditions = transpose(params.mets(:,1));

bound = [];
%Define bound for each organism and uptake reaction
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(Vmax)
        if (Vmax(upMet,organism)==0 & Ks(upMet,organism) ==0) | initialConditions(upMet)==0
            bound(upMet,organism)=0;
        else
            bound(upMet,organism) = abs((Vmax(upMet,organism) * initialConditions(upMet)) / (Ks(upMet,organism) + initialConditions(upMet)));
        end
    end
end

equationStrings = {};
%Construct equation strings for each exchange metabolite
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(metIdx)
        if metIdx(upMet,organism) ~= blocked(1,organism)
            equationStrings{upMet,organism}=constructEquations(superModel.subModels{organism},metIdx(upMet,organism));
        end
    end
end


lbub = {};
%Define lower bound/upper bound structure containing information about
%which bound needs to update: only for uptake metabolites
for organism = 1:length(superModel.organismID) %do this for each organism
    for upMet = 1:length(metIdx) % do this for each uptake metabolite
        if metIdx(upMet,organism) == blocked(1,organism) % check if model is not uptaking this metabolite
            lbub(upMet,organism) = {'x'};
        else
            if (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) > abs(superModel.subModels{organism}.ub(metIdx(upMet,organism))))&(strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))>0.5)
                lbub(upMet,organism) = {'l'};
            elseif (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) < abs(superModel.subModels{organism}.ub(metIdx(upMet,organism))))&(strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))<0.5)
                lbub(upMet,organism) = {'u'};
            elseif (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) == abs(superModel.subModels{organism}.ub(metIdx(upMet,organism)))) & (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) ~= 0)
                dispEM(['Unable to determine whether lower or upper bound should be modified for organism' superModel.organismID{organism} ' ,reaction involving metabolite number ' upMet '. Errors may arise if this occurs for reactions other than water.'],false)
            else
                lbub(upMet,organism) = {'x'};
            end
        end
    end
end

%Update bounds based previously defined objects 'bound','lbub', and 'equationStrings'
for organism = 1:length(superModel.organismID)
    for upMet = (length(superModel.organismID)+1):length(Vmax)
        if strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))>0.5 & lbub{upMet,organism}=='u'  %met on rxntnt side & lb<ub: secreted metabolite
            %do nothing
        elseif strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))<0.5 & lbub{upMet,organism}=='u'  %met on product side & lb<ub: uptake metabolite
            superModel.subModels{organism} = setParam(superModel.subModels{organism},'ub',superModel.subModels{organism}.rxns{metIdx(upMet,organism)},bound(upMet,organism));%update bound
        elseif strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))>0.5 & lbub{upMet,organism}=='l'  %met on rxntnt side & lb>ub: uptake metabolite
            superModel.subModels{organism} = setParam(superModel.subModels{organism},'lb',superModel.subModels{organism}.rxns{metIdx(upMet,organism)},-bound(upMet,organism));%update bound
        elseif strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))<0.5 & lbub{upMet,organism}=='l'  %met on product side & lb<ub: secreted metabolite
            %do nothing
        end
    end
end


ids = {};
exchangeRxns = {};
exchangeRxnsNames = {};
exchangeRxnsIndexes = {};
FBAsol = {};
exchangeFluxes = {};
bounds = {};

for q = 1:length(superModel.organismID) 
    if mets(q,1)~=0 | mets(q,2)~=0
        [exchangeRxns{q},exchangeRxnsIndexes{q}]=getExchangeRxns(superModel.subModels{q}); %get exchange reaction IDs and their indexes
        exchangeRxnsNames{1,q} = transpose({superModel.subModels{q}.rxnNames{exchangeRxnsIndexes{q}}}); % get exchange reaction names
        FBAsol{q} = solveLP(superModel.subModels{q},1); %solve initial FBA with default parameters in models
        exchangeFluxes{q} = FBAsol{q}.x(exchangeRxnsIndexes{q}); %get exchange flux values
        bounds{q} = {superModel.subModels{q}.lb superModel.subModels{q}.ub}; %get bounds for each exchange reaction
    else
        exchangeRxns{q}=0;
        exchangeRxnsIndexes{q}=0;
        exchangeRxnsNames{1,q} =0;
        FBAsol{q}=0;
        exchangeFluxes{q}=0;
        bounds{q} =0;       
    end
end
%populate superModel structure
superModel.rxnName = exchangeRxnsNames;
superModel.rxnID = exchangeRxns;
superModel.rxnIndex = exchangeRxnsIndexes;
superModel.exchangeFluxes = exchangeFluxes;
superModel.bounds = bounds;
superModel.FBAsol = FBAsol;


exch = {};
fluxes = {};
eqns = {}; 
lb = {};
ub = {};
idx = {};

for t=1:length(superModel.organismID)
    if mets(t,1)~=0 | mets(t,2)~=0
        eqns{t} = constructEquations(superModel.subModels{t},superModel.rxnID{t},true);
        lb{t} = superModel.subModels{t}.lb(superModel.rxnIndex{t});
        ub{t} = superModel.subModels{t}.ub(superModel.rxnIndex{t});
        idx{t} = superModel.rxnIndex{t};
        exch{t} = [ num2cell(idx{t}) superModel.rxnName{t} num2cell(superModel.exchangeFluxes{t}) superModel.rxnID{t} eqns{t} num2cell(lb{t}) num2cell(ub{t})];
        fluxes{t} = exch{t}(find(superModel.exchangeFluxes{t}),:);
    else
        eqns{t} = 0;
        lb{t} = 0;
        ub{t} = 0;
        idx{t} = 0;
        exch{t} = 0;
        fluxes{t} = 0;
    end
end

end