exportToExcelFormat(modelColon,'modelColon.xlsx');% Export to excel for inspection, no biomass function present! Obtain from HMR2 model

load('HMRdatabase2_00.mat');% HMR2.0

for q=1:length(ihuman.subSystems)
    ihuman.subSystems{q,1}={ihuman.subSystems{q,1}}; %Need to do this for exportToExcel() to work
end

exportToExcelFormat(ihuman,'ihuman.xlsx')%Try adding all 'artificial reactions' from HMR2.0 to colon

%super sketchy move
modelColon.genes = ihuman.genes

noGeneIdx=find(cellfun(@isempty,ihuman.grRules)); % Find rxns with no genes
rxnIdx=regexp(ihuman.rxnNames,'glucose'); % Find rxns with key words of interest in name
rxnIdx=find(~cellfun('isempty',rxnIdx)); % Find indices
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns=ihuman.rxns(rxnIdx); % Obtain reaction IDs
modelSbo=addRxnsGenesMets(modelSbo,modelSce,rxns,false,...
    'Modeling reaction required for growth, gene unknown',1); % Add reactions and metabolites to model
%undo super sketchy move
%modelColon = rmfield(modelColon,'genes')

%delete subSystems for exportToExcel() to work
modelColon = rmfield(modelColon,'subSystems')
exportToExcelFormat(modelColon,'modelColon.xlsx');

modelColon=setParam(modelColon,'obj','biomass_components',1);
sol=solveLP(modelColon,1);
printFluxes(modelColon,sol.x); % no growth, check bounds of exchange rxns

modelColon = setParam(modelColon,'ub',{'HMR_9034','HMR_9456'},1); %still no growth
