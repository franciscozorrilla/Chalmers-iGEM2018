%% Colon Model Curation

modelColon = importModel('colon.xml');
modelColon.id = 'Col';

modelColonMC = mergeCompartments(modelColon);

modelColon = removeReactions(modelColon,{'HMR_1696','HMR_1849','HMR_0006','HMR_0007','HMR_0008','HMR_0015','HMR_1919','HMR_7568','HMR_7569','HMR_2585'}) %these intercompartment transport rxns may become problematic after merging compartments
colRxns= getExchangeRxns(modelColon,'both'); modelColon = setParam(modelColon, 'eq', colRxns, 0);%CHANGE TO 'IN' INSTEAD OF 'BOTH'
modelColon = setParam(modelColon,'lb',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079'},-5);
modelColon = setParam(modelColon,'ub',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079'},5);
modelColon = setParam(modelColon,'eq','HMR_4149',1);%TCA cycle rxn
modelColon = setParam(modelColon,'obj','HMR_4394',1);%glycolysis pathway rxn

solColon = solveLP(modelColon,1);
printFluxes(modelColon,solColon.x);


%% BALANCE OXPHOS RXN ID HMR_6916

modelColonT = modelColon imported without exchange fluxes and then merged compartments

modelColonT = setParam(modelColonT,'obj','HMR_4149',1)

clear rxnToAdd metsToAdd % Make sure these structures are empty

rxnToAdd.rxns={'OXPHOS_BALANCED'};
rxnToAdd.rxnNames={'Oxidative Phosphorylation (Balanced)'};
rxnToAdd.equations={'ADP[s] + Pi[s] + 0.5 O2[s] + 2 H+[s] => ATP[s] + H2O[s]'};
rxnToAdd.lb=[0];
rxnToAdd.ub=[1000];
modelColonT=addRxns(modelColonT,rxnToAdd,3,'',true);

modelColonT= setParam(modelColonT,'eq','HMR_6916',0);

solColT = solveLP(modelColonT,1);
printFluxes(modelColonT,solColT.x);
%% Tried moving biomass equation from HMR2.0 to Colon model for growth, didnt work

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


