% Replace proteinPool with glutamine 

%exportToExcelFormat(modelCancerTest,'modelCancerTest.xlsx');

modelCancer=setParam(modelCancer,'eq','EXC_IN_m90000c',0);
modelCancer=setParam(modelCancer,'ub','HMR_9063',1000); %glutamine
modelCancer=setParam(modelCancer,'ub','HMR_9038',1000); %histidine
modelCancer=setParam(modelCancer,'ub','HMR_9041',1000); %lysine
modelCancer=setParam(modelCancer,'ub','HMR_9043',1000); %phenylalanine
modelCancer=setParam(modelCancer,'ub','HMR_9046',1000); %valine
modelCancer=setParam(modelCancer,'ub','HMR_9044',1000); %threonine
modelCancer=setParam(modelCancer,'ub','HMR_9045',1000); %tryptophan
modelCancer=setParam(modelCancer,'ub','HMR_9042',1000); %methionine
modelCancer=setParam(modelCancer,'ub','HMR_9040',1000); %leucine
modelCancer=setParam(modelCancer,'ub','HMR_9039',1000); %isoleucine

solCancerTest=solveLP(modelCancer,1);
printFluxes(modelCancer,solCancerTest.x);
checkRxn(modelCancer,'human_proteinPool')


%% CANCER MODEL CURATION

% Cancer model - growth achieved, but bounds seem to be unconstrained
modelCancer = load('HepG2.mat');
modelCancer = modelCancer.HepG2model;
modelCancer.id = 'Can'
for q=1:length(modelCancer.subSystems)
    modelCancer.subSystems{q,1}={modelCancer.subSystems{q,1}}; %Need to do this for exportToExcel() to work
end
modelCancer = rmfield(modelCancer,'unconstrained') % Need to remove these fields for exportToExcel() to work
modelCancer = rmfield(modelCancer,'metMiriams')
modelCancer = rmfield(modelCancer,'inchis')
modelCancer = rmfield(modelCancer,'rxnComps')
modelCancer = rmfield(modelCancer,'geneComps')
%exportToExcelFormat(modelCancer,'modelCancer.xlsx');
modelCancer=setParam(modelCancer,'obj','HumanGrowth',1);
%modelCancer.ub(8177:8192) = 5; % CONSTRAIN INNER FLUXES, MAYBE DO WITH FUNCTIONS?
%modelCancer=setParam(modelCancer,'eq',{'HMR_9729','HMR_9067','HMR_9044'},0);%try to constrain carbon source to glucse, maybe not needed
solCancer=solveLP(modelCancer,1);
printFluxes(modelCancer,solCancer.x);

modelCancerMC = mergeCompartments(modelCancer);

[modelCancerMC, addedRxns]=addExchangeRxns(modelCancerMC,'in','m90000c'); % add protein pool exchange rxn
canRxns= getExchangeRxns(modelCancerMC,'in'); modelCancerMC = setParam(modelCancerMC, 'eq', canRxns, 0); %constrain all exchange rxns goin in 
modelCancerMC = setParam(modelCancerMC,'lb',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079','HMR_9063','EXC_IN_m90000c'},-1);
modelCancerMC = setParam(modelCancerMC,'ub',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079','HMR_9063','EXC_IN_m90000c'},1);
modelCancerMC = setParam(modelCancerMC,'lb',{'HumanGrowth','human_proteinPool','humanGrowthOut'},0);
modelCancerMC = setParam(modelCancerMC,'ub',{'HumanGrowth','human_proteinPool','humanGrowthOut'},1);
modelCancerMC = removeRxns(modelCancerMC,{'HMR_1592','HMR_1696','HMR_4271','HMR_0006','HMR_0007','HMR_0008','HMR_0015','HMR_0016','HMR_0017','HMR_1919','HMR_7568','HMR_7569','HMR_2141','HMR_4184','HMR_2585'});%these intercompartment transport rxns may become problematic after merging compartments
solCancerMC=solveLP(modelCancerMC,1);
printFluxes(modelCancerMC,solCancerMC.x); % soln is a bit different than non merged compartment model
