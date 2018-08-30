%% Simulation environment for iGEM Chalmers 2018
%  Author: Francisco Zorrilla
%  2018-08-05
%
%  Community model of Human Gut + Colorectal Cancer + Engineered
%  S.boulardii. 

%{ 

  The human gut is modelled using an HMR2.0 derived
  tissue specific model for the colon. The model was obtained from: 
  http://www.metabolicatlas.org/downloads/tissue

  In addition, we include 3 representative human gut bateria: Bacteroides 
  thetaiotamicron (Bth), Eubacterium rectale (Ere) and Methanobrevibacter 
  smithii (Msi). These models were obtained from: 
  http://www.metabolicatlas.org/downloads/micro

  The Colorectal Cancer (CRC) model was obtained from a collegue, Jonathan
  Robbinson.

  The Saccharomyces boulardii model (Sbo) was derived from homology using a
  curated S.cerevisiae model as a template. Its reconstruction process is
  detailed in "reconstruction.m".

%}

cd('C:\Users\zorrilla\Desktop\iGEM simulation environment')

%% 1. Load & Check Genome Scale Models
% Please ensure that each model being loaded has had some degree of
% curation, and can be used for simulation (i.e. has objective function,
% set bounds, etc). This will be done in the following sections.

%% Saccharomyces boulardii model - growth achieved

load('modelSbo.mat'); 
modelSbo.id = 'Sbo';
modelSbo=setParam(modelSbo,'obj','r_4041',1); %set objective function

modelSbo=setParam(modelSbo,'ub','r_1654',0); % make ammonium uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_1714',0); % make glucose uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_1992',0); % make oxygen uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_2005',0); % make phosphate uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_2060',0); % make sulphate uptake irreversible

modelSbo=setParam(modelSbo,'lb','r_4066',0.001); %set MFalpha2 production to a minimum of 0.009, this seems to be the maximum it can produce
%however, once the model becomes constrained by competing with other
%organisms in the dFBA simulations, it seems like the solution become
%infeasible. Try instead with a more lenient value of 0.001

modelSbo=mergeCompartments(modelSbo);%merge compartments to make simulations run faster
solSbo=solveLP(modelSbo,1);
printFluxes(modelSbo,solSbo.x);
%exportToExcelFormat(modelSbo,'modelSbo.xlsx');

%% Bacteroides thetaiotamicron model - growth achieved

modelBth = importModel('iBth801 v1.00.xml');
modelBth.id = 'Bth';
modelBth=setParam(modelBth,'obj','biomass_red',1);
modelBth=setParam(modelBth,'lb','glcIn',0); % make glucose uptake irreversible
modelBth=mergeCompartments(modelBth);
solBth=solveLP(modelBth,1);
printFluxes(modelBth,solBth.x);
exportToExcelFormat(modelBth,'modelBth.xlsx');  

%% Eubacterium rectale model - growth achieved

modelEre = importModel('iEre400 v1.00.xml');
modelEre.id = 'Ere';
modelEre=setParam(modelEre,'obj','r55',1);
modelEre=setParam(modelEre,'lb','glcIn',0); % make glucose uptake irreversible
modelEre=setParam(modelEre,'ub','r113',0); % make water uptake irreversible
modelEre=setParam(modelEre,'lb','r114',0); % make phosphate uptake irreversible
modelEre=mergeCompartments(modelEre);
solEre=solveLP(modelEre,1);
printFluxes(modelEre,solEre.x);
exportToExcelFormat(modelEre,'modelEre.xlsx');

%% Cancer model - growth achieved, but bounds seem to be unconstrained

modelCancer = load('HepG2.mat');
modelCancer = modelCancer.HepG2model;
modelCancer.id = 'Can';

for q=1:length(modelCancer.subSystems)
    modelCancer.subSystems{q,1}={modelCancer.subSystems{q,1}}; %Need to do this for exportToExcel() to work
end
modelCancer = rmfield(modelCancer,'unconstrained'); % Need to remove these fields for exportToExcel() and printFluxes() functions to work
modelCancer = rmfield(modelCancer,'metMiriams'); 
modelCancer = rmfield(modelCancer,'inchis');
modelCancer = rmfield(modelCancer,'rxnComps');
modelCancer = rmfield(modelCancer,'geneComps');

modelCancer=setParam(modelCancer,'obj','HumanGrowth',1);
[modelCancer, addedRxns]=addExchangeRxns(modelCancer,'in','m90000c'); % add protein pool exchange rxn
canRxns= getExchangeRxns(modelCancer,'in'); %get exchange rxns going in
modelCancer = setParam(modelCancer, 'eq', canRxns, 0); %constrain all exchange rxns going in 

modelCancer = setParam(modelCancer,'lb',{'HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079','EXC_IN_m90000c'},-1); %allow for uptake/secretion of some mets
modelCancer = setParam(modelCancer,'ub',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079','EXC_IN_m90000c','HumanGrowth','human_proteinPool','humanGrowthOut'},1); % allow for excretion/consumption of some mets
modelCancer = setParam(modelCancer,'lb',{'HumanGrowth','human_proteinPool','humanGrowthOut','HMR_9047','HMR_9048','HMR_9058','HMR_9034','HMR_9072'},0);%constrain growth related rxns
modelCancer = setParam(modelCancer,'ub',{'HMR_9073','HMR_9078','HMR_9079'},0);

modelCancer = removeRxns(modelCancer,{'HMR_1592','HMR_1696','HMR_4271','HMR_0006','HMR_0007','HMR_0008','HMR_0015','HMR_0016','HMR_0017','HMR_1919','HMR_7568','HMR_7569','HMR_2141','HMR_4184','HMR_2585'});%these intercompartment transport rxns may become problematic after merging compartments
modelCancer = mergeCompartments(modelCancer);
solCancer=solveLP(modelCancer,1);
printFluxes(modelCancer,solCancer.x);
exportToExcelFormat(modelCancer,'modelCancer.xlsx');

%% Colon tissue model 

modelColon = importModel('colon.xml');
modelColon.id = 'Col';
modelColon = removeReactions(modelColon,{'HMR_1696','HMR_1849','HMR_0006','HMR_0007','HMR_0008','HMR_0015','HMR_1919','HMR_7568','HMR_7569','HMR_2585'}); %these intercompartment transport rxns may become problematic after merging compartments
colRxns= getExchangeRxns(modelColon,'in'); % get exchange fluxes going in
modelColon = setParam(modelColon, 'eq', colRxns, 0); %constrain fluxes going in
modelColon = setParam(modelColon,'lb',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079'},-1);
modelColon = setParam(modelColon,'ub',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079'},1);
modelColon = mergeCompartments(modelColon);
modelColon.genes = {''}; %insert dummy field to avoid errors with addRxns()

clear metsToAdd % Make sure this structure is empty
metsToAdd = {}; %Add proteinPool pseudo metabolite
metsToAdd.mets = {'proteinPool'};
metsToAdd.metNames = {'proteinPool'};
metsToAdd.compartments = {'s'};
modelColon=addMets(modelColon,metsToAdd);

clear rxnToAdd % Make sure this structure is empty
rxnToAdd.rxns={'OXPHOS_BALANCED','EXC_IN_m90000c','human_proteinPool','human_Maintenance'};%add balanced oxidative phosporylation rxn, AA exchange rxn and artificial AA pool rxn with same IDs as cancer model
rxnToAdd.rxnNames={'Oxidative Phosphorylation (Balanced)','EXC_IN_m90000c','human_proteinPool','human_Maintenance'};
rxnToAdd.equations={'ADP[s] + Pi[s] + 0.5 O2[s] + 2 H+[s] => ATP[s] + H2O[s]',' => proteinPool[s]','0.0937 alanine[s] + 0.0507 arginine[s] + 0.04 asparagine[s] + 0.04 aspartate[s] + 0.0142 cysteine[s] + 0.1118 glutamine[s] + 0.1831 glycine[s] + 0.0198 histidine[s] + 0.0309 isoleucine[s] + 0.0664 leucine[s] + 0.0571 lysine[s] + 0.0156 methionine[s] + 0.029 phenylalanine[s] + 0.0853 proline[s] + 0.0491 serine[s] + 0.0402 threonine[s] + 0.0072 tryptophan[s] + 0.019 tyrosine[s] + 0.0471 valine[s] <=> proteinPool[s]','ATP[s] + 6 H2O[s] + proteinPool[s] => ADP[s] + Pi[s] + H+[s]'};
rxnToAdd.lb=[0,0,-1,0];
rxnToAdd.ub=[1000,1,1,1000];
modelColon=addRxns(modelColon,rxnToAdd,3,'',true);

modelColon = setParam(modelColon,'obj','human_Maintenance',1);
modelColon= setParam(modelColon,'eq',{'HMR_6916','OXPHOS_BALANCED'},0); %block unbalanced oxphos rxn, also had to block balanced oxphos 
modelColon = rmfield(modelColon,'rxnComps');% only need to do after merging compartments for exportToExcel() to work

modelColon = setParam(modelColon,'ub',{'HMR_9047'},0)
modelColon = setParam(modelColon,'lb',{'HMR_9034','HMR_9048','HMR_9058','HMR_9078','HMR_9079'},0)

solColon = solveLP(modelColon,1);
printFluxes(modelColon,solColon.x);
exportToExcelFormat(modelColon,'modelColon.xlsx');

%% Methanobrevibacter smithii model
% This model gave me a lot of trouble. It seems to contain many
% gaps/isolated subnetworks, as can been seen from running
% gapReport(modelMsi). To try to make this a functioning model I had to
% unconstrain all reactions and create a biomass equation without the
% metabolites that could not be obtained. To get methane production, I had
% to introduce methane as a bi-product of growth.

modelMsi = importModel('iMsi385.xml');
modelMsi.id = 'Msi';
modelMsi=mergeCompartments(modelMsi);
%gapReport(modelMsi); Model seems to have many gaps. Press Control+C to quit gap report function if it gets stuck at ***Metabolite connectivity 
modelMsi=setParam(modelMsi,'obj','r449',1);

for q=1:length(modelMsi.rxnNames) % This is a fast and easy way to uncosntrain all reactions
   modelMsi.ub(q) = 1000;
   modelMsi.lb(q) = -1000;
end
checkRxn(modelMsi,'r449'); %this tells us that even when all reactions are unconstrained, model still cannot produce certain substrates needed for the biomass equation
%try recreating the biomass equation without these metabolites present as substrates

modelMsi.genes = {''}; %insert dummy field to avoid errors with addRxns()
clear rxnToAdd % Make sure this structure is empty
rxnToAdd.rxns={'biomassMsi'};
rxnToAdd.rxnNames={'biomassMsi'};
rxnToAdd.equations={'0.5621 L-Alanine[s] + 0.3237 L-Arginine[s] + 0.2638 L-Asparagine[s] + 0.2638 L-Aspartate[s] + 50.2 ATP[s] + 0.1002 L-Cysteine[s] + 0.0331 dATP[s] + 0.0215 dGTP[s] + 0.288 L-Glutamine[s] + 0.288 L-Glutamate[s] + 0.6704 Glycine[s] + 0.2222 GTP[s] + 50.2 H2O[s] + 0.3179 L-Isoleucine[s] + 0.493 L-Leucine[s] + 0.3755 L-Lysine[s] + 0.1682 L-Methionine[s] + 0.2027 L-Phenylalanine[s] + 0.2419 L-Proline[s] + 0.2361 L-Serine[s] + 0.2776 L-Threonine[s] + 0.1509 L-Tyrosine[s] + 0.4631 L-Valine[s] => 50.2 ADP[s] + 50.2 Phosphate[s] + Biomass_msi[s] + 0.5 Methane[s]'};
rxnToAdd.lb=[0];
rxnToAdd.ub=[1000];
modelMsi=addRxns(modelMsi,rxnToAdd,3,'',true);
% Note: M smithii produces methane gas, but the methane synthesis pathway
% seems to be disconnected from the biomass pathway. A fast way to fix this
% was to add methane as a product of biomass formation. An arbitrary
% stoiochiometric coefficient of 0.5 was used.

modelMsi = setParam(modelMsi,'lb','r225',0); % Need to constrain this reaction to be irreversible so that methane is excreted
msiRxns= getExchangeRxns(modelMsi,'in'); % get exchange fluxes going in
modelMsi = setParam(modelMsi, 'eq', msiRxns, 0); %constrain fluxes going in

modelMsi = setParam(modelMsi,'ub',{'actIn','NH4IN','piIN','H2Suptake','H2In','CO2In'},1);
modelMsi = setParam(modelMsi,'lb',{'H2ObTran'},-1);
modelMsi = setParam(modelMsi,'ub',{'H2ObTran'},0);
modelMsi = setParam(modelMsi,'lb',{'MethOut','actIn','H2In','CO2In'},0);

modelMsi = setParam(modelMsi,'obj','biomassMsi',1);
solMsi=solveLP(modelMsi,1);
printFluxes(modelMsi,solMsi.x); %now model can grow and produce methane
exportToExcelFormat(modelMsi,'modelMsi.xlsx');

%% Clean up workspace
clear addedRxns canRxns colRxns q rxnToAdd

%% 2. Initialize Models & Environment

%define cell strucutre of models to be initialized
models = {modelSbo modelBth modelEre modelMsi modelCancer modelColon};

%inspect superModel structure, refer to initModels() function file for details
[superModel] = initModels(models); 

% contains summary of all exchange rxns
out = superModelFluxes(superModel);


%% 3. Simulations

%{

curves look better when models have their respective compartments merged. 
when I do it without merging then Bth outgrows the shit out of everything 
else. also simulations take waaaaay longer if comps not merged, probably 
because doing so greatly simplifies the cancer and colon models. ode113 
seems to work the best for simulations, ode15s also works well sometimes.

BIOMASS
x(1):S.bo biomass
x(2):B.th biomass
x(3):E.re biomass
x(4):M.si biomass
x(5):Cancer biomass
x(6):Colon biomass
SUBSTRATES
x(7):Glucose
x(8):Amino Acid Pool
x(9):Water
x(10):Oxygen
x(11):Phosphate
x(12):Ammonium
PRODUCTS
x(13):Carbon dioxide
x(14):Propanoate
x(15):Butyrate
x(16):Methane
x(17):MFalpha2
x(18):Myrosinase

%}

%Glucose uptake kinetic parameters
%vmax_g1
%vmax_g2
%vmax_g3
%vmax_g5
%vmax_g6
%ks_g1
%ks_g2
%ks_g3
%ks_g5
%ks_g6

%


odeoptions = odeset('NonNegative',7:18);
tic
[t,xa] = ode15s(@(t,x)f(t,x,superModel),[0 10],[5 15 15 5 5 15 1 1 1 1 1 1 1 0 0 0 0 0],odeoptions); 
toc

subplot(3,1,1)
plot(t,xa(:,1:6))
title('Growth curves & Substrate/Product Concentrations')
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('S.bo Biomass','B.th Biomass','E.re Biomass','M.si Biomass', 'Cancer Biomass', 'Colon Biomass')
subplot(3,1,2)
plot(t,xa(:,7:12))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Glucose','Amino Acid Pool','Water','Oxygen','Phosphate','Ammonium')
subplot(3,1,3)
plot(t,xa(:,13:18))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Carbon dioxide','Propanoate','Butyrate','Methane','MFalpha2','Myrosinase')