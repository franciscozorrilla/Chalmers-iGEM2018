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

cd('C:\Users\zorrilla\Desktop\COM-dFBA\GEMs')

%To save time loading the models simply load a saved workspace in the
%/Scripts folder
cd('C:\Users\zorrilla\Desktop\COM-dFBA\Scripts')
load('cleanWorkspace.mat')

%% 1.0 Load & Check Genome Scale Models
% Please ensure that each model being loaded has had some degree of
% curation, and can be used for simulation (i.e. has objective function,
% set bounds, etc). This will be done in the following sections.
% To run these sections, go to the GEM folder. Make sure that the exchange
% reactions for each state variable, particularly uptake metabololites
% (eg x(7)=glucose), has one of its bounds equal to zero. For example, the 
% exchange reaction 'r_1714' 'D-glucose[s] <=> 'should have upper bound set
% to zero, since it will never secrete glucose. This is how the dFBA loop
% will determine which bound it needs to modify.


%% 1.1 Saccharomyces boulardii model

load('modelSbo.mat'); 
modelSbo.id = 'Sbo';
modelSbo=setParam(modelSbo,'obj','r_4041',1); %set objective function

modelSbo=setParam(modelSbo,'lb','r_1891',-1000); %glutamine
modelSbo=setParam(modelSbo,'lb','r_1893',-1000); %histidine
modelSbo=setParam(modelSbo,'lb','r_1900',-1000); %lysine
modelSbo=setParam(modelSbo,'lb','r_1903',-1000); %phenylalanine
modelSbo=setParam(modelSbo,'lb','r_1914',-1000); %valine
modelSbo=setParam(modelSbo,'lb','r_1911',-1000); %threonine
modelSbo=setParam(modelSbo,'lb','r_1912',-1000); %tryptophan
modelSbo=setParam(modelSbo,'lb','r_1902',-1000); %methionine
modelSbo=setParam(modelSbo,'lb','r_1899',-1000); %leucine
modelSbo=setParam(modelSbo,'lb','r_1897',-1000); %isoleucine

modelSbo=setParam(modelSbo,'eq',{'r_1904','r_1906','r_1913','r_1889','r_1879','r_1880','r_1881','r_1883'},0); %block non essential amino acid uptake
modelSbo=setParam(modelSbo,'ub',{'r_1891','r_1893','r_1900','r_1903','r_1914','r_1911','r_1912','r_1902','r_1899','r_1897','r_1810'},0); % make AA-exchange & glucose irreversible

modelSbo=setParam(modelSbo,'lb','r_4066',0.001); %set mfalpha2 production to minimum of 0.001
modelSbo=setParam(modelSbo,'lb','r_4067',0.0001); %set Myrosinase production to minimum of 0.0001

modelSbo=mergeCompartments(modelSbo);%merge compartments to make simulations run faster, this will be done for every model

modelSbo.metNames(223) = {'carbon dioxide'}; % CO2 carbon dioxide for consistency with other models
modelSbo.metNames(532) = {'oxygen'}; % rename O2 as oxygen for consistency with other models
modelSbo.metNames(354) = {'water'}; % rename H2O as water for consistency
modelSbo.metNames(1028) = {'biomass S.bo'}; % rename biomass as biomass S.bo

modelSbo=setParam(modelSbo,'ub','r_1654',0); % make ammonium uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_1714',0); % make glucose uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_1992',0); % make oxygen uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_2005',0); % make phosphate uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_2060',0); % make sulphate uptake irreversible
modelSbo=setParam(modelSbo,'ub','r_2100',0); % make water uptake irreveersible

solSbo=solveLP(modelSbo,1);
printFluxes(modelSbo,solSbo.x);

%Export model to excel format for inspection
%exportToExcelFormat(modelSbo,'modelSbo.xlsx');

%% 1.2 Bacteroides thetaiotamicron model curation

modelBth = importModel('iBth801 v1.00.xml');
modelBth.id = 'Bth';
modelBth=setParam(modelBth,'obj','biomass_red',1);
modelBth=setParam(modelBth,'lb','glcIn',0); % make glucose uptake irreversible
modelBth=mergeCompartments(modelBth);

[modelBth, addedRxnsBth]=addExchangeRxns(modelBth,'in',{'m358','m404','m450','m598','m754','m683','m418','m441','m468','m714'});%add exch rxns for glutamine, histidine, lysine, phenylalanine, valine, threonine, isoleucine, leucine, methionine and tryptophan

modelBth.metNames(192) = {'carbon dioxide'}; % CO2 carbon dioxide for consistency with other models
modelBth.metNames(321) = {'water'}; % rename H2O as water for consistency
modelBth.metNames(613) = {'biomass B.th'}; % rename biomass as biomass B.th

solBth=solveLP(modelBth,1);
printFluxes(modelBth,solBth.x);

%exportToExcelFormat(modelBth,'modelBth.xlsx');  

%% 1.3 Eubacterium rectale model curation

modelEre = importModel('iEre400 v1.00.xml');
modelEre.id = 'Ere';
modelEre=setParam(modelEre,'obj','r55',1);
modelEre=setParam(modelEre,'lb','glcIn',0); % make glucose uptake irreversible
modelEre=setParam(modelEre,'ub','r113',0); % make water uptake irreversible
modelEre=setParam(modelEre,'lb','r114',0); % make phosphate uptake irreversible
modelEre=mergeCompartments(modelEre);

[modelEre, addedRxnsEre]=addExchangeRxns(modelEre,'in',{'m218','m296','m171','m195','m208','m237','m318','m279','m200','m207','m271'});%add exch rxns for methionine and tryptophan, glutamine, histidine, lysine, phenylalanine, valine, threonine, isoleucine, leucine, succinate

modelEre.metNames(110) = {'carbon dioxide'}; % CO2 carbon dioxide for consistency with other models
modelEre.metNames(182) = {'water'}; % rename H2O as water for consistency
modelEre.metNames(327) = {'biomass E.re'}; % rename biomass as biomass E.re

solEre=solveLP(modelEre,1);
printFluxes(modelEre,solEre.x);

%exportToExcelFormat(modelEre,'modelEre.xlsx');

%% 1.4 Methanobrevibacter smithii model curation
% This model gave me a lot of trouble. It seems to contain many
% gaps/isolated subnetworks, as can been seen from running
% gapReport(modelMsi). To try to make this a functioning model I had to
% unconstrain all reactions and create a biomass equation without the
% metabolites that could not be obtained.

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
rxnToAdd.equations={'0.5621 L-Alanine[s] + 0.3237 L-Arginine[s] + 0.2638 L-Asparagine[s] + 0.2638 L-Aspartate[s] + 50.2 ATP[s] + 0.1002 L-Cysteine[s] + 0.0331 dATP[s] + 0.0215 dGTP[s] + 0.288 L-Glutamine[s] + 0.288 L-Glutamate[s] + 0.6704 Glycine[s] + 0.2222 GTP[s] + 50.2 H2O[s] + 0.3179 L-Isoleucine[s] + 0.493 L-Leucine[s] + 0.3755 L-Lysine[s] + 0.1682 L-Methionine[s] + 0.2027 L-Phenylalanine[s] + 0.2419 L-Proline[s] + 0.2361 L-Serine[s] + 0.2776 L-Threonine[s] + 0.1509 L-Tyrosine[s] + 0.4631 L-Valine[s] => 50.2 ADP[s] + 50.2 Phosphate[s] + Biomass_msi[s]'};
rxnToAdd.lb=[0];
rxnToAdd.ub=[1000];
modelMsi=addRxns(modelMsi,rxnToAdd,3,'',true);
% Note: M smithii produces methane gas, but the methane synthesis pathway
% seems to be disconnected from the biomass pathway.

modelMsi = setParam(modelMsi,'lb','r225',0); % Need to constrain this reaction to be irreversible so that methane is excreted
msiRxns= getExchangeRxns(modelMsi,'in'); % get exchange fluxes going in
modelMsi = setParam(modelMsi, 'eq', msiRxns, 0); %constrain fluxes going in

modelMsi = setParam(modelMsi,'ub',{'actIn','NH4IN','piIN','CO2In'},1);
modelMsi = setParam(modelMsi,'ub',{'H2Suptake','H2In','H2ObTran'},1000);
modelMsi = setParam(modelMsi,'lb',{'H2ObTran'},-1001);
modelMsi = setParam(modelMsi,'lb',{'MethOut','bmOut_msi'},0);

modelMsi = setParam(modelMsi,'obj','biomassMsi',1);
[modelMsi, addedRxnsMsi]=addExchangeRxns(modelMsi,'in',{'m305','m339','m370','m417','m512','m470','m496','m380','m367','m354','m582'}); %add AA & formate exchange

modelMsi.metNames(184) = {'carbon dioxide'}; % CO2 carbon dioxide for consistency with other models
modelMsi.metNames(301) = {'water'}; % rename H2O as water for consistency
modelMsi.metNames(476) = {'biomass M.si'}; % rename biomass as biomass M.si

solMsi=solveLP(modelMsi,1);
printFluxes(modelMsi,solMsi.x); %now model can grow and produce methane at rate of 0.6813 mmol/gDCW*hr

%exportToExcelFormat(modelMsi,'modelMsi.xlsx');

%% 1.5 Cancer model

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
%[modelCancer, addedRxns]=addExchangeRxns(modelCancer,'in','m90000c'); % add protein pool exchange rxn
canRxns= getExchangeRxns(modelCancer,'in'); %get exchange rxns going in
modelCancer = setParam(modelCancer, 'eq', canRxns, 0); %constrain all exchange rxns going in 

modelCancer = setParam(modelCancer,'lb',{'HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077'},-1); %allow for uptake/secretion of some mets
modelCancer = setParam(modelCancer,'ub',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HumanGrowth','humanGrowthOut'},1); % allow for excretion/consumption of some mets
modelCancer = setParam(modelCancer,'lb',{'HumanGrowth','humanGrowthOut','HMR_9047','HMR_9048','HMR_9058','HMR_9034','HMR_9072'},0);%constrain growth related rxns
modelCancer = setParam(modelCancer,'ub',{'HMR_9073','HMR_9078','HMR_9079'},0);

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

modelCancer=setParam(modelCancer,'ub','HMR_9808',1000); %propanoate uptake
modelCancer=setParam(modelCancer,'ub','HMR_9809',1000); %butyrate uptake
modelCancer=setParam(modelCancer,'ub','HMR_9086',1000); %acetate uptake

modelCancer = removeReactions(modelCancer,{'HMR_1592','HMR_1696','HMR_4271','HMR_0006','HMR_0007','HMR_0008','HMR_0015','HMR_0016','HMR_0017','HMR_1919','HMR_7568','HMR_7569','HMR_2141','HMR_4184','HMR_2585','HMR_9391','HMR_9372'});%these intercompartment transport rxns may become problematic after merging compartments
modelCancer = mergeCompartments(modelCancer);

modelCancer.metNames(1594) = {'carbon dioxide'}; % rename CO2 carbon dioxide for consistency with other models
modelCancer.metNames(2021) = {'water'}; % rename H2O as water for consistency
modelCancer.metNames(2607) = {'oxygen'}; % rename O2 as oxygen for consistency with other models
modelCancer.metNames(2555) = {'ammonium'}; % rename NH3 as ammonium for consistency with other models
modelCancer.metNames(2728) = {'phosphate'}; % rename Pi as phosphate for consistency with other models
modelCancer.metNames(3137) = {'biomass Can'}; % rename biomass as biomass Can

solCancer=solveLP(modelCancer,1);
printFluxes(modelCancer,solCancer.x);

%exportToExcelFormat(modelCancer,'modelCancer.xlsx');

%% 1.6 Colon tissue model 

modelColon = importModel('colon.xml');
modelColon.id = 'Col';
modelColon = removeReactions(modelColon,{'HMR_1696','HMR_1849','HMR_0006','HMR_0007','HMR_0008','HMR_0015','HMR_1919','HMR_7568','HMR_7569','HMR_2585'}); %these intercompartment transport rxns may become problematic after merging compartments
colRxns= getExchangeRxns(modelColon,'in'); % get exchange fluxes going in
modelColon = setParam(modelColon, 'eq', colRxns, 0); %constrain fluxes going in
modelColon = setParam(modelColon,'lb',{'HMR_9047'},-1000);
modelColon = setParam(modelColon,'ub',{'HMR_9034','HMR_9048','HMR_9072','HMR_9073'},1);
modelColon = mergeCompartments(modelColon);
modelColon.genes = {''}; %insert dummy field to avoid errors with addRxns()

clear metsToAdd % Make sure this structure is empty
metsToAdd = {}; %Add proteinPool pseudo metabolite
metsToAdd.mets = {'proteinPool'};
metsToAdd.metNames = {'proteinPool'};
metsToAdd.compartments = {'s'};
modelColon=addMets(modelColon,metsToAdd); 
modelColon.mets(2462)={'proteinPool'};%for some reason addMets fails to do this

clear rxnToAdd % Make sure this structure is empty
rxnToAdd.rxns={'human_proteinPool','human_Maintenance'};
rxnToAdd.rxnNames={'human_proteinPool','human_Maintenance'};
rxnToAdd.equations={'0.0937 alanine[s] + 0.0507 arginine[s] + 0.04 asparagine[s] + 0.04 aspartate[s] + 0.0142 cysteine[s] + 0.1118 glutamine[s] + 0.1831 glycine[s] + 0.0198 histidine[s] + 0.0309 isoleucine[s] + 0.0664 leucine[s] + 0.0571 lysine[s] + 0.0156 methionine[s] + 0.029 phenylalanine[s] + 0.0853 proline[s] + 0.0491 serine[s] + 0.0402 threonine[s] + 0.0072 tryptophan[s] + 0.019 tyrosine[s] + 0.0471 valine[s] <=> proteinPool[s]','ATP[s] + 6 H2O[s] + proteinPool[s] => ADP[s] + Pi[s] + H+[s]'};
rxnToAdd.lb=[-1000,0];
rxnToAdd.ub=[1000,1000];
modelColon=addRxns(modelColon,rxnToAdd,3,'',true);
modelColon=setParam(modelColon,'ub','HMR_9047',0); %make sure water uptake is irreversible

modelColon = setParam(modelColon,'obj','human_Maintenance',1);
modelColon = rmfield(modelColon,'rxnComps');% only need to do after merging compartments for exportToExcel() to work

modelColon=setParam(modelColon,'ub','HMR_9063',1000); %glutamine
modelColon=setParam(modelColon,'ub','HMR_9038',1000); %histidine
modelColon=setParam(modelColon,'ub','HMR_9041',1000); %lysine
modelColon=setParam(modelColon,'ub','HMR_9043',1000); %phenylalanine
modelColon=setParam(modelColon,'ub','HMR_9046',1000); %valine
modelColon=setParam(modelColon,'ub','HMR_9044',1000); %threonine
modelColon=setParam(modelColon,'ub','HMR_9045',1000); %tryptophan
modelColon=setParam(modelColon,'ub','HMR_9042',1000); %methionine
modelColon=setParam(modelColon,'ub','HMR_9040',1000); %leucine
modelColon=setParam(modelColon,'ub','HMR_9039',1000); %isoleucine

modelColon=setParam(modelColon,'ub','HMR_9062',1000); %aspargine NON ESSENTIAL, THESE AAs ARE NOT IN MASS BALANCES BUT ARE REQUIRED FOR COLON MAINTENANCE. SHOULD NOT AFFECT dFBA
modelColon=setParam(modelColon,'ub','HMR_9066',1000); %arginine
modelColon=setParam(modelColon,'ub','HMR_9067',1000); %glycine
modelColon=setParam(modelColon,'ub','HMR_9064',1000); %tyrosine
modelColon=setParam(modelColon,'ub','HMR_9065',1000); %cysteine

modelColon=setParam(modelColon,'ub','HMR_9808',1000); %propanoate uptake
modelColon=setParam(modelColon,'ub','HMR_9809',1000); %butyrate uptake
modelColon=setParam(modelColon,'ub','HMR_9086',1000); %acetate uptake

modelColon.metNames(1281) = {'carbon dioxide'}; %rename CO2 as carbon dioxide for consistency with other models
modelColon.metNames(1614) = {'water'}; % rename H2O as water for consistency
modelColon.metNames(2075) = {'oxygen'}; % rename O2 as oxygen for consistency with other models
modelColon.metNames(2040) = {'ammonium'}; % rename NH3 as ammonium for consistency with other models
modelColon.metNames(2157) = {'phosphate'}; % rename Pi as phosphate for consistency with other models

solColon = solveLP(modelColon,1);
printFluxes(modelColon,solColon.x);

%exportToExcelFormat(modelColon,'modelColon.xlsx');

%%

clear addedRxnsBth addedRxnsEre addedRxnsMsi ans canRxns colRxns metsToAdd msiRxns q rxnToAdd

%% 2.0 Initialize Models
% To run this section, go to the scripts folder

cd('C:\Users\zorrilla\Desktop\COM-dFBA\Scripts')

%Define cell strucutre of models to be initialized
models = {modelSbo modelBth modelEre modelMsi modelCancer modelColon};

%Define parameter that determines if all fluxes should be saved/printed, or
%just the exchange fluxes. Note: Do not change whichFluxes value in between
%running initModels and superModelFluxes
whichFluxes = 1; %To get only exchange fluxes
%whichFluxes = 2; %To get all fluxes. Mostly for the purpose of
%debugging/troubleshooting. 

%Inspect superModel structure, refer to initModels() function file for details
[superModel] = initModels(models,whichFluxes); % If errors are thrown it is likely due to models not being solvable. Make sure there is a feasible solution before initializing models

%Contains summary of all exchange rxns. Look through upper and lower bounds
%of reactions to make sure that uptake reactions of metabolites in metNames
%have one bound equal to zero (or at least make sure that they are not 
%equal in magnitude).
fluxes = superModelFluxes(superModel,whichFluxes);

%Create cell array with metabolite names in order to match them in each
%model, a mass balance equation will be created for each of these
%metabolites during the dFBA simulation. To add new metabolites to the
%model, simply add to the end of the list and check that the correct
%exchange reactions have been "matched". Also remember to add Vmax/Ks data
%for these metabolites if they are uptaken by models in the paramTables 
%excel file. Best to add new metabolites at the end so that the excel file 
%with parameters does not have to be extensively reformatted.
metNames = {'biomass S.bo',
'biomass B.th',
'biomass E.re',
'biomass M.si',
'biomass Can',
'biomass Col',
'Glucose',
'Phosphate',
'Ammonium',
'Glutamine',
'Histidine',
'Lysine',
'Phenylalanine',
'Valine',
'Threonine',
'Tryptophan',
'Methionine',
'Leucine',
'Isoleucine',
'Carbon dioxide',
'Acetate',
'Propanoate',
'Butyrate',
'Succinate',
'Formate',
'Ethanol',
'Methane',
'MFalpha2',
'Myrosinase',
'P28'};

%Please check these structures before proceeding to simulations.
%metEquations should be of dimensions (#of models)x(#mets in metNames).
%Ensure that metabolites are correctly matched across the different
%organisms/columns. It is possible that an incorrect exchange reaction was
%matched, and this may or may not have significant effects on subsequent
%simulations. metIdx contains the reaction index of the corresponding
%exchange reaction shown in metEquations, while blocked contains the
%reaction index of an arbitrary exchange reaction with both bounds set to
%zero. It would be good to double check that this reaction indeed is not
%used and stays at zero flux.
[metEquations, metIdx, bounds, blocked] = checkFluxes(superModel,metNames);


%% 3.0 Initialize Environment & Kinetic Parameters
% The kinetic paramters, initial conditions, flow rates, and other
% adjustable parameters are defined in the excel sheet 'parameters.xlsx'.
% To adjust parameters, change the values in the excel sheet and then 
% re-run this section to load the new values into matlab.

% Load excel file with kinetic parameters, initial concentrations, and flow
% rates

params = struct;
params.Vmax = xlsread('paramTables.xlsx',1);
params.Ks = xlsread('paramTables.xlsx',2); 
params.mets = xlsread('paramTables.xlsx',3);
params.other = xlsread('paramTables.xlsx',4);
params.metIdx = metIdx;
params.blocked = blocked;

% Initial conditions definition
initialConditions = transpose(params.mets(:,1));

% Run one test FBA using initial conditions and kinetic parameters 
[superModelConstrained,fluxesConstrained,boundsConstrained,lbubConstrained] = updateFluxes(superModel,params);

%% 4.0 Simulation
% Run this section to perform community simulations. To remove any of the
% models from the community, simply set the initial biomass and inflow of 
% an organism to zero (e.g. ere0 = 0, eref = 0). Make sure to re-run
% section 3.0 in order to update the adjusted parameters.

%{
curves look better when models have their respective compartments merged. 
when I do it without merging then Bth outgrows the shit out of everything 
else. also simulations take waaaaay longer if comps not merged, probably 
because doing so greatly simplifies the cancer and colon models. ode113 
seems to work the best for simulations, ode15s also works well sometimes.

BIOMASS
x(1):S.bo biomass (g/L)
x(2):B.th biomass (g/L)
x(3):E.re biomass (g/L)
x(4):M.si biomass (g/L)
x(5):Cancer biomass (g/L)
x(6):Colon biomass (g/L)

SUBSTRATES
x(7):Glucose (mmol/L)
x(8):Phosphate (mmol/L)
x(9):Ammonium (mmol/L)

AMINO ACIDS
x(10):Glutamine (mmol/L)
x(11):Histidine (mmol/L)
x(12):Lysine (mmol/L)
x(13):Phenylalanine (mmol/L)
x(14):Valine (mmol/L)
x(15):Threonine (mmol/L)
x(16):Tryptophan (mmol/L)
x(17):Methionine (mmol/L)
x(18):Leucine (mmol/L)
x(19):Isoleucine (mmol/L)

PRODUCTS
x(20):Carbon dioxide (mmol/L)
x(21):Acetate (mmol/L)
x(22):Propanoate (mmol/L)
x(23):Butyrate (mmol/L)
x(24):Succinate (mmol/L)
x(25):Formate (mmol/L)
x(26):Ethanol (mmol/L)
x(27):Methane (mmol/L)
x(28):MFalpha2 (mmol/L)
x(29):Myrosinase (mmol/L)
x(30):P28 (mmol/L)

%}

odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'NonNegative',1:length(metNames),'Stats','on','InitialStep',1e-3);

tic
[t,xa] = ode15s(@(t,x)f(t,x,superModel,params),[0 500], initialConditions ,odeoptions); 
toc

%% 5.0 Visualization of Results

figure(1)
plot(t,xa(:,1:6))
title('Growth curves')
xlabel('Time (hours)'), ylabel('Concentration (g/L)')
legend('S.bo Biomass','B.th Biomass','E.re Biomass','M.si Biomass', 'Cancer Biomass', 'Colon Biomass')

figure(2)
title('Substrates')
plot(t,xa(:,7:9))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Glucose','Phosphate','Ammonium')

figure(3)
title('Amino acids')
plot(t,xa(:,10:19))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Glutamine','Histidine','Lysine','Phenylalanine','Valine','Threonine','Tryptophan','Methionine','Leucine','Isoleucine')

figure(4)
title('Products')
plot(t,xa(:,20:30))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Carbon dioxide','Acetate','Propanoate','Butyrate','Succinate','Ethanol','Methane','MFalpha2','Myrosinase','P28')

%Plot for simulation with only Sbo and Cancer
figure(5)
subplot(2,1,1)
plot(t,xa(:,[1 5]))
legend('S.bo Biomass','Cancer Biomass')
xlabel('Time (hours)'), ylabel('Concentration (g/L)')
subplot(2,1,2)
plot(t,xa(:,[28 29]))
legend('MFalpha2','Myrosinase')
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
