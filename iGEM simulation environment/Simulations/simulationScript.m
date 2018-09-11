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

%% 1.0 Load & Check Genome Scale Models
% Please ensure that each model being loaded has had some degree of
% curation, and can be used for simulation (i.e. has objective function,
% set bounds, etc). This will be done in the following sections.
% To run these sections, go to the GEM folder

%% 1.1 Saccharomyces boulardii model

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

%% 1.2 Bacteroides thetaiotamicron model

modelBth = importModel('iBth801 v1.00.xml');
modelBth.id = 'Bth';
modelBth=setParam(modelBth,'obj','biomass_red',1);
modelBth=setParam(modelBth,'lb','glcIn',0); % make glucose uptake irreversible
modelBth=mergeCompartments(modelBth);
solBth=solveLP(modelBth,1);
printFluxes(modelBth,solBth.x);
%exportToExcelFormat(modelBth,'modelBth.xlsx');  

%% 1.3 Eubacterium rectale model - growth achieved

modelEre = importModel('iEre400 v1.00.xml');
modelEre.id = 'Ere';
modelEre=setParam(modelEre,'obj','r55',1);
modelEre=setParam(modelEre,'lb','glcIn',0); % make glucose uptake irreversible
modelEre=setParam(modelEre,'ub','r113',0); % make water uptake irreversible
modelEre=setParam(modelEre,'lb','r114',0); % make phosphate uptake irreversible
modelEre=mergeCompartments(modelEre);
solEre=solveLP(modelEre,1);
printFluxes(modelEre,solEre.x);
%exportToExcelFormat(modelEre,'modelEre.xlsx');

%% 1.4 Methanobrevibacter smithii model
% This model gave me a lot of trouble. It seems to contain many
% gaps/isolated subnetworks, as can been seen from running
% gapReport(modelMsi). To try to make this a functioning model I had to
% unconstrain all reactions and create a biomass equation without the
% metabolites that could not be obtained. To get methane production, I had
% to introduce methane as a bi-product of growth. This approach was 
% discarded, as it was later it was found that
% methane is produced at the same rate as aceteate is produced, therefore we
% define the methane production equal to the rate of acetate comsumption in
% the function f().

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

modelCancer = setParam(modelCancer,'lb',{'HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079'},-1); %allow for uptake/secretion of some mets
modelCancer = setParam(modelCancer,'ub',{'HMR_9034','HMR_9048','HMR_9047','HMR_9058','HMR_9072','HMR_9073','HMR_9074','HMR_9075','HMR_9076','HMR_9077','HMR_9078','HMR_9079','HumanGrowth','humanGrowthOut'},1); % allow for excretion/consumption of some mets
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

modelCancer = removeRxns(modelCancer,{'HMR_1592','HMR_1696','HMR_4271','HMR_0006','HMR_0007','HMR_0008','HMR_0015','HMR_0016','HMR_0017','HMR_1919','HMR_7568','HMR_7569','HMR_2141','HMR_4184','HMR_2585'});%these intercompartment transport rxns may become problematic after merging compartments
modelCancer = mergeCompartments(modelCancer);
solCancer=solveLP(modelCancer,1);
printFluxes(modelCancer,solCancer.x);
%exportToExcelFormat(modelCancer,'modelCancer.xlsx');

%% 1.6 Colon tissue model 

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
modelColon.mets(2462)={'proteinPool'};%for some reason addMets fails to do this

clear rxnToAdd % Make sure this structure is empty
rxnToAdd.rxns={'human_proteinPool','human_Maintenance'};
rxnToAdd.rxnNames={'human_proteinPool','human_Maintenance'};
rxnToAdd.equations={'0.0937 alanine[s] + 0.0507 arginine[s] + 0.04 asparagine[s] + 0.04 aspartate[s] + 0.0142 cysteine[s] + 0.1118 glutamine[s] + 0.1831 glycine[s] + 0.0198 histidine[s] + 0.0309 isoleucine[s] + 0.0664 leucine[s] + 0.0571 lysine[s] + 0.0156 methionine[s] + 0.029 phenylalanine[s] + 0.0853 proline[s] + 0.0491 serine[s] + 0.0402 threonine[s] + 0.0072 tryptophan[s] + 0.019 tyrosine[s] + 0.0471 valine[s] <=> proteinPool[s]','ATP[s] + 6 H2O[s] + proteinPool[s] => ADP[s] + Pi[s] + H+[s]'};
rxnToAdd.lb=[-1000,0];
rxnToAdd.ub=[1000,1000];
modelColon=addRxns(modelColon,rxnToAdd,3,'',true);

modelColon = setParam(modelColon,'obj','human_Maintenance',1);
modelColon = rmfield(modelColon,'rxnComps');% only need to do after merging compartments for exportToExcel() to work

modelColon = setParam(modelColon,'ub',{'HMR_9047'},0)%constrain exch rxns to be irreversible
modelColon = setParam(modelColon,'lb',{'HMR_9034','HMR_9048','HMR_9058','HMR_9073','HMR_9078','HMR_9079'},0)

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

modelColon=setParam(modelColon,'ub','HMR_9062',1000); %aspargine
modelColon=setParam(modelColon,'ub','HMR_9066',1000); %arginine
modelColon=setParam(modelColon,'ub','HMR_9067',1000); %glycine
modelColon=setParam(modelColon,'ub','HMR_9064',1000); %tyrosine
modelColon=setParam(modelColon,'ub','HMR_9065',1000); %cysteine

solColon = solveLP(modelColon,1);
printFluxes(modelColon,solColon.x);
%exportToExcelFormat(modelColon,'modelColon.xlsx');


%% 2.0 Initialize Models
% To run this section, go to the scripts folder

%define cell strucutre of models to be initialized
models = {modelSbo modelBth modelEre modelMsi modelCancer modelColon};

%inspect superModel structure, refer to initModels() function file for details
[superModel] = initModels(models); 

% contains summary of all exchange rxns
out = superModelFluxes(superModel);

%% 3.0 Define kinetic paramters, initial conditions & flow rates
% 

% Load excel file with kinetic parameters, initial concentrations, and flow
% rates
[num, txt, raw] = xlsread('parameters.xlsx');

% Define initial conditions
sbo0 = num(77);
bth0 = num(78);
ere0 = num(79);
msi0 = num(80);
can0 = num(81);
col0 = num(82);
glu0 = num(83);
wat0 = num(84);
oxy0 = num(85);
pho0 = num(86);
amm0 = num(87);
ace0 = num(88);
Q0 = num(89);
H0 = num(90);
K0 = num(91);
F0 = num(92);
V0 = num(93);
T0 = num(94);
W0 = num(95);
M0 = num(96);
L0 = num(97);
I0 = num(98);
car0 = num(99);
pro0 = num(100);
but0 = num(101);
met0 = num(102);
mfa0 = num(103);
myr0 = num(104);
p280 = num(105);

% Summary
initialConditions = [sbo0 bth0 ere0 msi0 can0 col0 glu0 wat0 oxy0 pho0 amm0 ace0 Q0 H0 K0 F0 V0 T0 W0 M0 L0 I0 car0 pro0 but0 met0 mfa0 myr0 p280];

%% 4.0 Simulations
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
x(8):Water (mmol/L)
x(9):Oxygen (mmol/L)
x(10):Phosphate (mmol/L)
x(11):Ammonium (mmol/L)
x(12):Acetate (mmol/L)

AMINO ACIDS
x(13):Glutamine (mmol/L)
x(14):Histidine (mmol/L)
x(15):Lysine (mmol/L)
x(16):Phenylalanine (mmol/L)
x(17):Valine (mmol/L)
x(18):Threonine (mmol/L)
x(19):Tryptophan (mmol/L)
x(20):Methionine (mmol/L)
x(21):Leucine (mmol/L)
x(22):Isoleucine (mmol/L)

PRODUCTS
x(23):Carbon dioxide (mmol/L)
x(24):Propanoate (mmol/L)
x(25):Butyrate (mmol/L)
x(26):Methane (mmol/L)
x(27):MFalpha2 (mmol/L)
x(28):Myrosinase (mmol/L)
x(29):P28 (mmol/L)

%}

odeoptions = odeset('NonNegative',1:29);
tic
[t,xa] = ode15s(@(t,x)f(t,x,superModel,num),[0 100], initialConditions ,odeoptions); 
toc

figure(1)
plot(t,xa(:,1:6))
title('Growth curves')
xlabel('Time (hours)'), ylabel('Concentration (g/L)')
legend('S.bo Biomass','B.th Biomass','E.re Biomass','M.si Biomass', 'Cancer Biomass', 'Colon Biomass')

figure(2)
title('Substrates')
plot(t,xa(:,7:12))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Glucose','Water','Oxygen','Phosphate','Ammonium','Acetate')

figure(3)
title('Amino acids')
plot(t,xa(:,13:22))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Glutamine','Histidine','Lysine','Phenylalanine','Valine','Threonine','Tryptophan','Methionine','Leucine','Isoleucine')

figure(4)
title('Products')
plot(t,xa(:,23:29))
xlabel('Time (hours)'), ylabel('Concentration (mmol/gDCW)')
legend('Carbon dioxide','Propanoate','Butyrate','Methane','MFalpha2','Myrosinase','P28')