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
% set bounds, etc)

% Saccharomyces boulardii model - growth achieved
load('modelSbo.mat'); 
modelSbo.id = 'Sbo';
exportToExcelFormat(modelSbo,'modelSbo.xlsx');
modelSbo=setParam(modelSbo,'obj','r_4041',1);
solSbo=solveLP(modelSbo,1);
printFluxes(modelSbo,solSbo.x);

% Bacteroides thetaiotamicron model - growth achieved
modelBth = importModel('iBth801 v1.00.xml');
modelBth.id = 'Bth';
exportToExcelFormat(modelBth,'modelBth.xlsx');  
modelBth=setParam(modelBth,'obj','bmOut',1);
solBth=solveLP(modelBth,1);
printFluxes(modelBth,solBth.x);

% Eubacterium rectale model - growth achieved
modelEre = importModel('iEre400 v1.00.xml');
modelEre.id = 'Ere';
exportToExcelFormat(modelEre,'modelEre.xlsx');
modelEre=setParam(modelEre,'obj','bmOut',1);
solEre=solveLP(modelEre,1);
printFluxes(modelEre,solEre.x);

% Methanobrevibacter smithii model - cannot carry flux through biomass rxn
modelMsi = importModel('iMsi385.xml');
modelMsi.id = 'Msi';
exportToExcelFormat(modelMsi,'modelMsi.xlsx');
modelMsi=setParam(modelMsi,'obj','r449',1);
modelMsi=setParam(modelMsi,'lb',{'NH4IN','piIN'},-1);
modelMsi=setParam(modelMsi,'lb','r58',-5);
modelMsi=setParam(modelMsi,'ub','r58',5);
solMsi=solveLP(modelMsi,1);
printFluxes(modelMsi,solMsi.x);

% Colon tissue model - no objective function
modelColon = importModel('colon.xml');
modelColon.id = 'Col';
exportToExcelFormat(modelColon,'modelColon.xlsx');

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
exportToExcelFormat(modelCancer,'modelCancer.xlsx');
modelCancer=setParam(modelCancer,'obj','humanGrowthOut',1);
modelCancer.ub(8177:8192) = 5; % CONSTRAIN INNER FLUXES, MAYBE DO WITH FUNCTIONS?
solCancer=solveLP(modelCancer,1);
printFluxes(modelCancer,solCancer.x);
 

%% 2. Initialize Models & Environment

models = {modelSbo modelBth modelEre modelCancer}; %define cell strucutre of models to be initialized

[superModel] = initModels(models); %inspect superModel structure, refer to initModels() function file for details

out = superModelFluxes(superModel); % contains summary of all exchange rxns

% NEXT NEED TO CREATE STRUCTURE FOR ENVIRONMENTAL CONCENTRATIONS 
