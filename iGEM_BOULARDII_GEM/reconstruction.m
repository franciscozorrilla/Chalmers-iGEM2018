%% iGEM CHALMERS 2018 Saccharomyces boulardii GEM RECONSTRUCTION FROM HOMOLOGY
%  Author: Francisco Zorrilla

%% 1. INSTALL RAVEN & SET UP WORKING DIRECTORY

%{

Download and install the raven toolbox from 
https://github.com/SysBioChalmers/RAVEN

%}

checkInstallation
cd C:\Users\zorrilla\Desktop\iGEM_BOULARDII_GEM % Change this to your working directory

%% 2. OBTAIN PROTEIN FASTA FILES & PERFORM BLAST

%{

Download Saccharomyces boulardii protein fasta using the NCBI database link 
below. Gunzip the file and rename it sbo_prot.faa. Place it in the working
directory.
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/413/975/GCA_001413975.1_ASM141397v1/GCA_001413975.1_ASM141397v1_protein.faa.gz

Download Saccharomyces cerevisiae protein fasta using the NCBI database
link below. Gunzip the file and rename it sce_prot.fasta. Place it in the
working directory.
https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz

%}

% Obtain bi-directional BLAST structure using FASTA files
sboBlastStructure = getBlast('sbo','sbo_prot.faa',{'sce'},{'sce_prot.fasta'});
% Save the BLAST structure
save('sboBlastStructure.mat','sboBlastStructure');
% In the future, load BLAST structure instead of re-running BLAST function
load('sboBlastStructure.mat');

%% 3. OBTAIN HIGH QUALITY Saccharomyces cerevisiae GEM

%{

Download the yeast gem from https://github.com/SysBioChalmers/yeast-GEM
Place the yeastGEM.xml file in the working directory

%}

% Load Saccharomyces cerevisiae GEM
modelSce=importModel('yeastGEM.xml');
% Change Saccharomyces cerevisiae GEM ID to 'sce' for getModelFromHomology() function to work properly 
modelSce.id='sce';
% Set objective function to growth
modelSce=setParam(modelSce,'obj','r_4041',1);
% Solve for fluxes
sol=solveLP(modelSce,1);
% Check fluxes
printFluxes(modelSce,sol.x);
% Export to excel for inspection if necessary
%exportToExcelFormat(modelSce,'modelSce.xlsx');
clear sol

%% 4. GET Saccharomyces boulardii HOMOLOGY MODEL

% To avoid keeping unneccesary old genes, the models should not have
% "OR" relations in their grRules. To remove these, use expandModel()
% to separate reactions that are annotated with isoenzymes into
% multiple reactions. Warnings can be ignored.
modelSceExpanded=expandModel(modelSce);

% Obtain homology model using getModelFromHomology() function with BLAST
% structure
modelSbo=getModelFromHomology({modelSceExpanded},sboBlastStructure,'sbo',{},1,false,10^-20,100,35);

% Contract model to combine the reactions were seperated by expandModel() 
% above. Isoenzymes are now annotated to the same reaction, with 'OR' gene 
% relationships.
modelSbo=contractModel(modelSbo);

% Remove artifact reactions with _EXP_ in their name. These are reactions 
% catalyzed by proteins corresponding to genes that were present 
% in the template model (Sce) but not in the target model (Sbo), therefore
% should be removed.
modelSbo.rxns=regexprep(modelSbo.rxns,'_EXP_.','')

% Export model to excel for inspection
% exportToExcelFormat(modelSbo,'modelSbo.xlsx');

clear sboBlastStructure modelSceExpanded

%% 5. Saccharomyces boulardii GEM CURATION

% Need to add transport, exchange, and any non-gene-annotated-reactions
% necessary for growth
noGeneIdx=find(cellfun(@isempty,modelSce.grRules)); % Find rxns with no genes
rxnIdx=regexp(modelSce.rxnNames,'(transport)|(diffusion)|(carrier)|(shuttle)|(growth)|(biomass)|(lipid)|(protein)|(carbohydrate)|(exchange)'); % Find rxns with key words of interest in name
rxnIdx=find(~cellfun('isempty',rxnIdx)); % Find indices
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns=modelSce.rxns(rxnIdx); % Obtain reaction IDs
modelSbo=addRxnsGenesMets(modelSbo,modelSce,rxns,false,...
    'Modeling reaction required for growth, gene unknown',1); % Add reactions and metabolites to model

% Ensure objective function set to growth
modelSbo=setParam(modelSbo,'obj','r_4041',1);

% Solve for fluxes
sol=solveLP(modelSbo,1);
%{
sol = 

  struct with fields:

       x: []
       f: []
    stat: -1
     msg: 'The problem is infeasible'

%}

% Check fluxes
printFluxes(modelSbo,sol.x); %ERROR: GEM cannot grow, need to do gap filling

% Force model to push flux through growth reaction by setting the lower 
% bound to an arbitrary positive number
modelSbo = setParam(modelSbo,'lb','r_4041',0.01);
% Run fillGaps() function to fill gaps in the model based on growth
% constraint
[newConnected, cannotConnect, addedRxns, modelSbo, exitFlag]=fillGaps(modelSbo,modelSce,false,true);
% Reset lower bound of growth reaction to zero in gap filled model to test
% if model can now grow
modelSbo = setParam(modelSbo,'lb','r_4041',0);

% Ensure objective function set to growth
modelSbo=setParam(modelSbo,'obj','r_4041',1);
% Solve for fluxes
sol=solveLP(modelSbo,1);
% Check fluxes
printFluxes(modelSbo,sol.x); % Growth achieved

% Export model for inspection
%exportToExcelFormat(modelSbo,'modelSbo.xlsx');

clear addedRxns cannotConnect exitFlag newConnected noGeneIdx rxnIdx rxns sol

%% 6. ADD PROTEINS OF INTEREST (MFalpha2, Myrosinase, and P28) TO S.bo MODEL

clear rxnToAdd metsToAdd % Make sure these structures are empty

metsToAdd = {};
metsToAdd.mets = {'MFalpha2','Myrosinase','P28'};
metsToAdd.metNames = {'MFalpha2','Myrosinase','P28'};
metsToAdd.compartments = {'e';'e';'e'};
modelSbo=addMets(modelSbo,metsToAdd);

protStruct = {}; %create empty protein structure
protStruct.name = {'MFalpha2','Myrosinase','P28'}; %add protein names
protStruct.seq = {'MKFISTFLTFILAAVSVTASSDEDIAQVPAEAIIGYLDFGGDHDIAFLPFSNATASGLLFINTTIAEAAEKEQNTTLAKREAVADAWHWLNLRPGQPMYKREANADAWHWLQLKPGQPMY', ...
    'MKHLGLILAFLLALATCKADEEITCEENLPFKCSQPDRLNSSSFEKDFIFGVASSAYQACCLGRGLNVWDGFTHRYPNKSGPDHGNGDTTCDSFSYWQKDIDVLDELNATGYRFSIAWSRIIPRGKRSRGVNKDGINYYHGLIDGLIDKGITPFVTLFHWDLPQVLQDEYEGFLDPQIIHDFKHYANLCFQEFGHKVKNWLTINQLYTVPTRGYGAGSDAPGRCSPMVDPTCYAGNSSTEPYIVAHNQLLAHATVVDLYRKNYSIGPVMITRWFLPYNDTDPDSIAATERMKEFFLGWFMGPLTNGTYPQIMIDTVGERLPSFSPEESNLVKGSYDYLGLNYYVTQYAQPSPNPVHWANHTAMMDAGAKLTFRGNSDETKNSYYYPKGIYYVMDYFKTKYYNPLIYVTENGISTPGNETRDESMLHYKRIEYLCSHLCFLSKVIKEKHVNVKGYFAWSLGDNYEFDKGFTVRFGLSYIDWNNVTDRDLKLSGKWYQKFISPAIKNPLKKDFLRSSLTFEKNKKFEDA', ...
    'LSTAADMQGVVTDGMASGLDKDYLKPDD'}; %add amino acid sequences

for x = 1:length(protStruct.name) %create an equation for each protein using makeProteinEqn()
   protStruct.eqn{x} = makeProteinEqn(protStruct.name{x},protStruct.seq{x});
end

clear rxnToAdd metsToAdd % Make sure these structures are empty

rxnToAdd.rxns={'r_4066','r_4067','r_4068','r_4069','r_4070','r_4071','r_4072','r_4073','r_4074'}; % These numbers were determined by finiding the highest rxn ID in the excel model
rxnToAdd.rxnNames={'MFalpha2 formation','Myrosinase formation','P28 formation','MFalpha2 transport','Myrosinase transport','P28 transport','MFalpha2 exchange','Myrosinase exchange','P28 exchange'};
rxnToAdd.equations={protStruct.eqn{1}, protStruct.eqn{2}, protStruct.eqn{3}, 'MFalpha2[c] <=> MFalpha2[e]', 'Myrosinase[c] <=> Myrosinase[e]', 'P28[c] <=> P28[e]', 'MFalpha2[e] <=> ', 'Myrosinase[e] <=> ', 'P28[e] <=> '};
%rxnToAdd.eccodes= {'NA','NA','NA'};
%rxnToAdd.grRules= {'NA','NA','NA'};
rxnToAdd.subSystems={{'Alpha Pheromone'},{'Myrosinase'},{'P28'}, {'Alpha Pheromone'},{'Myrosinase'},{'P28'},{'Alpha Pheromone'},{'Myrosinase'},{'P28'}}; %NOTE: each subsystem name has to be its own cell for exportToExcel() to work
%rxnToAdd.rxnReferences = {'NA','NA','NA'};
%rxnToAdd.confidenceScores = {'NA','NA','NA'};
rxnToAdd.lb=[0,0,0,-1000,-1000,-1000,0,0,0];
rxnToAdd.ub=[1000,1000,1000,1000,1000,1000,1000,1000,1000];
modelSbo=addRxns(modelSbo,rxnToAdd,3,'',true);

exportToExcelFormat(modelSbo,'modelSbo.xlsx');

clear rxnToAdd protStruct x

%% 7. SIMULATIONS

% Check if proteins can feasibly be produced by changing objective function

modelSbo=setParam(modelSbo,'obj','r_4041',1); %biomass production
sol=solveLP(modelSbo,1);
printFluxes(modelSbo,sol.x);
%{
EXCHANGE FLUXES:
r_2111	(growth):	0.074486
r_1654	(ammonium exchange):	-0.41022
r_1672	(carbon dioxide exchange):	3.283
r_1714	(D-glucose exchange):	-1
r_1861	(iron(2+) exchange):	-7.4486e-08
r_1992	(oxygen exchange):	-3.1396
r_2005	(phosphate exchange):	-0.018879
r_2060	(sulphate exchange):	-0.0056902
%}

modelSbo=setParam(modelSbo,'obj','r_4072',1); %MFalpha2 production
sol=solveLP(modelSbo,1);
printFluxes(modelSbo,sol.x);
%{
EXCHANGE FLUXES:
r_1654	(ammonium exchange):	-0.76809
r_1672	(carbon dioxide exchange):	2.9775
r_1714	(D-glucose exchange):	-1
r_1992	(oxygen exchange):	-2.8279
r_2060	(sulphate exchange):	-0.014963
r_4072	(MFalpha2 exchange):	0.0049876
%}

modelSbo=setParam(modelSbo,'obj','r_4073',1); %Myrosinase production
sol=solveLP(modelSbo,1);
printFluxes(modelSbo,sol.x);
%{
EXCHANGE FLUXES:
r_1654	(ammonium exchange):	-0.7732
r_1672	(carbon dioxide exchange):	2.9947
r_1714	(D-glucose exchange):	-1
r_1992	(oxygen exchange):	-2.8541
r_2060	(sulphate exchange):	-0.022966
r_4073	(Myrosinase exchange):	0.0010936
%}

modelSbo=setParam(modelSbo,'obj','r_4074',1); %P28 production
sol=solveLP(modelSbo,1);
printFluxes(modelSbo,sol.x);
%{
EXCHANGE FLUXES:
r_1654	(ammonium exchange):	-0.76725
r_1672	(carbon dioxide exchange):	2.9805
r_1714	(D-glucose exchange):	-1
r_1992	(oxygen exchange):	-2.8444
r_2060	(sulphate exchange):	-0.0495
r_4074	(P28 exchange):	0.02475
%}

%% 8. FINAL CLEAN UP OF Saccharomyces boulardii MODEL

modelSbo = modelSbo; % for clarity rename final model as modelSbo
modelSbo = deleteUnusedGenes(modelSbo); % delete genes that are not tied to any reaction
modelSbo.id='Sbo';
modelSbo.description='Genome-scale model of Saccharomyces boulardii';
modelSbo.annotation.familyName='Zorrilla';
modelSbo.annotation.givenName='Francisco';
modelSbo.annotation.email='zorrilla@chalmers.se';
modelSbo.annotation.organization='Chalmers University of Technology';
modelSbo.annotation.notes='This genome scale model reconstruction was performed for modeling of the Chalmers iGEM 2018 team, please cite.';

exportToExcelFormat(modelSbo,'modelSbo.xlsx');