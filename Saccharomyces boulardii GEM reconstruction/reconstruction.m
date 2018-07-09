%% iGEM CHALMERS 2018 Saccharomyces boulardii GEM RECONSTRUCTION FROM HOMOLOGY
%  Author: Francisco Zorrilla

%% 1. INSTALL RAVEN & SET UP WORKING DIRECTORY

%{

Download and install the raven toolbox from 
https://github.com/SysBioChalmers/RAVEN

%}

checkInstallation
cd C:\Users\zorrilla\Desktop\iGEM_BOULARDII_GEM

%% 2. OBTAIN PROTEIN FASTA FILES

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
exportToExcelFormat(modelSce,'modelSce.xlsx');

%% 4. GET Saccharomyces boulardii HOMOLOGY MODEL

% Obtain bi-directional BLAST structure using FASTA files
sboBlastStructure = getBlast('sbo','sbo_prot.faa',{'sce'},{'sce_prot.fasta'});
% Save the BLAST structure
save('sboBlastStructure.mat','sboBlastStructure');
% In the future, load BLAST structure instead of re-running BLAST function
load('sboBlastStructure.mat');

% To avoid keeping unneccesary old genes, the models should not have
% "OR" relations in their grRules. To remove these, use expandModel()
% to separate reactions that are annotated with isoenzymes into
% multiple reactions. Warnings will be shown.
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
exportToExcelFormat(modelSbo,'modelSbo.xlsx');

%% 5. Saccharomyces boulardii GEM CURATION

% Need to add transport, exchange, and any non-gene-annotated-reactions
% necessary for growth
noGeneIdx=find(cellfun(@isempty,modelSce.grRules)); % Find rxns with no genes
rxnIdx=regexp(modelSce.rxnNames,'(transport)|(diffusion)|(carrier)|(shuttle)|(growth)|(biomass)|(lipid)|(protein)|(carbohydrate)'); % Find rxns with key words of interest in name
rxnIdx=find(~cellfun('isempty',rxnIdx)); % Find indices
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns=modelSce.rxns(rxnIdx); % Obtain reaction IDs
modelSboTest=addRxnsGenesMets(modelSbo,modelSce,rxns,false,...
    'Modeling reaction required for growth, gene unknown',1); % Add reactions and metabolites to model

% Export model for inspection
exportToExcelFormat(modelSboTest,'modelSboTest.xlsx');

% Ensure objective function set to growth
modelSboTest=setParam(modelSboTest,'obj','r_4041',1);
% Solve for fluxes
sol=solveLP(modelSboTest,1);
% Check fluxes
printFluxes(modelSboTest,sol.x); %GEM cannot grow, need to do gap filling

% Force model to push flux through growth reaction by setting the lower 
% bound to an arbitrary positive number
modelSboTest = setParam(modelSboTest,'lb','r_4041',0.01);
% Run fillGaps() function to fill gaps in the model based on growth
% constraint
[newConnected, cannotConnect, addedRxns, modelSboGF, exitFlag]=fillGaps(modelSboTest,modelSce,false,true);
% Reset lower bound of growth reaction to zero in gap filled model to test
% if model can now grow
modelSboGF = setParam(modelSboGF,'lb','r_4041',0);

% Ensure objective function set to growth
modelSboGF=setParam(modelSboGF,'obj','r_4041',1);
% Solve for fluxes
sol=solveLP(modelSboGF,1);
% Check fluxes
printFluxes(modelSboGF,sol.x); % Growth achieved

%% 6. FINAL CLEAN UP OF Saccharomyces boulardii MODEL

modelSbo = modelSboGF; % for clarity rename final model as modelSbo
modelSbo = deleteUnusedGenes(modelSbo); % delete genes that are not tied to any reaction
modelSbo.id='Sbo';
modelSbo.description='Genome-scale model of Saccharomyces boulardii';
modelSbo.annotation.familyName='Zorrilla';
modelSbo.annotation.givenName='Francisco';
modelSbo.annotation.email='zorrilla@chalmers.se';
modelSbo.annotation.organization='Chalmers University of Technology';
modelSbo.annotation.notes='This genome scale model reconstruction was performed for modeling of the Chalmers iGEM 2018 team, please cite.';

exportToExcelFormat(modelSbo,'modelSbo.xlsx');

%% 7. MODEL EVALUATION

