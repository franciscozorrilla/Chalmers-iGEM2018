function out=exportForGit(model,prefix,path,formats,masterFlag)
% exportForGit
%   Generates a directory structure and populates this with model files, ready
%   to be commited to a Git(Hub) maintained model repository. Writes the model
%   as SBML L3V1 FBCv2 (both XML and YAML), COBRA text, Matlab MAT-file
%   orthologies in KEGG
%
%   model               model structure in RAVEN format that should be exported
%   prefix              prefix for all filenames (opt, default 'model')
%   path                path where the directory structure should be generated
%                       and populated with all files (opt, default to current
%                       working directory)
%   formats             cell array of strings specifying in what file formats
%                       the model should be exported (opt, default to all
%                       formats as {'mat', 'txt', 'xlsx', 'xml', 'yml'})
%   masterFlag          logical, if true, function will error if RAVEN (and
%                       COBRA if detected) is/are not on the master branch.
%                       (opt, default false)
%
%   Usage: exportForGit(model,prefix,path,formats,masterFlag)
%
%   Benjamin J. Sanchez, 2018-08-03
%
if nargin<5
    masterFlag=false;
end
if nargin<4
    formats={'mat', 'txt', 'xlsx', 'xml', 'yml'};
end
if ischar(formats)
    formats={formats};
end
if any(~ismember(formats, {'mat', 'txt', 'xlsx', 'xml', 'yml'}))
    EM='Unknown file format defined. Only mat, txt, xlsx, xml and yml are allowed file formats.';
    error(EM)
end
if nargin<3
    path='.';
end
if nargin<2
    prefix='model';
end

%Get RAVEN version:
RAVENpath = which('ravenCobraWrapper.m');
slashPos  = getSlashPos(RAVENpath);
RAVENpath = RAVENpath(1:slashPos(end-1));
checkIfMaster(RAVENpath,'RAVEN',masterFlag)
try
    fid      = fopen([RAVENpath 'version.txt'],'r');
    RAVENver = fscanf(fid,'%s');
    fclose(fid);
catch
    RAVENver = 'unknown';
end

%Retrieve latest COBRA commit:
COBRApath = which('initCobraToolbox.m');
if ~isempty(COBRApath)
    slashPos  = getSlashPos(COBRApath);
    COBRApath = COBRApath(1:slashPos(end));
    checkIfMaster(COBRApath,'COBRA',masterFlag)
    currentPath = pwd;
    cd(COBRApath)
    try
        COBRAcommit = git('log -n 1 --format=%H');
    catch
        disp('COBRA is not fully installed (including Git wrapper)')
        COBRAcommit = 'unknown';
    end
    cd(currentPath)
else
    disp('COBRA version cannot be found')
end

%Retrieve libSBML version:
try % 5.17.0 and newer
    libSBMLver=OutputSBML;
    libSBMLver=libSBMLver.libSBML_version_string;
catch % before 5.17.0
    fid = fopen('tempModelForLibSBMLversion.xml','w+');
    fclose(fid);
    evalc('[~,~,libSBMLver]=TranslateSBML(''tempModelForLibSBMLversion.xml'',0,0)');
    libSBMLver=libSBMLver.libSBML_version_string;
    delete('tempModelForLibSBMLversion.xml');
end

% Make ModelFiles folder, no warnings if folder already exists
[~,~,~]=mkdir(fullfile(path,'ModelFiles'));
for i = 1:length(formats)
    [~,~,~]=mkdir(fullfile(path,'ModelFiles',formats{i}));
end

% Write TXT format
if ismember('txt', formats)
    fid=fopen(fullfile(path,'ModelFiles','txt',strcat(prefix,'.txt')),'w');
    eqns=constructEquations(model,model.rxns,false,false,false,true);
    eqns=strrep(eqns,' => ','  -> ');
    eqns=strrep(eqns,' <=> ','  <=> ');
    eqns=regexprep(eqns,'> $','>');
    grRules=regexprep(model.grRules,'\((?!\()','( ');
    grRules=regexprep(grRules,'(?<!\))\)',' )');
    fprintf(fid, 'Rxn name\tFormula\tGene-reaction association\tLB\tUB\tObjective\n');
    for i = 1:numel(model.rxns)
        fprintf(fid, '%s\t', model.rxns{i});
        fprintf(fid, '%s \t', eqns{i});
        fprintf(fid, '%s\t', grRules{i});
        fprintf(fid, '%6.2f\t%6.2f\t%6.2f\n', model.lb(i), model.ub(i), model.c(i));
    end
    fclose(fid);
end

% Write YML format
if ismember('yml', formats)
    writeYaml(model,fullfile(path,'ModelFiles','yml',strcat(prefix,'.yml')));
end

% Write MAT format
if ismember('mat', formats)
    save(fullfile(path,'ModelFiles','mat',strcat(prefix,'.mat')),'model');
end

% Write XLSX format
if ismember('xlsx', formats)
    exportToExcelFormat(model,fullfile(path,'ModelFiles','xlsx',strcat(prefix,'.xlsx')));
end

% Write XML format
if ismember('xml', formats)
    exportModel(model,fullfile(path,'ModelFiles','xml',strcat(prefix,'.xml')));
end

%Save file with versions:
fid = fopen(fullfile(path,'ModelFiles','dependencies.txt'),'wt');
fprintf(fid,['MATLAB\t' version '\n']);
fprintf(fid,['libSBML\t' libSBMLver '\n']);
fprintf(fid,['RAVEN_toolbox\t' RAVENver '\n']);
if ~isempty(COBRApath)
    fprintf(fid,['COBRA_toolbox\tcommit ' COBRAcommit(1:7) '\n']);
end
if isfield(model,'modelVersion')
    fields = fieldnames(model.modelVersion);
    for i = 1:length(fields)
        value = model.modelVersion.(fields{i});
        fprintf(fid,[fields{i} '\t' num2str(value) '\n']);
    end
end
fclose(fid);
end

function slashPos = getSlashPos(path)
slashPos = strfind(path,'\');       %Windows
if isempty(slashPos)
    slashPos = strfind(path,'/');   %MAC/Linux
end
end

function checkIfMaster(toolboxPath,toolbox,masterFlag)
if masterFlag
    currentPath = pwd;
    cd(toolboxPath)
    currentBranch = git('rev-parse --abbrev-ref HEAD');
    cd(currentPath)
    if ~strcmp(currentBranch,'master')
        error(['ERROR: ' toolbox ' not in master. Check-out the master branch of ' toolbox ' before submitting model for Git.'])
    end
end
end
