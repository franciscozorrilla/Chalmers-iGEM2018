%% These are loops, code snippets, etc. that could be of use

%% find extracellular compartment in each model, rename consistently in
%  superModel.SubModel structure. IS THIS REALLY NECESSARY??

exComp = {};
for p=1:length(models) 
    for o=1:length(models{p}.compNames)
        if strcmpi(models{p}.compNames(o),'extracellular')
            exComp{p} = o;
        elseif strcmpi(models{p}.compNames(o),'e')
            exComp{p} = o;
        else
            exComp{p} = 0;
        end
    end
    if exComp{p} == 0
        EM=['The extracellular compartement in model' superModel.OrganismID{p} ' could not be identified'];
        dispEM(EM,false);
    else
        superModel.SubModels{p}.comps(exComp{p})= {'e'};
        superModel.SubModels{p}.compNames(exComp{p})= {'extracellular'};
    end
end

%% obtain extracellular metabolites in each model

exMets = {};
for v=1:length(models)
    exMets{v} = superModel.SubModels{v}.metNames(superModel.SubModels{v}.metComps==exComp{v});
end
