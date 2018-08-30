%% adding protein pool doesnt seem to help 

modelMsi.genes = {''}; %insert dummy field to avoid errors with addRxns()
clear metsToAdd % Make sure this structure is empty
metsToAdd = {}; %Add proteinPool pseudo metabolite
metsToAdd.mets = {'proteinPool'};
metsToAdd.metNames = {'proteinPool'};
metsToAdd.compartments = {'s'};
modelMsi=addMets(modelMsi,metsToAdd);

clear rxnToAdd % Make sure this structure is empty
rxnToAdd.rxns={'EXC_IN_m90000c','human_proteinPool'};% add AA exchange rxn and artificial AA pool rxn 
rxnToAdd.rxnNames={'EXC_IN_m90000c','human_proteinPool'};
rxnToAdd.equations={' => proteinPool[s]','0.0937 L-Alanine[s] + 0.0507 L-Arginine[s] + 0.04 L-Asparagine[s] + 0.04 L-Aspartate[s] + 0.0142 L-Cysteine[s] + 0.1118 L-Glutamine[s] + 0.1831 Glycine[s] + 0.0198 L-Histidine[s] + 0.0309 L-Isoleucine[s] + 0.0664 L-Leucine[s] + 0.0571 L-Lysine[s] + 0.0156 L-Methionine[s] + 0.029 L-Phenylalanine[s] + 0.0853 L-Proline[s] + 0.0491 L-Serine[s] + 0.0402 L-Threonine[s] + 0.0072 L-Tryptophan[s] + 0.019 L-Tyrosine[s] + 0.0471 L-Valine[s] <=> proteinPool[s]'};
rxnToAdd.lb=[-1000,-1000];
rxnToAdd.ub=[1000,1000];
modelMsi=addRxns(modelMsi,rxnToAdd,3,'',true);
%% doesnt help
unconstrain = [];
for q =1:length(modelMsi.rxns)
    if modelMsi.ub(q)==0 && modelMsi.lb(q)==0 && modelMsi.rev(q)==1
        unconstrain = [unconstrain;q];
    end
     modelMsi.ub(q) = 1000;
end
modelMsi = setParam(modelMsi,'ub',modelMsi.rxns(unconstrain),1000);
modelMsi = setParam(modelMsi,'lb',modelMsi.rxns(unconstrain),-1000);

%% doesnt work
constrain = [];
for q =1:length(modelMsi.rxns)
    if modelMsi.rev(q)==0
        constrain = [constrain;q];
    end
end
modelMsiConstrained = setParam(modelMsi,'lb',modelMsi.rxns(constrain),0);
solMsiConstrained=solveLP(modelMsiConstrained,1);
printFluxes(modelMsiConstrained,solMsiConstrained.x);