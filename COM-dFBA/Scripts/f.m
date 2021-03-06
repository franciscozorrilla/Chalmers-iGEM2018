function dx = f(t,x,superModel,params)

%1.KINETIC BLOCK

Vmax = params.Vmax;
Ks = params.Ks;
mets = params.mets;
metIdx = params.metIdx;
blocked = params.blocked;

%Define parameters based on params structure in order to make mass balances
%and other expressions legible/understandable 
sbo0 = mets(1,1);
bth0 = mets(2,1);
ere0 = mets(3,1);
msi0 = mets(4,1);
can0 = mets(5,1);
col0 = mets(6,1);

sbof = mets(1,2);
bthf = mets(2,2);
eref = mets(3,2);
msif = mets(4,2);
canf = mets(5,2);
colf = mets(6,2);

D = params.other(1);
ACE = params.other(2);
myrMW = params.other(3);
mfaThresh = params.other(4);

bound = [];
%Define bound for each organism and uptake reaction
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(Vmax)
        if (Vmax(upMet,organism)==0 & Ks(upMet,organism) ==0) | x(upMet)==0
            bound(upMet,organism)=0;
        else
            bound(upMet,organism) = abs((Vmax(upMet,organism) * x(upMet)) / (Ks(upMet,organism) + x(upMet)));
        end
    end
end

equationStrings = {};
%Construct equation strings for each exchange metabolite
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(metIdx)
        if metIdx(upMet,organism) ~= blocked(1,organism)
            equationStrings{upMet,organism}=constructEquations(superModel.subModels{organism},metIdx(upMet,organism));
        end
    end
end

lbub = {};
%Define lower bound/upper bound structure containing information about
%which bound needs to update: only for uptake metabolites
for organism = 1:length(superModel.organismID) %do this for each organism
    for upMet = 1:length(metIdx) % do this for each uptake metabolite
        if metIdx(upMet,organism) == blocked(1,organism) % check if model is not uptaking this metabolite
            lbub(upMet,organism) = {'x'};
        else
            if (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) > abs(superModel.subModels{organism}.ub(metIdx(upMet,organism))))&(strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))>0.5)
                lbub(upMet,organism) = {'l'};
            elseif (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) < abs(superModel.subModels{organism}.ub(metIdx(upMet,organism))))&(strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))<0.5)
                lbub(upMet,organism) = {'u'};
            elseif (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) == abs(superModel.subModels{organism}.ub(metIdx(upMet,organism)))) & (abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) ~= 0)
                dispEM(['Unable to determine whether lower or upper bound should be modified for organism' superModel.organismID{organism} ' ,reaction involving metabolite number ' upMet '. Errors may arise if this occurs for reactions other than water.'],false)
            else
                lbub(upMet,organism) = {'x'};
            end
        end
    end
end

%Update bounds based previously defined objects 'bound','lbub', and 'equationStrings'
for organism = 1:length(superModel.organismID)
    for upMet = (length(superModel.organismID)+1):length(Vmax)
        if strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))>0.5 & lbub{upMet,organism}=='u'  %met on rxntnt side & lb<ub: secreted metabolite
            %do nothing
        elseif strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))<0.5 & lbub{upMet,organism}=='u'  %met on product side & lb<ub: uptake metabolite
            superModel.subModels{organism} = setParam(superModel.subModels{organism},'ub',superModel.subModels{organism}.rxns{metIdx(upMet,organism)},bound(upMet,organism));%update bound
        elseif strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))>0.5 & lbub{upMet,organism}=='l'  %met on rxntnt side & lb>ub: uptake metabolite
            superModel.subModels{organism} = setParam(superModel.subModels{organism},'lb',superModel.subModels{organism}.rxns{metIdx(upMet,organism)},-bound(upMet,organism));%update bound
        elseif strfind(cell2mat(equationStrings{upMet,organism}),'=')/length(cell2mat(equationStrings{upMet,organism}))<0.5 & lbub{upMet,organism}=='l'  %met on product side & lb<ub: secreted metabolite
            %do nothing
        end
    end
end

%Rule to make S.bo start Myrosinase production: NOTE THAT INDEX IS
%HARDCODED, NEED TO ADJUST NUMBER IF METABOLITES ARE ADDED OR REMOVED
if x(28)>= mfaThresh
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0.001); 
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0.0001);
elseif x(28)< mfaThresh
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0.001); 
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0);
end

%2.FBA

%Perform FBA for modelSbo if organism is present in simulation
if sbo0 == 0 & sbof == 0
    solSbo.f = 0;
    solSbo.x = 0.*ones(length(superModel.subModels{1}.rxns),1);
    solSbo.skip = 1;
else
    solSbo = solveLP(superModel.subModels{1},1);
    FBAsol{1} = solSbo.f;
    solSbo.skip = 0;
end
%Need to a checking step to see if protein synthesis load is too heavy
%for S.bo once substrates become limiting
if isempty(solSbo.f) & solSbo.skip == 0
   dispEM('S.bo model not solvable likely due to MFalpha2/Myrosinase production constraints, try reducing protein synthesis load.',false)
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0.0001); 
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0.00001);
   solSbo = solveLP(superModel.subModels{1},1);
   FBAsol{1} = solSbo.f;
   %Check if no solution after decreasing protein production by one order
   %of magnitude
   if isempty(solSbo.f)
      dispEM('S.bo model still not solvable after lowering protein production by one order of magnitude!',false)
      superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0.00001); 
      superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0.000001);
      solSbo = solveLP(superModel.subModels{1},1);
      FBAsol{1} = solSbo.f;
      %Check if no solution after decreasing protein production by two
      %orders of magnitude
      if isempty(solSbo.f)
         dispEM('S.bo model still not solvable after lowering protein production by two orders of magnitude!',false)
         superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0); 
         superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0);
         solSbo = solveLP(superModel.subModels{1},1);
         FBAsol{1} = solSbo.f;
         %Check if no solution after stopping protein production
         if isempty(solSbo.f)
            dispEM('S.bo model still not solvable after stopping protein production.',false)
         else
            dispEM('S.bo model solvable after stopping protein production. Further simulation may be more computationally expensive if this message persists.',false)
         end
      else
         dispEM('S.bo model solvable after lowering protein production by 100 percent.',false)
      end
   else
      dispEM('S.bo model solvable after lowering protein production by 10 percent.',false)
   end
end

%Perform FBA for the rest of the models if they are present in simulation
if bth0 == 0 & bthf == 0
    solBth.f = 0;
    solBth.x = 0.*ones(length(superModel.subModels{2}.rxns),1);
    FBAsol{2} = {};
else
    solBth = solveLP(superModel.subModels{2},1);
    FBAsol{2} = solBth.f;
    if isempty(solBth.f)
        solBth.f = 0;
        solBth.x = 0.*ones(length(superModel.subModels{2}.rxns),1);
    end
    %printFluxes(superModel.subModels{2},solBth.x)
end

if ere0 == 0 & eref == 0
    solEre.f = 0; 
    solEre.x = 0.*ones(length(superModel.subModels{3}.rxns),1);
    FBAsol{3} = {};
else
    solEre = solveLP(superModel.subModels{3},1);
    FBAsol{3} = solEre.f;
end

if msi0 == 0 & msif == 0
    solMsi.f = 0; 
    solMsi.x = 0.*ones(length(superModel.subModels{4}.rxns),1);
    FBAsol{4} = {};
else
    solMsi = solveLP(superModel.subModels{4},1);
    FBAsol{4} = solMsi.f;
end

if can0 == 0 & canf == 0
    solCan.f = 0; 
    solCan.x = 0.*ones(length(superModel.subModels{5}.rxns),1);
    FBAsol{5} = {};
else
    solCan = solveLP(superModel.subModels{5},1);
    FBAsol{5} = solCan.f;
end

if col0 == 0 & colf == 0
    solCol.f = 0;
    solCol.x = 0.*ones(length(superModel.subModels{6}.rxns),1) ;
    FBAsol{6} = {};
else
    solCol = solveLP(superModel.subModels{6},1);
    FBAsol{6} = solCol.f;
end

%show progress
FBAsol

%Need to set solution to 0 if no solution, otherwise sets as blank by  
%default and causes errors. No need to run loop if model is not present in
%simulation (eg: solSbo.skip =1 as determined above)

if isempty(solSbo.f)
    solSbo.f = 0;
    solSbo.x = 0.*ones(length(superModel.subModels{1}.rxns),1);
    dispEM('S.bo model was not able to be solved',false)
end

if isempty(solBth.f)
    solBth.f = 0;
    solBth.x = 0.*ones(length(superModel.subModels{2}.rxns),1);
    dispEM('B.th model was not able to be solved',false)
end

if isempty(solEre.f)
    solEre.f = 0; 
    solEre.x = 0.*ones(length(superModel.subModels{3}.rxns),1);
    dispEM('E.re model was not able to be solved',false)
end

if isempty(solMsi.f)
    solMsi.f = 0; 
    solMsi.x = 0.*ones(length(superModel.subModels{4}.rxns),1);
    dispEM('M.si model was not able to be solved',false)
end

if isempty(solCan.f)
    solCan.f = 0; 
    solCan.x = 0.*ones(length(superModel.subModels{5}.rxns),1);
    dispEM('Cancer model was not able to be solved',false)
end

if isempty(solCol.f)
    solCol.f = 0;
    solCol.x = 0.*ones(length(superModel.subModels{6}.rxns),1) ;
    dispEM('Colon model was not able to be solved',false)
end

%3.DYNAMIC BLOCK

%Need to determine if flux should be adding or subtracting to each Mass
%Balance
plusmin = [];
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(metIdx)
        if lbub{upMet,organism} == 'x'
            plusmin(upMet,organism) = 1;
        elseif (lbub{upMet,organism} == 'l'|lbub{upMet,organism} == 'u')
            plusmin(upMet,organism) = -1;
        end    
    end
end

%Run loop to generate mass balances
%NOTE: Indexes are hardcoded, so this loop should be modified if more 
%models are to be included in simulation or order changes.
for metNum = 1:length(metIdx)
    if metNum == 1 & sbo0 == 0 & sbof == 0
        dx(metNum) = 0;
    elseif metNum == 2 & bth0 == 0 & bthf == 0
        dx(metNum) = 0;    
    elseif metNum == 3 & ere0 == 0 & eref == 0
        dx(metNum) = 0;
    elseif metNum == 4 & msi0 == 0 & msif == 0
        dx(metNum) = 0;
    elseif metNum == 5 & can0 == 0 & canf == 0
        dx(metNum) = 0;
    elseif metNum == 5 & (can0 ~= 0 | canf ~= 0)
        dx(metNum)= abs(x(metNum)*solCan.f) - ACE*myrMW*x(29);%Cancer biomass mass balance: NOTE THAT INDEX IS HARDCODED, NEED TO ADJUST NUMBER IF METABOLITES ARE ADDED OR REMOVED
    elseif metNum == 6
        dx(metNum)= 0;%Colon biomass mass balance should not change
    else
        dx(metNum)=  abs(solSbo.x(metIdx(metNum,1))*x(1))*plusmin(metNum,1) + ... 
                     abs(solBth.x(metIdx(metNum,2))*x(2))*plusmin(metNum,2) + ... 
                     abs(solEre.x(metIdx(metNum,3))*x(3))*plusmin(metNum,3) + ... 
                     abs(solMsi.x(metIdx(metNum,4))*x(4))*plusmin(metNum,4) + ...
                     abs(solCan.x(metIdx(metNum,5))*x(5))*plusmin(metNum,5) + ...
                     abs(solCol.x(metIdx(metNum,6))*x(6))*plusmin(metNum,6) - ...
                     D*(x(metNum)- mets(metNum,2));
    end
end

%This is just because the last metabolite,p28, is not used in these
%simulations. For future usage this last line can be removed.
dx(length(metIdx) ) = 0 

%Do this to avoid errors with function
dx = dx(:);

%Show current time
t

%It is very helpful to print out exchnage fluxes for an organism when
%debuggning/troubleshooting. Uncomment one or multiple of the following
%lines to see exchange fluxes during simulations:

printFluxes(superModel.subModels{1},solSbo.x)
%printFluxes(superModel.subModels{2},solBth.x)
%printFluxes(superModel.subModels{3},solEre.x)
%printFluxes(superModel.subModels{4},solMsi.x)
%printFluxes(superModel.subModels{5},solCan.x)
%printFluxes(superModel.subModels{6},solCol.x)

end