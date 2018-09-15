function dx = f(t,x,superModel,params)

%1.KINETIC BLOCK

Vmax = params.Vmax;
Ks = params.Ks;
mets = params.mets;
metIdx = params.metIdx;


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
gluf = mets(7,2);
watf = mets(8,2);
oxyf = mets(9,2);
phof = mets(10,2);
ammf = mets(11,2);
acef = mets(12,2);
Qf = mets(13,2);
Hf = mets(14,2);
Kf = mets(15,2);
Ff = mets(16,2);
Vf = mets(17,2);
Tf = mets(18,2);
Wf = mets(19,2);
Mf = mets(20,2);
Lf = mets(21,2);
If = mets(22,2);
carf = mets(23,2);
prof = mets(24,2);
butf = mets(25,2);
sucf = mets(26,2);
ethf = mets(27,2);
metf = mets(28,2);
mfaf = mets(29,2);
myrf = mets(30,2);
p28f = mets(31,2);

D = params.other(1);
ACE = params.other(2);
myrMW = params.other(3);
mfaThresh = params.other(4);

bound = [];
%Define bound for each organism and uptake reaction
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(Vmax)
        if Vmax(upMet,organism)==0 & Ks(upMet,organism) ==0
            bound(upMet,organism)=0;
        else
            bound(upMet,organism) = abs((Vmax(upMet,organism) * x(upMet)) / (Ks(upMet,organism) + x(upMet)));
        end
    end
end

lbub = {};
%Define lower bound/upper bound structure containing information about
%which bound needs to update
for organism = 1:length(superModel.organismID) %do this for each organism
    for upMet = 1:length(metIdx) % do this for each uptake metabolite
        if bound(upMet,organism) ==0 | metIdx(upMet,organism) == 0 % check if model is not uptaking this metabolite
            lbub(upMet,organism) = {'0'};
        else
            if abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) > abs(superModel.subModels{organism}.ub(metIdx(upMet,organism)))
                lbub(upMet,organism) = {'l'};
            elseif abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) < abs(superModel.subModels{organism}.ub(metIdx(upMet,organism)))
                lbub(upMet,organism) = {'u'};
            elseif abs(superModel.subModels{organism}.lb(metIdx(upMet,organism))) == abs(superModel.subModels{organism}.ub(metIdx(upMet,organism)))
                dispEM(['Unable to determine whether lower or upper bound should be modified for organism' superModel.organismID{organism} ' ,reaction ' superModel.subModels{organism}.rxnNames{metIdx(upMet,organism)} '. Errors may arise if this occurs for reactions other than water.'],false)
            end
        end
    end
end

%plusmin = [];
%for organism = 1:length(superModel.organismID)
%    for upMet = 1:length(metIdx)
%        if str2num(lbub{upMet,organism}) == 0
%            plusmin(upMet,organism) = 1;
%        else
%            plusmin(upMet,organism) = -1;
%        end    
%    end
%end

%Update bounds based previously defined objects 'bound' and 'lbub'
for organism = 1:length(superModel.organismID)
    for upMet = 1:length(Vmax)
        if bound(upMet,organism)==0
            %no uptake of this metabolite
        else
            if lbub{upMet,organism} == 'l'
                superModel.subModels{organism} = setParam(superModel.subModels{organism},'lb',superModel.subModels{organism}.rxns{metIdx(upMet,organism)},-bound(upMet,organism));
            elseif lbub{upMet,organism} == 'u'
                superModel.subModels{organism} = setParam(superModel.subModels{organism},'ub',superModel.subModels{organism}.rxns{metIdx(upMet,organism)},bound(upMet,organism));
            end
        end
    end
end

%Rule to make S.bo start Myrosinase production and adjust MFalpha2
%production to account for limited protein pool
if x(29)>= mfaThresh
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0.0005); 
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0.0001);
else
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
    %if isempty(solSbo.f)
    %    solSbo.f = 0;
    %    solSbo.x = 0.*ones(length(superModel.subModels{1}.rxns),1);
    %end
    %printFluxes(superModel.subModels{1},solSbo.x)
end
%Need to a checking step to see if protein synthesis load is too heavy
%for S.bo once substrates become limiting
if isempty(solSbo.f) & solSbo.skip == 0
   dispEM('S.bo model not solvable likely due to MFalpha2/Myrosinase production constraints, try reducing protein synthesis load.',false)
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0.00005); 
   superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0.00001);
   solSbo = solveLP(superModel.subModels{1},1);
   FBAsol{1} = solSbo.f;
   %Check if no solution after decreasing protein production by 10%
   if isempty(solSbo.f)
      dispEM('S.bo model still not solvable after lowering protein production by 10 percent!',false)
      superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0.000005); 
      superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0.000001);
      solSbo = solveLP(superModel.subModels{1},1);
      FBAsol{1} = solSbo.f;
      %Check if no solution after decreasing protein production by 100%
      if isempty(solSbo.f)
         dispEM('S.bo model still not solvable after lowering protein production by 100 percent!',false)
         superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4066',0); 
         superModel.subModels{1} = setParam(superModel.subModels{1},'lb','r_4067',0);
         solSbo = solveLP(superModel.subModels{1},1);
         FBAsol{1} = solSbo.f;
         %Check if no solution after stopping protein production
         if isempty(solSbo.f)
            dispEM('S.bo model still not solvable after stopping protein production!',false)
         else
            dispEM('S.bo model solvable after stopping protein production!',false)
         end
      else
         dispEM('S.bo model solvable after lowering protein production by 100 percent!',false)
      end
   else
      dispEM('S.bo model solvable after lowering protein production by 10 percent!',false)
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
    %if isempty(solBth.f)
    %    solBth.f = 0;
    %    solBth.x = 0.*ones(length(superModel.subModels{2}.rxns),1);
    %end
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

%for metNum = 1:length(metIdx) 
%    if metNum == 5
%        dx(metNum)= abs(x(metNum)*solCan.f) - ACE*myrMW*x(30);
%    elseif metNum == 6
%        dx(metNum)= 0;
%    else
%        dx(metNum)=  abs(solSbo.x(metIdx(metNum,1))*x(1))*plusmin(metNum,1) + abs(solBth.x(metIdx(metNum,2))*x(2))*plusmin(metNum,2) + abs(solEre.x(metIdx(metNum,3))*x(3))*plusmin(metNum,3) + abs(solMsi.x(metIdx(metNum,4))*x(4))*plusmin(metNum,4) + abs(solCan.x(metIdx(metNum,5))*x(5))*plusmin(metNum,5) + abs(solCol.x(metIdx(metNum,6))*x(6))*plusmin(metNum,6) - D*(x(metNum)- mets(metNum,2))
%    end
%end


%Biomass mass balances (g/L*hour)
dx(1)= abs(x(1)*solSbo.f) - D*(x(1) - sbof);
dx(2)= abs(x(2)*solBth.f) - D*x(2);
dx(3)= abs(x(3)*solEre.f) - D*x(3);
dx(4)= abs(x(4)*solMsi.f) - D*x(4); 
dx(5)= abs(x(5)*solCan.f) - ACE*myrMW*x(30); % Second term accounts for biomass destroyed by myrosinase action, define as the product of
                                             % ACE: Anti-cancer efficiency (1/hr)
                                             % myrMW: Myrosinase molecular weight(g/mmol) 
                                             % x(30): Myrosinase concentration (mmol/L)
dx(6)= 0;

%Glucose mass balance (mmol/L*hour)
dx(7) = - abs(solSbo.x(1832)*x(1)) - abs(solBth.x(944)*x(2))- abs(solEre.x(108)*x(3)) - abs(solCan.x(3996)*x(5)) - abs(solCol.x(3113)*x(6)) - D*(x(7) - gluf);

%Water mass balance
dx(8) = abs(solSbo.x(1835)*x(1)) - abs(solBth.x(945)*x(2)) - abs(solEre.x(115)*x(3)) + abs(solMsi.x(418)*x(4)) -abs(solCol.x(3127)*x(6)) - D*(x(8)- watf); 

%Oxygen mass balance
dx(9)= - abs(solSbo.x(1903)*x(1)) - abs(solCan.x(4010)*x(5)) - abs(solCol.x(3126)*x(6)) - D*(x(9) - oxyf);

%Phosphate mass balance
dx(10)= -abs(solSbo.x(1909)*x(1)) - abs(solBth.x(948)*x(2)) - abs(solEre.x(116)*x(3)) - abs(solMsi.x(414)*x(4)) - abs(solCan.x(4031)*x(5)) - D*(x(10) - phof );

%Ammonium mass balance
dx(11)= - abs(solSbo.x(1819)*x(1)) - abs(solBth.x(947)*x(2)) - abs(solEre.x(114)*x(3)) - abs(solMsi.x(413)*x(4)) + abs(solCan.x(4032)*x(5)) - D*(x(11) - ammf );

%Acetate mass balance
dx(12)= - abs(solMsi.x(412)*x(4)) - D*(x(12) - acef );

%Glutamine mass balance
dx(13)= - abs(solSbo.x(1878)*x(1)) - abs(solCan.x(4022)*x(5)) - D*(x(13) - Qf );

%Histidine mass balance
dx(14)= - abs(solSbo.x(1879)*x(1)) - abs(solCan.x(4000)*x(5)) - abs(solCol.x(3117)*x(6))- D*(x(14) - Hf );

%Lysine mass balance
dx(15)= - abs(solSbo.x(1883)*x(1)) - abs(solCan.x(4003)*x(5)) - abs(solCol.x(3120)*x(6)) - D*(x(15) - Kf);

%Phenylalanine mass balance
dx(16)= - abs(solSbo.x(1885)*x(1)) - abs(solCan.x(4005)*x(5)) - abs(solCol.x(3122)*x(6)) - D*(x(16) - Ff);

%Valine mass balance
dx(17)= - abs(solSbo.x(1892)*x(1)) - abs(solCan.x(4008)*x(5)) - abs(solCol.x(3125)*x(6)) - D*(x(17) - Vf);

%Threonine mass balance
dx(18)= - abs(solSbo.x(1889)*x(1)) - abs(solCan.x(4006)*x(5)) - abs(solCol.x(3123)*x(6)) - D*(x(18) - Tf);

%Tryptophan mass balance
dx(19)= - abs(solSbo.x(1890)*x(1)) - abs(solBth.x(963)*x(2)) - abs(solEre.x(404)*x(3)) - abs(solCan.x(4007)*x(5)) - abs(solCol.x(3124)*x(6)) - D*(x(19) -  Wf);

%Methionine mass balance
dx(20)= - abs(solSbo.x(1884)*x(1)) - abs(solBth.x(962)*x(2)) - abs(solEre.x(403)*x(3))  - abs(solCan.x(4004)*x(5)) - abs(solCol.x(3121)*x(6)) - D*(x(20) - Mf);

%Leucine mass balance
dx(21)= - abs(solSbo.x(1882)*x(1)) - abs(solCan.x(4002)*x(5)) - abs(solCol.x(3119)*x(6)) - D*(x(21) - Lf);

%Isoleucine mass balance
dx(22)= - abs(solSbo.x(1881)*x(1)) - abs(solCan.x(4001)*x(5)) - abs(solCol.x(3118)*x(6)) - D*(x(22) - If );

%CO2 mass balance
dx(23)=  abs(solSbo.x(1822)*x(1)) + abs(solEre.x(109)*x(3)) - abs(solMsi.x(417)*x(4)) + abs(solCol.x(3128)*x(6)) - D*x(23);

%Propanoate mass balance
dx(24)= abs(solBth.x(954)*x(2)) - D*x(24);

%Butyrate mass balance
dx(25)= abs(solEre.x(111)*x(3)) - D*x(25);

%Succinate mass balance
dx(26)= abs(solSbo.x(1922)*x(1)) - D*x(26);

%Ethanol mass balance
dx(27)= abs(solSbo.x(1841)) - D*x(27);

%Methane mass balance
dx(28)= abs(solMsi.x(412)*x(4))  - D*x(28);%412 index refers to acetate uptake, since 1 mol of methane is produced for 1 mol of acetate

%MFalpha2 mass balance
dx(29)= abs(solSbo.x(1972)*x(1)) - D*(x(29) - mfaf);

%Myrosinase mass balance
dx(30)= abs(solSbo.x(1973)*x(1)) - D*x(30);

%P28 mass balance
dx(31) = abs(solSbo.x(1974)*x(1)) - D*x(31)


dx = dx(:);

end