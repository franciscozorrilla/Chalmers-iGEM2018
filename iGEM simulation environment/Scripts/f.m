function dx = f(t,x,superModel,num)

%Define kinetic parameters
vmax_g1 = num(1); % mmol/gDCW*hr
vmax_g2 = num(2); % mmol/gDCW*hr
vmax_g3 = num(3); % mmol/gDCW*hr
vmax_g5 = num(4); % mmol/gDCW*hr
vmax_g6 = num(5); % mmol/gDCW*hr
ks_g1 = num(6); % mmol/L
ks_g2 = num(7); % mmol/L
ks_g3 = num(8); % mmol/L
ks_g5 = num(9); % mmol/L
ks_g6 = num(10); % mmol/L
vmax_o1 = num(11); % mmol/gDCW*hr
vmax_o5 = num(12); % mmol/gDCW*hr
vmax_o6 = num(13); % mmol/gDCW*hr
ks_o1 = num(14); % mmol/L
ks_o5 = num(15); % mmol/L
ks_o6 = num(16); % mmol/L
vmax_p1 = num(17); % mmol/gDCW*hr
vmax_p2 = num(18); % mmol/gDCW*hr
vmax_p3 = num(19); % mmol/gDCW*hr
vmax_p4 = num(20); % mmol/gDCW*hr
vmax_p5 = num(21); % mmol/gDCW*hr
ks_p1 = num(22); % mmol/L
ks_p2 = num(23); % mmol/L
ks_p3 = num(24); % mmol/L
ks_p4 = num(25); % mmol/L
ks_p5 = num(26); % mmol/L
vmax_ac4 = num(27); % mmol/gDCW*hr
ks_ac4 = num(28); % mmol/L
vmax_am1 = num(29); % mmol/gDCW*hr
vmax_am2 = num(30); % mmol/gDCW*hr
vmax_am3 = num(31); % mmol/gDCW*hr
vmax_am4 = num(32); % mmol/gDCW*hr
ks_am1 = num(33); % mmol/L
ks_am2 = num(34); % mmol/L
ks_am3 = num(35); % mmol/L
ks_am4 = num(36); % mmol/L
vmax_Q5 = num(37); % mmol/gDCW*hr
ks_Q5 = num(38); % mmol/L
vmax_H5 = num(39); % mmol/gDCW*hr
vmax_H6 = num(40); % mmol/gDCW*hr
ks_H5 = num(41); % mmol/L
ks_H6 = num(42); % mmol/L
vmax_K5 = num(43); % mmol/gDCW*hr
vmax_K6 = num(44); % mmol/gDCW*hr
ks_K5 = num(45); % mmol/L
ks_K6 = num(46); % mmol/L
vmax_F5 = num(47); % mmol/gDCW*hr
vmax_F6 = num(48); % mmol/gDCW*hr
ks_F5 = num(49); % mmol/L
ks_F6 = num(50); % mmol/L
vmax_V5 = num(51); % mmol/gDCW*hr
vmax_V6 = num(52); % mmol/gDCW*hr
ks_V5 = num(53); % mmol/L
ks_V6  = num(54); % mmol/L
vmax_T5 = num(55); % mmol/gDCW*hr
vmax_T6 = num(56); % mmol/gDCW*hr
ks_T5 = num(57); % mmol/L
ks_T6 = num(58); % mmol/L
vmax_W5 = num(59); % mmol/gDCW*hr
vmax_W6 = num(60); % mmol/gDCW*hr
ks_W5 = num(61); % mmol/L
ks_W6 = num(62); % mmol/L
vmax_M5 = num(63); % mmol/gDCW*hr
vmax_M6 = num(64); % mmol/gDCW*hr
ks_M5 = num(65); % mmol/L
ks_M6 = num(66); % mmol/L
vmax_L5 = num(67); % mmol/gDCW*hr
vmax_L6 = num(68); % mmol/gDCW*hr
ks_L5 = num(69); % mmol/L
ks_L6 = num(70); % mmol/L
vmax_I5 = num(71); % mmol/gDCW*hr
vmax_I6 = num(72); % mmol/gDCW*hr
ks_I5 = num(73); % mmol/L
ks_I6 = num(74); % mmol/L

%Define flow related terms
sbof = num(104); % g/L
bthf = num(105); % g/L
eref = num(106); % g/L
msif = num(107); % g/L
canf = num(108); % g/L
colf = num(109); % g/L
gluf = num(110); % mmol/L
watf = num(111); % mmol/L
oxyf = num(112); % mmol/L
phof = num(113); % mmol/L
ammf = num(114); % mmol/L
acef = num(115); % mmol/L
Qf = num(116); % mmol/L
Hf = num(117); % mmol/L
Kf = num(118); % mmol/L
Ff = num(119); % mmol/L
Vf = num(120); % mmol/L
Tf = num(121); % mmol/L
Wf = num(122); % mmol/L
Mf = num(123); % mmol/L
Lf = num(124); % mmol/L
If = num(125); % mmol/L
carf = num(126); % mmol/L
prof = num(127); % mmol/L
butf = num(128); % mmol/L
metf = num(129); % mmol/L
mfaf = num(130); % mmol/L
myrf = num(131); % mmol/L
p28f = num(132); % mmol/L
D = num(133); % 1/hr
ACE = num(134);% 1/hr
myrMW = num(135); % g/mmol
mfaThresh = num(136); % mmol/L

%Glucose uptake kinetics
g1 =(vmax_g1*x(7)/(ks_g1+x(7)));
g2 =(vmax_g2*x(7)/(ks_g2+x(7)));
g3 =(vmax_g3*x(7)/(ks_g3+x(7)));
g5 =(vmax_g5*x(7)/(ks_g5+x(7)));
g6 =(vmax_g6*x(7)/(ks_g6+x(7)));

%Water uptake: abundant metabolite, therefore not limiting
h2 =100;
h3 =100;
h6 =100;

%Oxygen uptake kinetics
o1 =(vmax_o1*x(9)/(ks_o1+x(9)));
o5 =(vmax_o5*x(9)/(ks_o5+x(9)));
o6 =(vmax_o6*x(9)/(ks_o6+x(9)));

%Phosphate uptake kinetics
p1 =(vmax_p1*x(10)/(ks_p1+x(10)));
p2 =(vmax_p2*x(10)/(ks_p2+x(10)));
p3 =(vmax_p3*x(10)/(ks_p3+x(10)));
p4 =(vmax_p4*x(10)/(ks_p4+x(10)));
p5 =(vmax_p5*x(10)/(ks_p5+x(10)));

%Ammonium uptake kinetics
am1 =(vmax_am1*x(11)/(ks_am1+x(11)));
am2 =(vmax_am2*x(11)/(ks_am2+x(11)));
am3 =(vmax_am3*x(11)/(ks_am3+x(11)));
am4 =(vmax_am4*x(11)/(ks_am4+x(11)));

%Acetate uotake kinetics
ac4 =(vmax_ac4*x(12)/(ks_ac4+x(12)));

%Glutamine uptake kinetics
Q5 = (vmax_Q5*x(13)/(ks_Q5+x(13)));

%Histidine uptake kinetics
H5 = (vmax_H5*x(14)/(ks_H5+x(14)));
H6 = (vmax_H6*x(14)/(ks_H6+x(14)));

%Lysine uptake kinetics
K5 = (vmax_K5*x(15)/(ks_K5+x(15)));
K6 = (vmax_K6*x(15)/(ks_K6+x(15)));

%Phenylalanine uptake kinetics
F5 = (vmax_F5*x(16)/(ks_F5+x(16)));
F6 = (vmax_F6*x(16)/(ks_F6+x(16)));

%Valine uptake kinetics
V5 = (vmax_V5*x(17)/(ks_V5+x(17)));
V6 = (vmax_V6*x(17)/(ks_V6+x(17)));

%Threonine uptake kinetics
T5 = (vmax_T5*x(18)/(ks_T5+x(18)));
T6 = (vmax_T6*x(18)/(ks_T6+x(18)));

%Tryptophan uptake kinetics
W5 = (vmax_W5*x(19)/(ks_W5+x(19)));
W6 = (vmax_W6*x(19)/(ks_W6+x(19)));

%Methionine uptake kinetics
M5 = (vmax_M5*x(20)/(ks_M5+x(20)));
M6 = (vmax_M6*x(20)/(ks_M6+x(20)));

%Leucine uptake kinetics
L5 = (vmax_L5*x(21)/(ks_L5+x(21)));
L6 = (vmax_L6*x(21)/(ks_L6+x(21)));

%Isoleucine uptake kinetics
I5 = (vmax_I5*x(22)/(ks_I5+x(22)));
I6 = (vmax_I6*x(22)/(ks_I6+x(22)));

%Carbon dioxide uptake kinetics
c4 = (5*x(23)/(5+x(23)));

%Rule to make S.bo start Myrosinase production and adjust MFalpha2
%production to account for limited protein pool
if x(27)>= mfaThresh
   superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_4066',0.0001); 
   superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_4067',0.0001);
end

%Update glucose consumption
superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_1714',-g1); %g1 takes a negative sign to account for the way this uptake reaction is defined
superModel.SubModels{2} = setParam(superModel.SubModels{2},'ub','glcIn',g2);
superModel.SubModels{3} = setParam(superModel.SubModels{3},'ub','glcIn',g3);
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9034',g5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9034',g6);

%Update water consumption
superModel.SubModels{2} = setParam(superModel.SubModels{2},'ub','EX_h2o_LPAREN_e_RPAREN_Up',h2);
superModel.SubModels{3} = setParam(superModel.SubModels{3},'lb','r113',-h3);% h3 should be negative to account for uptake reaction definition
superModel.SubModels{6} = setParam(superModel.SubModels{6},'lb','HMR_9034',-h6);%h6 should be negative to account for uptake reaction definition

%Update Oxygen consumption
superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_1992',-o1);
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9048',o5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9048',o6);

%Update phosphate consumption
superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_2005',-p1);
superModel.SubModels{2} = setParam(superModel.SubModels{2},'ub','EX_pi_LPAREN_e_RPAREN_',p2);
superModel.SubModels{3} = setParam(superModel.SubModels{3},'ub','r114',p3);
superModel.SubModels{4} = setParam(superModel.SubModels{4},'ub','piIN',p4);
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9072',p5);

%Update ammonium consumption
superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_1654',-am1);
superModel.SubModels{2} = setParam(superModel.SubModels{2},'ub','EX_nh4_LPAREN_e_RPAREN_',am2);
superModel.SubModels{3} = setParam(superModel.SubModels{3},'ub','r112',am3);
superModel.SubModels{4} = setParam(superModel.SubModels{4},'ub','NH4IN',am4);

%Update acetate consumption
superModel.SubModels{4} = setParam(superModel.SubModels{4},'ub','actIn',ac4);

%Update Glutamine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9063',Q5);

%Update Histidine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9038',H5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9038',H6);

%Update Lysine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9041',K5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9041',K6);

%Update Phenylalanine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9043',F5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9043',F6);

%Update Valine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9046',V5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9046',V6);

%Update Threonine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9044',T5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9044',T6);

%Update Tryptophan consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9045',W5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9045',W6);

%Update Methionine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9042',M5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9042',M6);

%Update Leucine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9040',L5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9040',L6);

%Update Isoleucine consumption
superModel.SubModels{5} = setParam(superModel.SubModels{5},'ub','HMR_9039',I5);
superModel.SubModels{6} = setParam(superModel.SubModels{6},'ub','HMR_9039',I6);

%Update carbon dioxide consumption
superModel.SubModels{4} = setParam(superModel.SubModels{4},'ub','CO2In',c4);

%FBA for each model
solSbo = solveLP(superModel.SubModels{1},1);
FBAsol{1} = solSbo.f;
%Need to include this check to see if protein synthesis load is too heavy
%for S.bo once substrates become limiting
if isempty(FBAsol{1})
   dispEM('S.bo model was not able to be solved likely due to MFalpha2/Myrosinase production constraints, production stopped to check if model is solvable',false)
   superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_4066',0); 
   superModel.SubModels{1} = setParam(superModel.SubModels{1},'lb','r_4067',0);
   FBAsol{1} = solSbo.f;
end

solBth = solveLP(superModel.SubModels{2},1);
solEre = solveLP(superModel.SubModels{3},1);
solMsi = solveLP(superModel.SubModels{4},1);
solCan = solveLP(superModel.SubModels{5},1);
solCol = solveLP(superModel.SubModels{6},1);

FBAsol{2} = solBth.f;
FBAsol{3} = solEre.f;
FBAsol{4} = solMsi.f;
FBAsol{5} = solCan.f;
FBAsol{6} = solCol.f;

FBAsol %to show progress

%Need to set solution to 0 if unfeasible, sets as blank by default and
%causes errors otherwise

if isempty(solSbo.f)
    solSbo.f = 0;
    solSbo.x = 0.*ones(length(superModel.SubModels{1}.rxns),1);
    dispEM('S.bo model was not able to be solved',false)
end
if isempty(solBth.f)
    solBth.f = 0;
    solBth.x = 0.*ones(length(superModel.SubModels{2}.rxns),1);
    dispEM('B.th model was not able to be solved',false)
end
if isempty(solEre.f)
   solEre.f = 0; 
   solEre.x = 0.*ones(length(superModel.SubModels{3}.rxns),1);
   dispEM('E.re model was not able to be solved',false)
end
if isempty(solMsi.f)
   solMsi.f = 0; 
   solMsi.x = 0.*ones(length(superModel.SubModels{4}.rxns),1);
   dispEM('M.si model was not able to be solved',false)
end
if isempty(solCan.f)
   solCan.f = 0; 
   solCan.x = 0.*ones(length(superModel.SubModels{5}.rxns),1);
   dispEM('Cancer model was not able to be solved',false)
end
if isempty(solCol.f)
    solCol.f = 0;
    solCol.x = 0.*ones(length(superModel.SubModels{6}.rxns),1) ;
    dispEM('Colon model was not able to be solved',false)
end

%Biomass mass balances (g/L*hour)
dx(1)= abs(x(1)*solSbo.f) - D*(x(1) - sbof);
dx(2)= abs(x(2)*solBth.f) - D*x(2);
dx(3)= abs(x(3)*solEre.f) - D*x(3);
dx(4)= abs(x(4)*solMsi.f) - D*x(4); 
dx(5)= abs(x(5)*solCan.f) - ACE*myrMW*x(28); % Second term accounts for biomass destroyed by myrosinase action, define as the product of
                                             % ACE: Anti-cancer efficiency (1/hr)
                                             % myrMW: Myrosinase molecular weight(g/mmol) 
                                             % x(28): Myrosinase concentration (mmol/L)
dx(6)= 0;

%Glucose mass balance (mmol/L*hour)
dx(7) = - abs(solSbo.x(1808)*x(1)) - abs(solBth.x(944)*x(2))- abs(solEre.x(108)*x(3)) - abs(solCan.x(3996)*x(5)) - abs(solCol.x(3113)*x(6)) - D*(x(7) - gluf);

%Water mass balance
dx(8) = -abs(solBth.x(945)*x(2)) - abs(solEre.x(115)*x(3)) + abs(solMsi.x(418)*x(4)) -abs(solCol.x(3127)*x(6)) - D*(x(8)- watf); 

%Oxygen mass balance
dx(9)= - abs(solSbo.x(1811)*x(1)) - abs(solCan.x(4010)*x(5)) - abs(solCol.x(3126)*x(6)) - D*(x(9) - oxyf);

%Phosphate mass balance
dx(10)= -abs(solSbo.x(1812)*x(1)) - abs(solBth.x(948)*x(2)) - abs(solEre.x(116)*x(3)) - abs(solMsi.x(414)*x(4)) - abs(solCan.x(4031)*x(5)) - D*(x(10) - phof );

%Ammonium mass balance
dx(11)= - abs(solSbo.x(1806)*x(1)) - abs(solBth.x(947)*x(2)) - abs(solEre.x(114)*x(3)) - abs(solMsi.x(413)*x(4)) + abs(solCan.x(4032)*x(5)) - D*(x(11) - ammf );

%Acetate mass balance
dx(12)= - abs(solMsi.x(412)*x(4)) - D*(x(12) - acef );

%Glutamine mass balance
dx(13)= - abs(solCan.x(4022)*x(5)) - D*(x(13) - Qf );

%Histidine mass balance
dx(14)= - abs(solCan.x(4000)*x(5)) - abs(solCol.x(3117)*x(6))- D*(x(14) - Hf );

%Lysine mass balance
dx(15)= - abs(solCan.x(4003)*x(5)) - abs(solCol.x(3120)*x(6)) - D*(x(15) - Kf);

%Phenylalanine mass balance
dx(16)= - abs(solCan.x(4005)*x(5)) - abs(solCol.x(3122)*x(6)) - D*(x(16) - Ff);

%Valine mass balance
dx(17)= - abs(solCan.x(4008)*x(5)) - abs(solCol.x(3125)*x(6)) - D*(x(17) - Vf);

%Threonine mass balance
dx(18)= - abs(solCan.x(4006)*x(5)) - abs(solCol.x(3123)*x(6)) - D*(x(18) - Tf);

%Tryptophan mass balance
dx(19)= - abs(solCan.x(4007)*x(5)) - abs(solCol.x(3124)*x(6)) - D*(x(19) -  Wf);

%Methionine mass balance
dx(20)= - abs(solCan.x(4004)*x(5)) - abs(solCol.x(3121)*x(6)) - D*(x(20) - Mf);

%Leucine mass balance
dx(21)= - abs(solCan.x(4002)*x(5)) - abs(solCol.x(3119)*x(6)) - D*(x(21) - Lf);

%Isoleucine mass balance
dx(22)= - abs(solCan.x(4001)*x(5)) - abs(solCol.x(3118)*x(6)) - D*(x(22) - If );

%CO2 mass balance
dx(23)=  abs(solSbo.x(1807)*x(1)) + abs(solEre.x(109)*x(3)) - abs(solMsi.x(417)*x(4)) + abs(solCol.x(3128)*x(6)) - D*x(23);

%Propanoate mass balance
dx(24)= abs(solBth.x(954)*x(2)) - D*x(24);

%Butyrate mass balance
dx(25)= abs(solEre.x(111)*x(3)) - D*x(25);

%Methane mass balance
dx(26)= abs(solMsi.x(412)*x(4))  - D*x(26);%412 index refers to acetate uptake, since 1 mol of methane is produced for 1 mol of acetate

%MFalpha2 mass balance
dx(27)= abs(solSbo.x(1818)*x(1)) - D*(x(27) - mfaf);

%Myrosinase mass balance
dx(28)= abs(solSbo.x(1819)*x(1)) - D*x(28);

%P28 mass balance
dx(29) = abs(solSbo.x(1820)*x(1)) - D*x(29)

dx = dx(:);

end