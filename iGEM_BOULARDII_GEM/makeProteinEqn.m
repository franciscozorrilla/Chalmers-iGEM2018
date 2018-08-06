function eqn = makeProteinEqn(proteinName,inputSeq)
% makeProteinEqn
%   Parses protein forming equation for proteinName using inputSeq
%   
%   proteinName     String containing protein name
%   inputSeq        String containing amino acid sequence of protein
%
%   Example: If multiple proteins are being added
%
%   protStruct = {}; %create empty protein structure
%   protStruct.name = {'MFalpha2','Myrosinase','P28'}; %add protein names
%   protStruct.seq = {'MKFISTFLTFILAAVSVTASSDEDIAQVPAEAIIGYLDFGGDHDIAFLPFSNATASGLLFINTTIAEAAEKEQNTTLAKREAVADAWHWLNLRPGQPMYKREANADAWHWLQLKPGQPMY', ...
%    'MKHLGLILAFLLALATCKADEEITCEENLPFKCSQPDRLNSSSFEKDFIFGVASSAYQACCLGRGLNVWDGFTHRYPNKSGPDHGNGDTTCDSFSYWQKDIDVLDELNATGYRFSIAWSRIIPRGKRSRGVNKDGINYYHGLIDGLIDKGITPFVTLFHWDLPQVLQDEYEGFLDPQIIHDFKHYANLCFQEFGHKVKNWLTINQLYTVPTRGYGAGSDAPGRCSPMVDPTCYAGNSSTEPYIVAHNQLLAHATVVDLYRKNYSIGPVMITRWFLPYNDTDPDSIAATERMKEFFLGWFMGPLTNGTYPQIMIDTVGERLPSFSPEESNLVKGSYDYLGLNYYVTQYAQPSPNPVHWANHTAMMDAGAKLTFRGNSDETKNSYYYPKGIYYVMDYFKTKYYNPLIYVTENGISTPGNETRDESMLHYKRIEYLCSHLCFLSKVIKEKHVNVKGYFAWSLGDNYEFDKGFTVRFGLSYIDWNNVTDRDLKLSGKWYQKFISPAIKNPLKKDFLRSSLTFEKNKKFEDA', ...
%    'LSTAADMQGVVTDGMASGLDKDYLKPDD'}; %add amino acid sequences
%
%   for x = 1:length(protStruct.name) %create an equation for each protein using makeProteinEqn()
%      protStruct.eqn{x} = makeProteinEqn(protStruct.name{x},protStruct.seq{x});
%   end
%
%   protStruct.eqn{1:3}
%
%   '19 Ala-tRNA(Ala)[c] + 3 Arg-tRNA(Arg)[c] + 5 Asn-tRNA(Asn)[c] + 7 Asp-tRNA(Asp)[c] + 0 Cys-tRNA(Cys)[c] + 5 Gln-tRNA(Gln)[c] + 7 Glu-tRNA(Glu)[c] + 6 Gly-tRNA(Gly)[c] + 3 His-tRNA(His)[c] + 8 Ile-tRNA(Ile)[c] + 11 Leu-tRNA(Leu)[c] + 5 Lys-tRNA(Lys)[c] + 3 Met-tRNA(Met)[c] + 7 Phe-tRNA(Phe)[c] + 6 Pro-tRNA(Pro)[c] + 6 Ser-tRNA(Ser)[c] + 8 Thr-tRNA(Thr)[c] + 4 Trp-tRNA(Trp)[c] + 3 Tyr-tRNA(Tyr)[c] + 4 Val-tRNA(Val)[c] + 479 ATP[c] + 479 H2O[c] => 19 tRNA(Ala)[c] + 3 tRNA(Arg)[c] + 5 tRNA(Asn)[c] + 7 tRNA(Asp)[c] + 0 tRNA(Cys)[c] + 5 tRNA(Gln)[c] + 7 tRNA(Glu)[c] + 6 tRNA(Gly)[c] + 3 tRNA(His)[c] + 8 tRNA(Ile)[c] + 11 tRNA(Leu)[c] + 5 tRNA(Lys)[c] + 3 tRNA(Met)[c] + 7 tRNA(Phe)[c] + 6 tRNA(Pro)[c] + 6 tRNA(Ser)[c] + 8 tRNA(Thr)[c] + 4 tRNA(Trp)[c] + 3 tRNA(Tyr)[c] + 4 tRNA(Val)[c] + 479 ADP[c] + 479 phosphate[c] + 479 H+[c] + MFalpha2[c]'
%   '26 Ala-tRNA(Ala)[c] + 20 Arg-tRNA(Arg)[c] + 31 Asn-tRNA(Asn)[c] + 36 Asp-tRNA(Asp)[c] + 11 Cys-tRNA(Cys)[c] + 13 Gln-tRNA(Gln)[c] + 24 Glu-tRNA(Glu)[c] + 39 Gly-tRNA(Gly)[c] + 15 His-tRNA(His)[c] + 28 Ile-tRNA(Ile)[c] + 44 Leu-tRNA(Leu)[c] + 35 Lys-tRNA(Lys)[c] + 10 Met-tRNA(Met)[c] + 30 Phe-tRNA(Phe)[c] + 28 Pro-tRNA(Pro)[c] + 35 Ser-tRNA(Ser)[c] + 31 Thr-tRNA(Thr)[c] + 11 Trp-tRNA(Trp)[c] + 36 Tyr-tRNA(Tyr)[c] + 24 Val-tRNA(Val)[c] + 2107 ATP[c] + 2107 H2O[c] => 26 tRNA(Ala)[c] + 20 tRNA(Arg)[c] + 31 tRNA(Asn)[c] + 36 tRNA(Asp)[c] + 11 tRNA(Cys)[c] + 13 tRNA(Gln)[c] + 24 tRNA(Glu)[c] + 39 tRNA(Gly)[c] + 15 tRNA(His)[c] + 28 tRNA(Ile)[c] + 44 tRNA(Leu)[c] + 35 tRNA(Lys)[c] + 10 tRNA(Met)[c] + 30 tRNA(Phe)[c] + 28 tRNA(Pro)[c] + 35 tRNA(Ser)[c] + 31 tRNA(Thr)[c] + 11 tRNA(Trp)[c] + 36 tRNA(Tyr)[c] + 24 tRNA(Val)[c] + 2107 ADP[c] + 2107 phosphate[c] + 2107 H+[c] + Myrosinase[c]'
%   '3 Ala-tRNA(Ala)[c] + 0 Arg-tRNA(Arg)[c] + 0 Asn-tRNA(Asn)[c] + 6 Asp-tRNA(Asp)[c] + 0 Cys-tRNA(Cys)[c] + 1 Gln-tRNA(Gln)[c] + 0 Glu-tRNA(Glu)[c] + 3 Gly-tRNA(Gly)[c] + 0 His-tRNA(His)[c] + 0 Ile-tRNA(Ile)[c] + 3 Leu-tRNA(Leu)[c] + 2 Lys-tRNA(Lys)[c] + 2 Met-tRNA(Met)[c] + 0 Phe-tRNA(Phe)[c] + 1 Pro-tRNA(Pro)[c] + 2 Ser-tRNA(Ser)[c] + 2 Thr-tRNA(Thr)[c] + 0 Trp-tRNA(Trp)[c] + 1 Tyr-tRNA(Tyr)[c] + 2 Val-tRNA(Val)[c] + 111 ATP[c] + 111 H2O[c] => 3 tRNA(Ala)[c] + 0 tRNA(Arg)[c] + 0 tRNA(Asn)[c] + 6 tRNA(Asp)[c] + 0 tRNA(Cys)[c] + 1 tRNA(Gln)[c] + 0 tRNA(Glu)[c] + 3 tRNA(Gly)[c] + 0 tRNA(His)[c] + 0 tRNA(Ile)[c] + 3 tRNA(Leu)[c] + 2 tRNA(Lys)[c] + 2 tRNA(Met)[c] + 0 tRNA(Phe)[c] + 1 tRNA(Pro)[c] + 2 tRNA(Ser)[c] + 2 tRNA(Thr)[c] + 0 tRNA(Trp)[c] + 1 tRNA(Tyr)[c] + 2 tRNA(Val)[c] + 111 ADP[c] + 111 phosphate[c] + 111 H+[c] + P28[c]'
%
%   Usage: eqn = makeProteinEqn(proteinName,inputSeq)
%
%   NOTE: This function was written for the metabolite nomenclature of the
%   Chalmers Sysbio Saccharomyces cereviae GEM available here:
%   https://github.com/SysBioChalmers/yeast-GEM. The amino acid structure
%   (AAstruct) should be modified accordingly if a different amino acid 
%   nomenclature is being used.
%
%   Francisco Zorrilla, 2018-08-02 

AAstruct = {'Ala-tRNA(Ala)[c]' 'A'
'Arg-tRNA(Arg)[c]' 'R'
'Asn-tRNA(Asn)[c]' 'N'
'Asp-tRNA(Asp)[c]' 'D'
'Cys-tRNA(Cys)[c]' 'C'
'Gln-tRNA(Gln)[c]' 'Q'
'Glu-tRNA(Glu)[c]' 'E'
'Gly-tRNA(Gly)[c]' 'G'
'His-tRNA(His)[c]' 'H'
'Ile-tRNA(Ile)[c]' 'I'
'Leu-tRNA(Leu)[c]' 'L'
'Lys-tRNA(Lys)[c]' 'K'
'Met-tRNA(Met)[c]' 'M'
'Phe-tRNA(Phe)[c]' 'F'
'Pro-tRNA(Pro)[c]' 'P'
'Ser-tRNA(Ser)[c]' 'S'
'Thr-tRNA(Thr)[c]' 'T'
'Trp-tRNA(Trp)[c]' 'W'
'Tyr-tRNA(Tyr)[c]' 'Y'
'Val-tRNA(Val)[c]' 'V'
};

for q = 1:length(AAstruct) %create third column with all zeros
    AAstruct{q,3}=0;
end
    
ATPcost = num2str((length(inputSeq)*4)-1); %approximately 4 ATP/amino acid

for c = 1:length(inputSeq) %cycle through each AA in input sequence
    for d = 1:length(AAstruct) %cycle through each possible AA
       if inputSeq(c) == AAstruct{d,2} %check if match
           AAstruct{d,3}= AAstruct{d,3}+1; %add count if match
       end
    end
end

for w = 1:length(AAstruct) %convert numbers to strings
    AAstruct{w,3} = num2str(AAstruct{w,3});
end

rxtntSide = ''; %initialize vectors to hold eqn
prodSide = ''; 
eqn = ''; 

for r = 1:length(AAstruct) % generate rxtns side of eqation
    rxtntSide = [rxtntSide ' + ' AAstruct{r,3} ' ' AAstruct{r,1}];
    prodSide = [prodSide ' + ' AAstruct{r,3} ' ' AAstruct{r,1}(5:length(AAstruct{r,1}))];
end

eqn = [rxtntSide(4:length(rxtntSide)) ' + ' ATPcost ' ATP[c]' ' + ' ATPcost ' H2O[c] => ' prodSide(4:length(prodSide)) ' + ' ATPcost ' ADP[c]' ' + ' ATPcost ' phosphate[c]' ' + ' ATPcost ' H+[c]' ' + ' proteinName '[c]' ]; %add ATP

end