%This program creates a cobra model for the toycon1 network. Reactions,
%genes, and metabolites are added into the model, and the model is ouputed
%as a .mat file that is used in the FigureCreationScript.m program that
%creates the figures associated with this model.

%initCobraToolbox  %uncommented the first time running this program
clear all

%Reaction formulas <=> means reversible -> means an irreverisble reaction.
%Metabolite names are formatted as follows: metName[metCompartment] where
%[c] is the cytoplasm and [m] is the mitochondria. 
E1 = '1 m04[c] <=>';
E2 = '1 m07[c] ->';
E3 = '1 m08[c] <=>';
E4 = '1 m06[c] <=>';
E5 = '1 m03[c] <=>';
R1 = '2 m01[c] + 1 m04[c] + 2 m09[c] -> 2 m02[c] + 2 m06[c] + 2 m07[c]';
R2 = '1 m01[c] + 56 m05[m] + 1 m07[c] + 3 m08[c] + 1 m09[c] -> 1 m02[c] + 3 m03[c] + 56 m05[c] + 4 m06[c]';
R3 = '1 m01[c] + 4 m05[c] + 1 m09[c] -> 1 m02[c] + 4 m05[m] + 1 m06[c]';
R4 = '1 m02[c] + 1 m06[c] -> 1 m01[c] + 1 m09[c]';

%Define reaction IDs for the model
reactionIDs = {'E1','E2','E3','E4','E5','R1','R2','R3','R4'};
%define reversibility
rev = [1,0,1,1,1,0,0,0,0];
%Define reaction names
reactionNames = {'glucose exchange','lactate exchange','O2 exchange','H20 exchange','CO2 exchange','glycolysis','respiration','ATP synthase','ATP demand'};
%Gather formulas
reactionForm = {E1,E2,E3,E4,E5,R1,R2,R3,R4};
%Define lower reaction flux bounds
lb = double([-1,0,-1000,-1000,-1000,0,0,0,0]);
%Define upper reaction flux bounds
ub = double([1000,1000,1000,1000,1000,1000,1000,1000,1000]);
%Define gene reaction rules for the above reactions
grRule = {'HK','LDH','ETC','AQP','CO2','PFK','PDH','ATPV','MYH2'}';
%gather gene names
geneNames = {'HK','LDH','ETC','AQP','CO2','PFK','PDH','ATPV','MYH2'}';
%create model with above arguments
model = createModel(reactionIDs,reactionNames,reactionForm);
model.lb = lb';
model.ub = ub';
model.grRules = grRule;
model.genes = geneNames;
model.rev = rev;
%add metabolite names
model.metNames = {'glucose','lactate','O2','H20','CO2','ADP','Pi','ATP','H','H'}';
%set objective reaction which will be optimized
model = changeObjective(model,'R4');
model.c = model.c;
%save model
save('toycon1.mat','model');

%output summary of the model's reactions
disp('--------------------------Summary of Model Reactions------------------------')
disp([{'ReactionName','Reaction ID','Lower Bound','Upper Bound','Formula','Gene Rule','Objective'};...
    model.rxnNames,model.rxns,num2cell(model.lb),num2cell(model.ub),printRxnFormula(model),model.grRules,num2cell(model.c)])