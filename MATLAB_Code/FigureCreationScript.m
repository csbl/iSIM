%%
clear all
close all
load toycon1.mat; %load model
changeCobraSolver('gurobi'); %change solver (change to whatever solver you are using)

%%
toycon1_rxn_decompositon = struct(); %create structure
smat = full(model.S)  %grab smatrix
[numMet,numReact] = size(smat); %get number of reactions and metabolites
count =1 ; %counter for structure indexing
for j = 1:numReact
   k = find(smat(:,j)); %find indexes of non zero entries in a given column of smatrix (reaction)
   for l = 1:length(k)
        met_id(count) = model.mets(k(l));
        rxn_id(count) = model.rxns(j); %assign characteristics
        coeff(count) = smat(k(l),j);
        temp = char(met_id(count));
        met_compartment(count) = string(temp(5));  
        met_compound(count) = string(temp(1:3));
        cpd_name(count) = model.metNames(k(l));
        met_name(count) = strcat(model.metNames(k(l)), '[' , met_compartment(count) , ']'); 
        count = count + 1;
   end
end

%Create structure fields and observations
toycon1_rxn_decompositon.met_id = met_id;
toycon1_rxn_decompositon.rxn_id = rxn_id;
toycon1_rxn_decompositon.coeff = coeff;
toycon1_rxn_decompositon.met_compartment = met_compartment;
toycon1_rxn_decompositon.met_compount = met_compound;
toycon1_rxn_decompositon.cpd_name = cpd_name;
toycon1_rxn_decompositon.met_name = met_name;
disp(toycon1_rxn_decompositon)

%Output into text file
file = fopen('toycon1_rxn_decompositon.txt','w');
fprintf(file,'%s ',string(fieldnames(toycon1_rxn_decompositon)));
fprintf(file,'\n%s %s %s %s %s %s %s',[string(met_id);string(rxn_id);floor(coeff);string(met_compartment);string(met_compound)...
    ;string(cpd_name);string(met_name)]);
fclose(file);
%% Make S matrices txt files
%output S matrix with names
file = fopen('toycon1_smatrix_names.txt','w');
smatrix_names = [string(model.metNames),smat];
smatrix_names = [" ",string(model.rxnNames)';smatrix_names];
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',smatrix_names');
fclose(file);

%output S matrix with IDs
file = fopen('toycon1_smatrix_id.txt','w');
smatrix_id = [string(model.mets),smat];
smatrix_id = [" ",string(model.rxns)';smatrix_id];
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',smatrix_id');
%% Make rxn info struct
toycon1_rxn_info = struct();
%assign rxn_id,rxn_name,lower and upper bound fields
toycon1_rxn_info.rxn_id = model.rxns;
toycon1_rxn_info.rxn_name = model.rxnNames;
toycon1_rxn_info.lb = model.lb;
toycon1_rxn_info.ub = model.ub;
%construct reactions
for i = 1 :numReact
    sym = "-->";
    if model.rev(i) == 1%determine directionality symbol
        sym = "<==>";
    end
    temp = sym;
    k = find(smat(:,i)); 
    fcount = 0;
    rcount = 0;
    for j = 1 : length(k)
        comp = char(model.mets(k(j)));
        comp  = strcat('[',comp(5),']'); %get compartment
        if smat(k(j),i) < 0 %add onto reactants
            if fcount > 0
                temp = strcat(strcat(string(abs(smat(k(j),i))),model.metNames(k(j)),comp),strcat(' + ',temp));
            else
                temp = strcat(strcat(string(abs(smat(k(j),i))),model.metNames(k(j)),comp),temp);
            end
            fcount = fcount + 1;
        else  %add to products
            if rcount > 0
                temp = strcat(strcat(temp,' + '),strcat(string(abs(smat(k(j),i))),model.metNames(k(j))),comp);
            else
                temp = strcat(temp,strcat(string(abs(smat(k(j),i))),model.metNames(k(j))),comp);
            end
            rcount = rcount + 1;
        end
    end
    reactForm(i) = strrep(temp,'1',''); %remove 1 coefficent and add to vector
end
toycon1_rxn_info.rxn_formula = reactForm'; %assign reaction formula fields
disp(toycon1_rxn_info)
file = fopen('toycon1_rxn_info.txt','w');
fprintf(file,'%s %s %s %s %s\n',[string(toycon1_rxn_info.rxn_id),string(toycon1_rxn_info.rxn_name),string(toycon1_rxn_info.lb),string(toycon1_rxn_info.ub),...
    string(toycon1_rxn_info.rxn_formula)]');
fclose(file);

%% make S matrix pdf
fig = figure;
smat = fliplr(smat)
temp = smat;
smat(:,6:9) = temp(:,1:4);
smat(:,1:5) = temp(:,5:9);
[Posrow,Poscol] = find(smat < 0 );%find negative coefficient matrix indices
ycoord = numMet - Posrow; %transform to xy coordinates
xcoord = numReact - Poscol;
labels = string(smat(smat < 0)); %make coefficent labels
labels = char(labels);
histogram2(xcoord,ycoord,'DisplayStyle','tile','ShowEmptyBins','off');%plot negative coefficient values in histogram
hold on;
text(xcoord,ycoord,labels,'color','white');%add labels
[Posrow,Poscol] = find(smat > 0);%repeat for positive coefficients
ycoord = numMet - Posrow;
xcoord = numReact - Poscol;
labels = char(string(smat(smat > 0)));
histogram2([xcoord;xcoord],[ycoord;ycoord],'DisplayStyle','tile','ShowEmptyBins','off');
text(xcoord,ycoord,labels,'color','white');
map = [1,1,0;0,0,1;1,0,0];%make color map
colormap(map);%apply color map
xlabels = flip(model.rxnNames);
temp = xlabels;
xlabels(6:9) = temp(1:4);
xlabels(1:5) = temp(5:9);
xlabels = flip(xlabels)
xticklabels(xlabels)%create xaxis tick labels
xtickangle(45); %rotate labels
for x = 1:numMet
    temp = char(model.mets{x}) ;
    metCompartments{x} = temp(length(temp)-2:length(temp));%create ylabels
end
yticklabels(strcat(flip(model.metNames),flip(metCompartments')))%add y tick labels
set(fig,'Units','Inches');
pos = get(fig,'Position'); %https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(fig,'toycon1_smatrix','pdf')%save as pdf
smat = fliplr(smat);
%%      Perform flux variabilitiy analysis with percentage of objective function %%%%%%%%

fva_pct_result = ef_tbl_fva(0,model,toycon1_rxn_info,0);%perform initial fva
for i = 1:20
   fva_pct_result = ef_tbl_fva(i*5,model,fva_pct_result,1); %perform for all percentages
end
disp(fva_pct_result)
file = fopen('toycon1_fva_result_percentage.txt','w');%open file
fprintf(file,'rxn_id fva_lb fva_ub rxn_name lb ub rxn_formula fva_pct fva_req fva_on\n')%print headers
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',[string(fva_pct_result.rxn_id),string(fva_pct_result.fva_lb),...%print file
    string(fva_pct_result.fva_ub),string(fva_pct_result.rxn_name),string(fva_pct_result.lb),string(fva_pct_result.ub),...
    string(fva_pct_result.rxn_formula),string(fva_pct_result.fva_pct),string(fva_pct_result.fva_req),string(fva_pct_result.fva_on)]')

fclose(file);%close file

%%
%Plot percentage
fig = figure;
j = 1;
for rx = ["R1","R2"] %reactions to plot
    subplot(1,2,j)
    map = [.2,.2,.2;0,0,1;1,0,0];%make color map
    colormap(map);%apply color map
    k = strcmp(char(fva_pct_result.rxn_id) ,rx); %extract indexes for matching reactions
    xcoordl = fva_pct_result.fva_lb(k) ; %gather bounds for x coordinates
    xcoordu = fva_pct_result.fva_ub(k) ;
    ycoord = fva_pct_result.fva_pct(k); %gather pct for y coordinates
    color = coloring(fva_pct_result.fva_on(k),fva_pct_result.fva_req(k)); %make coloring based on fva_req and fva_off
    scatter(xcoordu,ycoord,[100],color,'filled','d') %Plot points
    hold on
    scatter(xcoordl,ycoord,[100],color,'filled','s')
    for i = 1:length(ycoord)
        switch color(i)
            case 10
                temp = 'red';
            case 0
                temp = 'black';
            case 5
                temp = 'blue';
        end
        line([xcoordl(i),xcoordu(i)],[ycoord(i),ycoord(i)],'color',temp) %Draw connecting lines
    end
    j = j +1;
    xlabel('Units of Flux Through Reaction')
    ylabel('Required Flux Through the Objective Function')%label
    xlim([0,2])  %set xaxis limits
    temp = fva_pct_result.rxn_name(k);
    t = title(string(temp(1)));
    pos = get(t,'Position');
    pos(2) = pos(2) + 3;
    set(t,'Position',pos);
    grid on;
end
%save figure
set(fig,'Units','Inches');
pos = get(fig,'Position'); %https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
saveas(fig,'toycon1_fva_percentage','pdf')

%%
% Peform FVA again but with ATP increments this time
fva_inc_result = ef_tbl_fva(0,model,toycon1_rxn_info,0);%perform initial fva
for i = 1:16
   fva_inc_result = ef_tbl_fva(100*i/16,model,fva_inc_result,1); %perform for all increments
end
disp(fva_inc_result)
file = fopen('toycon1_fva_result_increment.txt','w');%open file
fprintf(file,'rxn_id fva_lb fva_ub rxn_name lb ub rxn_formula fva_pct fva_req fva_on\n')%print headers
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',[string(fva_inc_result.rxn_id),string(fva_inc_result.fva_lb),...%print file
    string(fva_inc_result.fva_ub),string(fva_inc_result.rxn_name),string(fva_inc_result.lb),string(fva_inc_result.ub),...
    string(fva_inc_result.rxn_formula),string(fva_inc_result.fva_pct),string(fva_inc_result.fva_req),string(fva_inc_result.fva_on)]')

fclose(file);%close file

%%
%Plot increment
fig = figure;
j = 1;
for rx = ["R1","R2"] %reactions to plot
    subplot(1,2,j)
    map = [.2,.2,.2;0,0,1;1,0,0];%make color map
    colormap(map);%apply color map
    k = strcmp(char(fva_inc_result.rxn_id) ,rx); %extract indexes for matching reactions
    xcoordl = fva_inc_result.fva_lb(k) ; %gather bounds for x coordinates
    xcoordu = fva_inc_result.fva_ub(k) ;
    ycoord = fva_inc_result.fva_pct(k).*32/100; %gather pct for y coordinates
    color = coloring(fva_inc_result.fva_on(k),fva_inc_result.fva_req(k)); %make coloring based on fva_req and fva_off
    scatter(xcoordu,ycoord,[100],color,'filled','d') %Plot points
    hold on
    scatter(xcoordl,ycoord,[100],color,'filled','s')
    for i = 1:length(ycoord)
        switch color(i)
            case 10
                temp = 'red';
            case 0
                temp = 'black';
            case 5
                temp = 'blue';
        end
        line([xcoordl(i),xcoordu(i)],[ycoord(i),ycoord(i)],'color',temp) %Draw connecting lines
    end
    j = j +1;
    xlabel('Units of Flux Through Reaction')
    ylabel('Required # of ATP produced')%label
    xlim([0,2])  %set xaxis limits
    temp = fva_inc_result.rxn_name(k);
    t = title(string(temp(1)));
    pos = get(t,'Position');
    pos(2) = pos(2) + 1;
    set(t,'Position',pos);
    grid on;
end
%save figure
set(fig,'Units','Inches');
pos = get(fig,'Position'); %https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
saveas(fig,'toycon1_fva_increment','pdf')
%%
% Create figure for all reactions
fig = figure;
j = 1;
rxnlist = string(model.rxns);
for rx = 1:length(rxnlist) %reactions to plot
    ax1 = subplot(3,3,j);
    k = strcmp(char(fva_inc_result.rxn_id) ,rxnlist(rx)); %extract indexes for matching reactions
    xcoordl = fva_inc_result.fva_lb(k) ; %gather bounds for x coordinates
    xcoordu = fva_inc_result.fva_ub(k) ;
    ycoord = fva_inc_result.fva_pct(k).*32/100; %gather pct for y coordinates
    color = coloring(fva_inc_result.fva_on(k),fva_inc_result.fva_req(k)); %make coloring based on fva_req and fva_off
    if max(color) ~= 10
        disp('changed')
        map = [0,0,0;0,1,0;0,0,1];%make color map    
    else
        map = [0,0,0;0,0,1;1,0,0];%make color map    
    end
    colormap(ax1,map);%apply color map 
    scatter(xcoordu,ycoord,25,color,'filled','d') %Plot points
    hold on
    scatter(xcoordl,ycoord,25,color,'filled','s')
    for i = 1:length(ycoord)
        switch color(i)
            case 10
                temp = 'red';
            case 0
                temp = 'black';
            case 5
                temp = 'blue';
        end
        line([xcoordl(i),xcoordu(i)],[ycoord(i),ycoord(i)],'Color',temp) %Draw connecting lines
    end
    j = j +1;
    xlim([min(xcoordl),max(xcoordu)])
    ylim([0,max(ycoord)])
    temp = fva_inc_result.rxn_name(k);        
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times','fontsize',4)
    a = get(gca,'YTickLabel');
    set(gca,'YTickLabel',a,'FontName','Times','fontsize',4)
    xlabel('Units of Flux Through Reaction','FontSize',6)
    ylabel('Required # of ATP produced','FontSize',6)%label
    t = title(string(temp(1)),'FontWeight','normal','FontSize',8);
    pos = get(t,'Position');
    pos(2) = pos(2) ;
    set(t,'Position',pos);
    grid on;
end
%save figure
set(fig,'Units','Inches');
pos = get(fig,'Position'); %https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
saveas(fig,'toycon1_fva_increment_all','pdf')


%%%%%%%Perform Single Gene Deletions %%%%%%
%%
model = buildRxnGeneMat(model);%create RxnGene relationship field in model
[~,~,~,~,~,sol] = singleGeneDeletion(model);%perform single gene simulation
toycon1_gene_ko = struct(toycon1_rxn_info); %create new data struct
toycon1_gene_ko.gene_id = model.grRules; % add gene_ids
dele = zeros(length(model.grRules),1);
for i =1:length(model.grRules)
    k = find(strcmp(model.grRules,model.genes(i)));
    dele(k) = sol(length(sol(:,1)),i);
end
toycon1_gene_ko.gene_ko_atp =dele; %add result of simulation (Flux through objective function)
disp(toycon1_gene_ko)
file = fopen('toycon1_gene_ko_screen.txt','w');%output result
fprintf(file,'gene_id gene_ko_atp rxn_id rxn_name lb ub rxn_formula\n')
fprintf(file,'%s %s %s %s %s %s %s\n',[string(toycon1_gene_ko.gene_id),string(toycon1_gene_ko.gene_ko_atp),string(toycon1_gene_ko.rxn_id),string(toycon1_gene_ko.rxn_name),string(toycon1_gene_ko.lb),string(toycon1_gene_ko.ub),...
    string(toycon1_gene_ko.rxn_formula)]');
fclose(file);
%%
%%%%%%Perform double gene deletion simulation%%%%%%%%%%%%%%%%%%%%%
[~,res,~] = doubleGeneDeletion(model);%perform double gene deletion
rowVal = model.genes;%gather row and col gene names
colVal = model.genes;
[numRows,numCols] = size(res);
z = 1;
for i = 1:numRows
    for j = 1:i-1
        gene1(z) = rowVal(i);
        gene2(z) = rowVal(j);%get corresponding genes deleted
        gene12(z) = strcat(gene1(z),'_',gene2(z));
        k1 = find(strcmp(model.grRules,gene1(z)));%find matching reactions
        k2 = find(strcmp(model.grRules,gene2(z)));
        rxnID1(z) = model.rxns(k1);
        rxnID2(z) = model.rxns(k2);%get matching reaction IDS and names
        rxnID12(z) = strcat(rxnID1(z),'_',rxnID2(z));
        rxnN1(z) = model.rxnNames(k1);
        rxnN2(z) = model.rxnNames(k2);
        rxnN12(z) = strcat(rxnN1(z), ' / ' , rxnN2(z));
        atp1(z) = toycon1_gene_ko.gene_ko_atp(k1);%get corresponding single gene deletion atp result
        atp2(z) = toycon1_gene_ko.gene_ko_atp(k2);
        atp12(z) = res(i,j);%get double deletion result
        z = z + 1;
    end
end

gene_ko2_rxns = struct();%add vectors to structure
gene_ko2_rxns.genes = gene12;
gene_ko2_rxns.rxn1 = rxnID1;
gene_ko2_rxns.rxn2 = rxnID2;
gene_ko2_rxns.name1 = rxnN1;
gene_ko2_rxns.name2 = rxnN2;
gene_ko2_rxns.rxns = rxnID12;
gene_ko2_rxns.names = rxnN12;
disp(gene_ko2_rxns)
file = fopen('toycon1_gene_2ko_rxns.txt','w');%output result
fprintf(file,'genes rxn1 rxn2 name1 name2 rxns names\n');
fprintf(file,'%s %s %s %s %s %s %s\n',[string(gene_ko2_rxns.genes);string(gene_ko2_rxns.rxn1);...
    string(gene_ko2_rxns.rxn2);string(gene_ko2_rxns.name1);string(gene_ko2_rxns.name2);string(gene_ko2_rxns.rxns);string(gene_ko2_rxns.names)]);
fclose(file);

gene_ko2_tbl = struct(gene_ko2_rxns);%create new structure from existing
gene_ko2_tbl.atp = atp12;%add atp results into structure
gene_ko2_tbl.atp1 = atp1;
gene_ko2_tbl.atp2 = atp2;
gene_ko2_tbl.atp12 = atp12;

disp(gene_ko2_tbl)

file = fopen('toycon1_gene_2ko_tbl.txt','w');%output result
fprintf(file,'genes atp rxn1 rxn2 name1 name2 rxns names atp1 atp2 atp12\n');
fprintf(file,'%s %s %s %s %s %s %s %s %s %s %s\n',[string(gene_ko2_tbl.genes);string(gene_ko2_tbl.atp);string(gene_ko2_tbl.rxn1);...
    string(gene_ko2_tbl.rxn2);string(gene_ko2_tbl.name1);string(gene_ko2_tbl.name2);string(gene_ko2_tbl.rxns);string(gene_ko2_tbl.names);...
    string(gene_ko2_tbl.atp1);string(gene_ko2_tbl.atp2);string(gene_ko2_tbl.atp12)]);
fclose(file);

%%
%find double gene deletions where the double deletion causes a larger
%negative effect on growth than both of the constituient single deletions
cond1 = (gene_ko2_tbl.atp1 - gene_ko2_tbl.atp12) > 0;
cond2 = (gene_ko2_tbl.atp2 - gene_ko2_tbl.atp12) > 0;%create logical matrices for the two conditons
z= 1;
for i =1 : length(cond1)
   if cond1(i) && cond2(i)%find indices that match both conditions
       indices(z) = i;
       z = z  + 1;
   end
end
fprintf('genes atp rxn1 rxn2 name1 name2 rxns names atp1 atp2 atp12\n'); %extract filtered results
for i = indices
   fprintf('%s %s %s %s %s %s %s %s %s %s %s\n',[string(gene_ko2_tbl.genes(i));string(gene_ko2_tbl.atp(i));string(gene_ko2_tbl.rxn1(i));...
    string(gene_ko2_tbl.rxn2(i));string(gene_ko2_tbl.name1(i));string(gene_ko2_tbl.name2(i));string(gene_ko2_tbl.rxns(i));string(gene_ko2_tbl.names(i));...
    string(gene_ko2_tbl.atp1(i));string(gene_ko2_tbl.atp2(i));string(gene_ko2_tbl.atp12(i))]);  
end
