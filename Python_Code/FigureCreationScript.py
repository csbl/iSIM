import cobra.io
import numpy
import pandas
import csv
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages


model = cobra.io.read_sbml_model('iSIM.xml')#Read in iSIM model
model.optimize() #optimize model
Amodel = model.to_array_based_model() #change to arraybased model to extract S matrix
smatrix = numpy.matrix(Amodel.S) #extract S matrix
sdataName = pandas.DataFrame(Amodel.S.toarray(),[x.name for x in model.metabolites],[x.name for x in model.reactions])#create data frame object
pandas.DataFrame.to_csv(sdataName,'iSIM_smatrix_names.txt',index = True,sep = ' ',quoting=csv.QUOTE_NONE,header=True,escapechar=' ',float_format='%d') #print S matrix as txt file

sdataID = pandas.DataFrame(Amodel.S.toarray(),[x.id for x in model.metabolites],[x.id for x in model.reactions])#create data frame object
pandas.DataFrame.to_csv(sdataID,'iSIM_smatrix_id.txt',index = True,sep = ' ',quoting=csv.QUOTE_NONE,header=True,escapechar=' ',float_format='%d') #print S matrix as txt file
rxnDecompDict = {}
z = 0
for x in model.reactions:
    for y in x.metabolites: #create dictionary
        rxnDecompDict[z] = str.split('%s %s %d %s %s %s %s' % (y.id,x.id,x.get_coefficient(y.id),y.compartment,y.id[:-1],y.name,y.name+'['+y.compartment+']'),' ')
        rxnDecompDict[z][2] = int(rxnDecompDict[z][2])
        z += 1

iSIM_rxn_decomposition = pandas.DataFrame.from_dict(rxnDecompDict,orient='index') #create datafrane
iSIM_rxn_decomposition.columns=['met_id','rxn_id','coeff','met_compartment','met_compound','cpd_name','met_name']#add column labels
iSIM_rxn_decomposition.to_csv('iSIM_rxn_decomposition.txt',sep = ' ',float_format='%d')#output txt file
print ('The first 6 entries of the iSIM_rxn_decomposition')
print (iSIM_rxn_decomposition.head(6)) #print the first few rows
rxnInfoDict = {} #create new dictionary
z=0
for x in model.reactions:
    temp = str.replace(x.build_reaction_string(True),'.0','')
    temp = temp.replace('-','=')
    rxnInfoDict[z] = [x.id,x.name,x.lower_bound,x.upper_bound,temp] #add key value pairs
    z+=1
iSIM_rxn_info = pandas.DataFrame.from_dict(rxnInfoDict,orient = 'index') #create data frame
iSIM_rxn_info.columns = ['rxn_id','rxn_name','lb','ub','rxn_formula'] #write columns
iSIM_rxn_info.to_csv('iSIM_rxn_info.txt',sep = ' ',float_format='%d',quoting=csv.QUOTE_NONE,escapechar=' ')#output formatted txt file
print ('The first 6 entries of the iSIM_rxn_info')
print (iSIM_rxn_info.head(6)) #output first 6 lines
print ('Objective Reaction is: ')
print (iSIM_rxn_info[iSIM_rxn_info.rxn_id == 'R4'].rxn_id + '  ' + iSIM_rxn_info[iSIM_rxn_info.rxn_id == 'R4'].rxn_formula) ##Print objective reaction
def sign(x):
    ma = max(x)
    mi = min(x)
    result = []
    for y in x:
        if(y > 0):
            result.append(1)
        else:
            result.append(.5)
    return result
rxnUniques,indexRxn,invRxn=numpy.unique(iSIM_rxn_decomposition['rxn_id'].values,return_inverse=True,return_index=True,)
metUniques,indexMet,invMet=numpy.unique(iSIM_rxn_decomposition['met_name'].values,return_inverse=True,return_index=True)

#Create dictionaries for formatting table
rxnID2Name = pandas.Series(iSIM_rxn_info.rxn_name.values,index = iSIM_rxn_info.rxn_id).to_dict()
metName2int = pandas.Series([(len(model.metabolites)-1) - x for x in range(len(model.metabolites))],index = [x.name+'['+x.compartment+']' for x in model.metabolites]).to_dict()
int2metName = pandas.Series([x.name+'['+x.compartment+']' for x in model.metabolites],index = [(len(model.metabolites)-1) - x for x in range(len(model.metabolites))]).to_dict()

#create y coordinates for table from metabolic names
yCor = [metName2int[x] for x in iSIM_rxn_decomposition['met_name'].values]
#set coloring scheme
coloring = {.5 : 'dodgerblue',1:'firebrick'}
temp = invRxn
for i in range(len(invRxn)):
    if invRxn[i] < 5:
        invRxn[i] += 4
    else:
        invRxn[i]  -= 5
temp = list(rxnUniques)
rxnUniques = list(rxnUniques)
rxnUniques[0:4] = temp[5:9]
rxnUniques[4:9] = temp[0:5]

#create data
iSIM_rxn_decomposition2 = pandas.DataFrame({
    'w' : sign(iSIM_rxn_decomposition['coeff'].values),
    'y' : yCor,
    'x' : invRxn,
    'l' : iSIM_rxn_decomposition['coeff'].values
})
fig = plt.figure()  ##Create S matrix figure and save as a pdf
#Create 2D histogram from metabolites and reactions
counts, xedge, yedge, imag = plt.hist2d(bins = [len(rxnUniques),len(metUniques)],x=iSIM_rxn_decomposition2['x'].values,y=iSIM_rxn_decomposition2['y'].values,normed = False,weights = iSIM_rxn_decomposition2['w'].values,cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom",[(0,'White'),(.5,'dodgerblue'),(1,'firebrick')]))
#format axis labels and ticks
plt.yticks(range(len(metUniques)),[int2metName[x] for x in range(len(metUniques))])
plt.xticks([x * .95 for x in range(len(rxnUniques))],[rxnID2Name[x] for x in rxnUniques],rotation = 50)

#add the coefficient as the label
[plt.text(xedge[x]+.5,yedge[y]+.5,l,color = 'White',ha = 'center', va = 'center') for x,y,l in zip(iSIM_rxn_decomposition2['x'].values,iSIM_rxn_decomposition2['y'].values,iSIM_rxn_decomposition2['l'].values)]
plt.subplots_adjust(bottom = .25)
plt.tick_params(axis = u'both',which = u'both',length = 0)
pp = PdfPages('iSIM_smatrix.pdf') #save figure
pp.savefig(fig)
pp.close()
def ef_tbl_fva(pct,model,tbl = None):
    tbl = tbl.copy()
    fvaRes = cobra.flux_analysis.flux_variability_analysis(model,model.reactions[:],fraction_of_optimum=pct) #perform FVA
    df = pandas.DataFrame(data = fvaRes)
    df = df.transpose()
    ub = list(df["maximum"])
    lb = list(df["minimum"]) #gather bounds
    n = len(ub)
    tbl.insert(1,'fva_lb', lb)
    tbl.insert(2,'fva_ub', ub) #insert bounds into the dataframe
    tbl['fva_pct'] = list([int(z+pct*100) for z in numpy.zeros((n,1),numpy.int64)]) #insert percentage of obj. fun. into dataframe
    def fva_req(u,l): #Function to determine if FVA is required(if one of the upper or lower bounds in either greater than 10^9 or -10^9 respectivly
        result = []
        for x,y in zip(u,l):
            if y > 1*10**-9 or x < -1*10**-9:
                result.append(True)
            else:
                result.append(False)
        return result
    tbl['fva_req'] = fva_req(ub,lb) #insert column with the determination
    def fva_on(u, l): #determine if FVA is on (one of the upper or lower bound must have an absolute value greater than 10^-9
        result = []
        for x, y in zip(u, l):
            if abs(y) > 1 * 10 ** -9 or abs(x) > 1 * 10 ** -9:
                result.append(True)
            else:
                result.append(False)
        return result
    tbl['fva_on'] = fva_on(ub,lb)
    return tbl

fva_pct_result = ef_tbl_fva(0,model,iSIM_rxn_info)#perform unconstrained FBA for initial point
for x in [(y+1)*.05 for y in range(20)]:
    fva_pct_result = fva_pct_result.append(ef_tbl_fva(x,model,iSIM_rxn_info),ignore_index = True) #perform FVA for 5-100% of obj. fun.
fva_pct_result = fva_pct_result.sort_values(by = 'rxn_id') #sort results by rxn_id
fva_pct_result = fva_pct_result.reset_index(drop=True)
print(fva_pct_result.head()) #print first few lines
#output result to text file
fva_pct_result.to_csv('iSIM_fva_result_percentage.txt',sep= ' ',float_format='%.2f', quoting = csv.QUOTE_NONE,escapechar=' ',index = False)#output fva result
def coloring(on,req):
    if on and req:
        return "firebrick"
    if not req and on:
        return "Grey"
    else:
        return "dodgerblue"
fig = plt.figure()
toPlot = ['R1','R2']#the two reaction to plot the FVA result dependency on the percentage of the objective function
i =1
for z in toPlot:
    plt.subplot(1,2,i)
    req = fva_pct_result.loc[fva_pct_result['rxn_id'] == z]['fva_req'].values#grab fva_req value
    on = fva_pct_result.loc[fva_pct_result['rxn_id'] == z]['fva_on'].values#grab fva_on values
    xcoordl = fva_pct_result.loc[fva_pct_result['rxn_id'] == z]['fva_lb'].values
    xcoordu = fva_pct_result.loc[fva_pct_result['rxn_id'] == z]['fva_ub'].values#Grab upper, lower bounds,  and perctage
    ycoord = fva_pct_result.loc[fva_pct_result['rxn_id'] == z]['fva_pct'].values
    plt.xlabel("Units of Flux Through Reaction")
    plt.ylabel("Required Flux Through Obj. Fun.") #plot results
    plt.xlim(-0.2,2.2)
    plt.title(rxnID2Name[z])
    plt.scatter(xcoordl,ycoord,color = [coloring(x,y) for x,y in zip(on,req)],marker="s")
    plt.scatter(xcoordu,ycoord,color = [coloring(x,y) for x,y in zip(on,req)],marker = "D")
    plt.hlines(ycoord,xcoordl,xcoordu,color = [coloring(x,y) for x,y in zip(on,req)])
    i +=1
fig.tight_layout()#adjust layout
pp = PdfPages('iSIM_fva_percentage.pdf')#save figure
pp.savefig(fig)
pp.close()
fva_inc_result = ef_tbl_fva(0,model,iSIM_rxn_info) #perform initial simulation
for x in [(y+1)/16. for y in range(16)]:
    fva_inc_result = fva_inc_result.append(ef_tbl_fva(x,model,iSIM_rxn_info),ignore_index = True) #perform rest of simulation incrementing by 2 atp
fva_inc_result = fva_inc_result.sort_values(by = 'rxn_id')
fva_inc_result = fva_inc_result.reset_index(drop=True) #sort results
print(fva_inc_result.head()) #print first few lines
#save results
fva_inc_result.to_csv('iSIM_fva_result_increment.txt',sep= ' ',float_format='%.2f', quoting = csv.QUOTE_NONE,escapechar=' ',index = False)#output fva result
fig = plt.figure()
toPlot = ['R1','R2']
i =1
for z in toPlot:
    plt.subplot(1,2,i)
    req = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_req'].values
    on = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_on'].values#grab values as above
    xcoordl = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_lb'].values
    xcoordu = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_ub'].values
    ycoord = [y*32/100. for y in fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_pct'].values] #plot results
    plt.xlim(-0.2,2.2)
    plt.xlabel("Units of Flux Through Reaction")
    plt.ylabel("Required # of ATP Produced.")
    plt.title(rxnID2Name[z])
    plt.scatter(xcoordl,ycoord,color = [coloring(x,y) for x,y in zip(on,req)],marker="s")
    plt.scatter(xcoordu,ycoord,color = [coloring(x,y) for x,y in zip(on,req)],marker = "D")
    plt.hlines(ycoord,xcoordl,xcoordu,color = [coloring(x,y) for x,y in zip(on,req)])
    i +=1
fig.tight_layout()
pp = PdfPages('iSIM_fva_increment.pdf')
pp.savefig(fig) #save figure
pp.close()

#create fva_inc plot for all reactions
fig = plt.figure()
toPlot = [y.id for y in model.reactions]#grab all reactions
i =1
matplotlib.rc('xtick', labelsize=4)
matplotlib.rc('ytick', labelsize=4)
for z in toPlot:
    plt.subplot(3,3,i)
    #grab appropriate data, and plot the results
    req = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_req'].values
    on = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_on'].values
    xcoordl = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_lb'].values
    xcoordu = fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_ub'].values
    ycoord = [y*32/100. for y in fva_inc_result.loc[fva_inc_result['rxn_id'] == z]['fva_pct'].values]
    plt.title(rxnID2Name[z],fontsize = 10)
    plt.xlabel("Units of Flux Through Reaction",fontsize  = 5)
    plt.ylabel("Required # of ATP Produced.",fontsize = 5)
    plt.scatter(xcoordl,ycoord,color = [coloring(x,y) for x,y in zip(on,req)],marker="s",s = 3)
    plt.scatter(xcoordu,ycoord,color = [coloring(x,y) for x,y in zip(on,req)],marker = "D",s = 3)
    plt.hlines(ycoord,xcoordl,xcoordu,color = [coloring(x,y) for x,y in zip(on,req)],linewidths = .2)
    i +=1
fig.tight_layout()#format layout
pp = PdfPages('iSIM_fva_increment_all.pdf')#save result
pp.savefig(fig)
pp.close()
res = cobra.flux_analysis.single_gene_deletion(model)
df2 = pandas.DataFrame(list(res))
df2 = df2.transpose()
rxnIDs = []
for x in list(df2.index):
    for y in model.reactions:
        if y.gene_name_reaction_rule == x:
            rxnIDs.append(y.id) #match reactions and flux changes, to create same ordering
test = {'rxn_id': rxnIDs}
res = pandas.DataFrame(list(res))
res = res.transpose()
res.insert(2,'rxn_id',rxnIDs) #insert matching rxn IDs
res = res.sort_values('rxn_id') #sort values
iSIM_gene_ko = iSIM_rxn_info.copy() #make copy of rxn_info dataframe
iSIM_gene_ko.insert(0,'gene_ko_atp',list(res[0].values)) #insert results into dataframe
iSIM_gene_ko.insert(0,'gene_id',list(res.index))
print(iSIM_gene_ko) #print dataframe
#save result
iSIM_gene_ko.to_csv('iSIM_gene_knockout_screen.txt',sep= ' ',float_format='%d', quoting = csv.QUOTE_NONE,escapechar=' ',index = False)#output ko result
res2ko = cobra.flux_analysis.double_gene_deletion(model,return_frame=True,number_of_processes = 1)#perform double gene deletions (# of processes is dependent on system's cpu count
print(res2ko) #print result
g1 = res2ko.index.values.tolist()#make gene lists
g2 = res2ko.columns.values.tolist()
gene_ko2_rxnsT = pandas.DataFrame(columns = ['genes','rxn1','rxn2','name1','name2','rxns','names','atp'])#create new data frame
i = 0
for x in range(len(g1)):
    for y in range(x):
        temp = list()
        temp.append(g1[x]+'_'+g2[y])#gather genes deleted
        for z in model.reactions:
            if z.gene_name_reaction_rule == g2[y] or z.gene_name_reaction_rule == g1[x]:
                temp.append(z.id) #get matching reaction
        for z in model.reactions:
            if z.gene_name_reaction_rule == g2[y] or z.gene_name_reaction_rule == g1[x]:
                temp.append(z.name) #get matching reaction
        temp.append(temp[1]+'_'+temp[2]) #append combined gene names
        temp.append(temp[3]+' / '+temp[4]) #append combines rxn names
        temp.append(round(res2ko.iloc[x,y],1)) #append result of deletion
        gene_ko2_rxnsT.loc[i] = temp #add to dataframe
        i += 1
gene_ko2_rxns = gene_ko2_rxnsT.copy().drop('atp',1).sort_values("genes").reset_index(drop=True)#sort values
gene_ko2_rxns.to_csv('iSIM_gene_2ko_rxns.txt',sep= ' ',float_format='%d', quoting = csv.QUOTE_NONE,escapechar=' ',index = False)#output ko result
print(gene_ko2_rxns.head())#print first few lines of dataframe
gene_ko2_tbl = gene_ko2_rxnsT.copy().drop('atp',1)
gene_ko2_tbl.insert(1,'atp',gene_ko2_rxnsT['atp'].values)#move ordering of dataframe to make new dataframe
atp1 = []
atp2 = []
for x,y in zip(gene_ko2_tbl['rxn1'].values,gene_ko2_tbl['rxn2'].values):
    atp1.append(iSIM_gene_ko.loc[iSIM_gene_ko["rxn_id"]== x]['gene_ko_atp'].values[0])#find matching atp affect for individual ko
    atp2.append(iSIM_gene_ko.loc[iSIM_gene_ko["rxn_id"]== y]['gene_ko_atp'].values[0])
gene_ko2_tbl.insert(8,"atp1",atp1)
gene_ko2_tbl.insert(9,"atp2",atp2)#insert results
gene_ko2_tbl.insert(10,"atp12",gene_ko2_tbl['atp'].values) #insert double ko results
gene_ko2_tbl = gene_ko2_tbl.sort_values("rxns").reset_index(drop = True) #sort dataframe
print(gene_ko2_tbl.head()) #print first few liens
gene_ko2_tbl.to_csv('iSIM_gene_2ko_tbl.txt',sep= ' ',float_format='%d', quoting = csv.QUOTE_NONE,escapechar= ' ',index = False) #save the result
filtered_ko2 = gene_ko2_tbl.query('atp12 < atp1 and atp12 < atp2').reset_index(drop=True) #filter results to show reactions where the double deletion had a greater effect than either single deletion
print(filtered_ko2) #print filtered results
plt.show() #display figures

