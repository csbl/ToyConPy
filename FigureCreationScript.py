import cobra.io
import cobra.core.arraybasedmodel
import numpy
import pandas
import csv
import ggplot
from ggplot import *
import matplotlib.pyplot as plt
import matplotlib.axes
import matplotlib.colors
from matplotlib.backends.backend_pdf import PdfPages

model = cobra.io.read_sbml_model('ToyCon.xml')#Read in toycon model
model.solver = 'gurobi' #set solver
model.optimize() #optimize model

Amodel = model.to_array_based_model() #change to arraybased model to extract S matrix
smatrix = numpy.matrix(Amodel.S) #extract S matrix

#make S matrix with names
sdataName = pandas.DataFrame(Amodel.S.toarray(),[x.name for x in model.metabolites],[x.name for x in model.reactions])#create data frame object
pandas.DataFrame.to_csv(sdataName,'toycon1_smatrix_names.txt',index = True,sep = ' ',quoting=csv.QUOTE_NONE,header=True,escapechar=' ',float_format='%d') #print S matrix as txt file

#make S matrix with IDs
sdataID = pandas.DataFrame(Amodel.S.toarray(),[x.id for x in model.metabolites],[x.id for x in model.reactions])#create data frame object
pandas.DataFrame.to_csv(sdataID,'toycon1_smatrix_id.txt',index = True,sep = ' ',quoting=csv.QUOTE_NONE,header=True,escapechar=' ',float_format='%d') #print S matrix as txt file

#creat a reaction decomposition of the S matrix where one reaction is represented as multiple rows
rxnDecompDict = {}
z = 0
for x in model.reactions:
    for y in x.metabolites: #create dictionary
        rxnDecompDict[z] = str.split('%s %s %d %s %s %s %s' % (y.id,x.id,x.get_coefficient(y.id),y.compartment,y.id[:-1],y.name,y.name+'['+y.compartment+']'),' ')
        rxnDecompDict[z][2] = int(rxnDecompDict[z][2])
        z += 1

toycon1_rxn_decomposition = pandas.DataFrame.from_dict(rxnDecompDict,orient='index') #create datafrane
toycon1_rxn_decomposition.columns=['met_id','rxn_id','coeff','met_compartment','met_compound','cpd_name','met_name']#add column labels
toycon1_rxn_decomposition.to_csv('toycon1_rxn_decomposition.txt',sep = ' ',float_format='%d')#output txt file
print 'The first 6 entries of the toycon1_rxn_decomposition'
print toycon1_rxn_decomposition.head(6) #print the first few rows

#create reaction formulas
rxnInfoDict = {} #create new dictionary
z=0
for x in model.reactions:
    temp = str.replace(x.build_reaction_string(True),'.0','')
    temp = temp.replace('-','=')
    rxnInfoDict[z] = [x.id,x.name,x.lower_bound,x.upper_bound,temp] #add key value pairs
    z+=1

toycon1_rxn_info = pandas.DataFrame.from_dict(rxnInfoDict,orient = 'index') #create data frame
toycon1_rxn_info.columns = ['rxn_id','rxn_name','lb','ub','rxn_formula'] #write columns
toycon1_rxn_info.to_csv('toycon1_rxn_info.txt',sep = ' ',float_format='%d',quoting=csv.QUOTE_NONE,escapechar=' ')#output formatted txt file
print 'The first 6 entries of the toycon1_rxn_info'
print toycon1_rxn_info.head(6) #output first 6 lines

########################Perform FBA###############################
print 'Objective Reaction is: '
print toycon1_rxn_info[toycon1_rxn_info.rxn_id == 'R4'].rxn_id + '  ' + toycon1_rxn_info[toycon1_rxn_info.rxn_id == 'R4'].rxn_formula##Print objective reaction

#### Function to transform positive and negative coefficient values into 1 for positive and .5 for negative.
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
#toycon1_rxn_decomposition.sort_values(by=[])
#create indexing for rxn_id and met_name
rxnUniques,indexRxn,invRxn=numpy.unique(toycon1_rxn_decomposition['rxn_id'].values,return_inverse=True,return_index=True,)
metUniques,indexMet,invMet=numpy.unique(toycon1_rxn_decomposition['met_name'].values,return_inverse=True,return_index=True)

#Create dictionaries for formatting table
rxnID2Name = pandas.Series(toycon1_rxn_info.rxn_name.values,index = toycon1_rxn_info.rxn_id).to_dict()
metName2int = pandas.Series([(len(model.metabolites)-1) - x for x in range(len(model.metabolites))],index = [x.name+'['+x.compartment+']' for x in model.metabolites]).to_dict()
int2metName = pandas.Series([x.name+'['+x.compartment+']' for x in model.metabolites],index = [(len(model.metabolites)-1) - x for x in range(len(model.metabolites))]).to_dict()

#create y coordinates for table from metabolic names
yCor = [metName2int[x] for x in toycon1_rxn_decomposition['met_name'].values]
#set coloring scheme
coloring = {.5 : 'b',1:'r'}

#create data
toycon1_rxn_decomposition2 = pandas.DataFrame({
    'w' : sign(toycon1_rxn_decomposition['coeff'].values),
    'y' : yCor,
    'x' : invRxn,
    'l' : toycon1_rxn_decomposition['coeff'].values
})

fig = plt.figure()  ##Create S matrix figure and save as a pdf
counts, xedge, yedge, imag = plt.hist2d(bins = [len(rxnUniques),len(metUniques)],x=toycon1_rxn_decomposition2['x'].values,y=toycon1_rxn_decomposition2['y'].values,normed = False,weights = toycon1_rxn_decomposition2['w'].values,cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom",[(0,'White'),(.5,'Blue'),(1,'Red')]))
plt.yticks(range(len(metUniques)),[int2metName[x] for x in range(len(metUniques))])
plt.xticks(range(len(rxnUniques)),[rxnID2Name[x] for x in rxnUniques],rotation = 45)
[plt.text(xedge[x]+.5,yedge[y]+.5,l,color = 'White',ha = 'center', va = 'center') for x,y,l in zip(toycon1_rxn_decomposition2['x'].values,toycon1_rxn_decomposition2['y'].values,toycon1_rxn_decomposition2['l'].values)]
plt.subplots_adjust(bottom = .25)
plt.tick_params(axis = u'both',which = u'both',length = 0)
#plt.show() #TODO: uncomment
pp = PdfPages('toycon1_smatrix.pdf')
pp.savefig(fig)
pp.close()
##To line 164 on R

#function to run FVA with inputs of a required percentage, model, and table to add results

def ef_tbl_fva(pct,model,tbl = None):
    tbl = tbl.copy()
    fvaRes = cobra.flux_analysis.flux_variability_analysis(model,model.reactions[:],fraction_of_optimum=pct)
    ub = fvaRes.maximum.values
    lb = fvaRes.minimum.values
    n = len(ub)
    tbl.insert(1,'fva_lb', ub)
    tbl.insert(2,'fva_ub', lb)
    tbl['fva_pct'] = list([int(z+pct*100) for z in numpy.zeros((n,1),numpy.int64)])
    def fva_req(u,l):
        result = []
        for x,y in zip(u,l):
            if y > 1*10**-9 or x < -1*10**-9:
                result.append(True)
            else:
                result.append(False)
        return result
    tbl['fva_req'] = fva_req(ub,lb)
    def fva_on(u, l):
        result = []
        for x, y in zip(u, l):
            if abs(y) > 1 * 10 ** -9 or abs(x) > 1 * 10 ** -9:
                result.append(True)
            else:
                result.append(False)
        return result
    tbl['fva_on'] = fva_on(ub,lb)
    return tbl

fva_pct_result = ef_tbl_fva(0,model,toycon1_rxn_info)
for x in [(y+1)*.5 for y in range(20)]:
    fva_pct_result = fva_pct_result.append(ef_tbl_fva(0,model,toycon1_rxn_info),ignore_index = True)
fva_pct_result = fva_pct_result.sort_values(by = 'rxn_id')
fva_pct_result = fva_pct_result.reset_index(drop=True)
print fva_pct_result
fva_pct_result.to_csv('toycon1_fva_result_percentage.txt',float_format='%d', quoting = csv.QUOTE_NONE,escapechar=' ',index = False)#output fva result

