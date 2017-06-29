import cobra.io
import cobra.core.arraybasedmodel
import numpy
import pandas
import csv
import ggplot
from ggplot import *

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
print toycon1_rxn_info[toycon1_rxn_info.rxn_id == 'R4'].rxn_id + '  ' + toycon1_rxn_info[toycon1_rxn_info.rxn_id == 'R4'].rxn_formula

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
def rxnIDCon(x):
    result = []
    value = 0
    i = 0
    for y in x:
        if y.startswith('R') == 1:
            value +=4


rxnUniques,temp=numpy.unique(toycon1_rxn_decomposition['rxn_id'].values,return_inverse=True)
metUniques,temp2=numpy.unique(toycon1_rxn_decomposition['met_name'].values,return_inverse=True)
rxnID2Name = pandas.Series(toycon1_rxn_info.rxn_name.values,index = toycon1_rxn_info.rxn_id).to_dict()
print rxnID2Name
print rxnUniques

toycon1_rxn_decomposition2 = pandas.DataFrame({
    'w' : sign(toycon1_rxn_decomposition['coeff'].values),
    'y' : temp2,
    'x' : temp,
    'l' : toycon1_rxn_decomposition['coeff'].values
})


p = ggplot(toycon1_rxn_decomposition2,aes(x='x',y='y',fill='w',label = 'l'))+geom_tile(xbins = len(rxnUniques)+1, ybins = len(metUniques)+1,fill='w')
p = p + xlim(0,max(toycon1_rxn_decomposition2['x'])) + scale_x_continuous(breaks = range(len(rxnUniques)),labels = list(rxnUniques))
p = p + ylim(0,max(toycon1_rxn_decomposition2['y'])) + scale_y_continuous(breaks = range(len(metUniques)),labels = list(metUniques))
p = p + xlab(' ') + ylab(' ') +theme(x_axis_text=element_text(angle=45,hjust=1))
p= p #+ scale_fill_identity()
p.xtick_labels = [rxnID2Name[x] for x in rxnUniques]
#p.apply_scales()
p.xbreaks = range(len(rxnUniques))
p.apply_axis_labels
print p




"""
plotPar = aes(x = 'rxn_name', y = 'met_name', fill = numpy.sign(coef), label = coef)) + geom_tile() +
geom_text(size = 2, color = "#FFFFFF") + scale_fill_gradient2(low = "#4F81BD", mid = "#FFFFFF", high = "#C0504D") +
theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + xlab(NULL) + ylab(NULL)
#ggplot(
numpy.sign"""

