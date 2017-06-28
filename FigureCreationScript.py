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
    result = []
    for y in x:
        if not y == 0:
            result.append(numpy.random.normal(0,1))
        else:
            result.append(0)
    return numpy.random.uniform(0,1,len(x))
def rxnIDCon(x):
    result = []
    value = 0
    i = 0
    for y in x:
        if y.startswith('R') == 1:
            value +=4

toycon1_rxn_decomposition2 = toycon1_rxn_decomposition.copy()
toycon1_rxn_decomposition2['sign'] = pandas.Series(sign(toycon1_rxn_decomposition['coeff'].values),index = toycon1_rxn_decomposition2.index)
rxnUniques,temp=numpy.unique(toycon1_rxn_decomposition['rxn_id'].values,return_inverse=True)
toycon1_rxn_decomposition2['rxn_id'] = pandas.Series(temp.copy(),index = toycon1_rxn_decomposition.index)
metUniques,temp=numpy.unique(toycon1_rxn_decomposition['met_name'].values,return_inverse=True)
toycon1_rxn_decomposition2['met_name'] = pandas.Series(temp.copy(),index = toycon1_rxn_decomposition.index)

print toycon1_rxn_decomposition2['sign']
print toycon1_rxn_decomposition2['rxn_id']
print toycon1_rxn_decomposition2['met_name']

p = ggplot(toycon1_rxn_decomposition2,aes(x='rxn_id',y='met_name',fill='sign'))+geom_tile(fill='sign')
#p = p+geom_text(size = 2,color = "black")+scale_color_gradient(low= 'grey',high = 'blue')
p.show()

df = pandas.DataFrame(dict(
    x=numpy.random.normal(0, 1, 1000),
    y=numpy.random.normal(0, 1, 1000),
    w=numpy.random.uniform(0, 1, 1000)
))
print toycon1_rxn_decomposition2.head()
print df['w']
p = ggplot(df, aes(x='x', y='y', fill='w')) + geom_tile(fill= 'w')
p.show()
#p = p + geom_tile



"""
plotPar = aes(x = 'rxn_name', y = 'met_name', fill = numpy.sign(coef), label = coef)) + geom_tile() +
geom_text(size = 2, color = "#FFFFFF") + scale_fill_gradient2(low = "#4F81BD", mid = "#FFFFFF", high = "#C0504D") +
theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + xlab(NULL) + ylab(NULL)
#ggplot(
numpy.sign"""

