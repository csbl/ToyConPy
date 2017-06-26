import cobra.io
import cobra.core.arraybasedmodel
import numpy
import pandas

model = cobra.io.read_sbml_model('ToyCon.xml')
model.solver = 'gurobi'
model.optimize()


Amodel = model.to_array_based_model()
smatrix = numpy.matrix(Amodel.S)


sdata = pandas.DataFrame(Amodel.S.toarray(),[x.name for x in model.metabolites],[x.name for x in model.reactions])
pandas.DataFrame.to_csv(sdata,'SMatrix.csv')
#print sdata