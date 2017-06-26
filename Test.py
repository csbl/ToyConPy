import cobra.io
import cobra.flux_analysis
from cobra import Model, Reaction, Metabolite
import cobra.solvers.gurobi_solver



toyconmodel = Model('ToyCon')
reaction = Reaction('E1')
reaction.name = 'glucose exchange'
reaction.lower_bound = -1
reaction.upper_bound = 1000
reaction.subsystem='Metabolism'

glc = Metabolite('m04c',name= 'glucose',compartment='c')
lct = Metabolite('m07c',name = 'lactate',compartment='c')
O2 = Metabolite('m08c',name = '02',compartment='c')
h20 = Metabolite('m06c',name = 'H20',compartment='c')
CO2 = Metabolite('m03c',name = 'CO2',compartment='c')
adp = Metabolite('m01c',name = 'CO2',compartment='c')
pi = Metabolite('m09c',name = 'Pi',compartment='c')
atp = Metabolite('m02c',name = 'ATP',compartment='c')
Hm = Metabolite('m05m',name = 'H',compartment='m')
Hc = Metabolite('m05c',name = 'H',compartment='c')


reaction.add_metabolites({
    glc : -1,
})


reaction.gene_reaction_rule = '( HK )'

toyconmodel.add_reactions(reaction)


reaction = Reaction('E2','lactate exchange',subsystem='Metabolism',lower_bound=0,upper_bound=1000)
reaction.add_metabolites({lct:-1})
reaction.gene_reaction_rule = '( LDH )'
toyconmodel.add_reactions(reaction)

reaction = Reaction('E3','O2 exchange',subsystem='Metabolism',lower_bound=-1000,upper_bound=1000)
reaction.add_metabolites({O2:-1})
reaction.gene_reaction_rule = '( ETC )'
toyconmodel.add_reactions(reaction)

reaction = Reaction('E4','H20 exchange',subsystem='Metabolism',lower_bound=-1000,upper_bound=1000)
reaction.add_metabolites({h20:-1})
reaction.gene_reaction_rule = '( AQP )'
toyconmodel.add_reactions(reaction)

reaction = Reaction('E5','CO2 exchange',subsystem='Metabolism',lower_bound=-1000,upper_bound=1000)
reaction.add_metabolites({CO2:-1})
reaction.gene_reaction_rule = '( CO2 )'
toyconmodel.add_reactions(reaction)

reaction = Reaction('R1','glycolysis',subsystem='Metabolism',lower_bound=0,upper_bound=1000)
reaction.add_metabolites({glc:-1,lct:2,h20:2,adp:-2,pi:-2,atp:2})
reaction.gene_reaction_rule = '( PFK )'
toyconmodel.add_reactions(reaction)

reaction = Reaction('R2','respiration',subsystem='Metabolism',lower_bound=0,upper_bound=1000)
reaction.add_metabolites({lct:-1,O2:-3,h20:4,CO2:3,adp:-1,pi:-1,atp:1,Hm:-56,Hc:56})
reaction.gene_reaction_rule = '( PDH )'
toyconmodel.add_reactions(reaction)

reaction = Reaction('R3','ATP synthase',subsystem='Metabolism',lower_bound=0,upper_bound=1000)
reaction.add_metabolites({h20:1,adp:-1,pi:-1,atp:1,Hm:4,Hc:-4})
reaction.gene_reaction_rule = '( ATPV )'
toyconmodel.add_reactions(reaction)

reaction = Reaction('R4','ATP demand',subsystem='Metabolism',lower_bound=0,upper_bound=1000)
reaction.add_metabolites({h20:-1,adp:1,pi:1,atp:-1})
reaction.gene_reaction_rule = '( MYH2 )'
toyconmodel.add_reactions(reaction)

toyconmodel.objective = 'R4'

toyconmodel.solver = 'gurobi'

# Iterate through the the objects in the toyconmodel
print("Reactions")
print("---------")
for x in toyconmodel.reactions:
    print("%s : %s" % (x.id, x.reaction))

print("")
print("Metabolites")
print("-----------")
for x in toyconmodel.metabolites:
    print('%9s : %s' % (x.id, x.formula))

print("")
print("Genes")
print("-----")
for x in toyconmodel.genes:
    associated_ids = (i.id for i in x.reactions)
    print("%s is associated with reactions: %s" %
          (x.id, "{" + ", ".join(associated_ids) + "}"))

print(toyconmodel.objective.expression)
print(toyconmodel.objective.direction)


toyconmodel.optimize()

toyconmodel.summary()

