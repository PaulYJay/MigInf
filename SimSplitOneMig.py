import msprime
import numpy as np
import ast

def LCA(tree,pop1,pop2):
	#d=sorted(ast.literal_eval(str(tree)))
	d=ast.literal_eval(str(tree))
	r={}
	ind=pop1+pop2
	for i in range(0,pop1):
		d[i]="A"
	for i in range(pop1,ind):
		d[i]="B"
	for key in sorted(d.keys()):
		if key > (ind-1) and d[key] != -5:
			nod=tree.get_children(key)
			if d[nod[0]] == d[nod[1]] and d[nod[0]] !=-5:
				d[key]=d[nod[0]]
				del d[nod[0]]
				del d[nod[1]]
			else:
				r[key]=d[key]
				for k in d.keys():
					if k > key:
						d[k]=-5
	#print(tree.get_time(min(r)))
simul=msprime.simulate(sample_size=20, Ne=10)
tree=next(simul.trees())                                                                                                                    
LCA(tree,10,10)

def sim1():
	TS=6e6 #concedering 3 generation per years and a split 2 millions years ago
	TM=2e6 # Time of start of migration (backward : if O, mean speciation until now)
	TF=6e6 # Time of end of migration
	NA=2e5 # ancestral pop size
	N1=1e5 # pop 1 size
	N2=1e5 # pop 2 size
	S1=10 # Number of sampled individual in pop1
	S2=10# Number of sampled individual in pop2
	M1=0.001 # Migration value
	d=3
	Le=100000 # Length of the sequence
	Mu=2.9e-9 # mutation rate 
	Ru=2.6e-8 # recombination rate. Value for insect such Drosophila

	t.write("#"+"\t"+str(S1)+"\t"+str(S2)+"\t"+str(Le)+"\n") # write to the output the input value of the simul

	population_configurations = [
	msprime.PopulationConfiguration( sample_size=S1, initial_size=N1),
	msprime.PopulationConfiguration( sample_size=S2, initial_size=N2)
	]

	#migration_matrix = [
	#	[0, M1],
	#	[M1, 0],
	#	]

	demographic_events = [
	msprime.MigrationRateChange(time=TM, rate=M1, matrix_index=(0,1)),
	msprime.MigrationRateChange(time=TM, rate=M1, matrix_index=(1,0)),
	msprime.MigrationRateChange(time=TF, rate=0),
	msprime.MassMigration(time=TS, source=0, destination=1, proportion=1.0),
	msprime.PopulationParametersChange(time=TS, initial_size=NA, growth_rate=0, population_id=1)
	]

	#dp=msprime.DemographyDebugger(
	#	Ne=NA,
	#	population_configurations=population_configurations,
	#	migration_matrix=migration_matrix,
	#	demographic_events=demographic_events)

	#dp.print_history()
#	num_replicates=1000
	simulation=msprime.simulate(
		population_configurations=population_configurations,
		#migration_matrix=migration_matrix,
		demographic_events=demographic_events,
		#num_replicates=num_replicates,
		length=Le, recombination_rate=Ru, mutation_rate=Mu
		#num_replicates=num_replicates, length=1e6, recombination_rate=2e-6, mutation_rate=2e-6
		)
	#T= np.zeros(num_replicates)
	for tree in simulation.trees():
		
		#tree=next(tree_sequence.trees())
		#print(tree.get_interval())
		#print(tree.get_time(tree.get_root())) 
		#tree.draw("test.svg", width=1000)
		#print(tree)
		LCA(tree,S1,S2)
	for variant in simulation.variants():
		t.write(str(variant.index) + "\t" +str(variant.position) + "\t" + str(variant.genotypes) + "\n")

		#print(tree.get_interval(), str(tree))
	#analytical = d / 2 + (d - 1) / (2 * M1)
	#print("Observed  =", np.mean(T))
	#print("Predicted =", analytical)
t = open("/media/pjay/DATAPART1/MigInf/test.txt","w")
sim1()
t.close()
