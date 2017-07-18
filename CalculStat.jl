#!/usr/bin/julia
tic()
using Combinatorics

function pairwisedifference(n,s)
   egal=binomial(s,2) + binomial((n-s),2)
   tot=binomial(n,2)
   val=(tot-egal)/tot
    return val
end

f = open(ARGS[1])

S = 0 # segregating site Number
Fst=0 # Fst value
PiA=0 # Fst value
PiB=0 # Fst value
PiO=0 # Fst value
varA=0 # Number of site where all A sample have different genotype
varB=0

for l in eachline(f) #Reading msprime file line per line
	if startswith(l,"#") # Reading simulation input value
		a=split(l,r"\t|\n")
		global pop1=parse(Int,a[2]) # Number of sample in pop 1
		global pop2=parse(Int,a[3]) # Number of sample in pop 2
		global tot=pop1+pop2 # total number of sample
		deb=pop1+1 #the index of the pop2 sample start a pop1 index +1
		global lsize=parse(Float64,a[4]) # length of the sequence simulated
		global compt=Array{Any}(binomial(tot,2),2) ## Create an array of the size of the combination of n sample in group of 2
		compt[1:end,1]=0  # first dimension of array contains pairwise difference between sample
		corres=Array{Any}(tot) # a vector for ecah sample which contains the population affiliation of the sample
		corres[1:pop1]="A" 
		corres[deb:tot]="B"
		corresAB=collect(combinations(corres,2)) # create an array of the possible combination of the sample --> same size of the "compt" array, in the same order
		for i in 1:size(compt,1)
			if corresAB[i][1]==corresAB[i][2] # if the two letter stored in the array are the same, the sample belong to the same pop
				compt[i,2]=corresAB[i][2] # in this case, the second dimension of the array have the name of the popilation
			else
				compt[i,2]="O" # if not, it the combination of individual from the two pop, and it is named "O"
			end
		end
	
	else # reading the simulation result
		S += 1 # each line is a segrating site
		a=split(l, r"\t|\n|$") ## split string to array
		t=include_string(a[3])## convert a string of the form [a,b] to an array
		p1=0
		for i in 1:pop1
			p1 += t[i]
		end
	
		p2=0
		for i in pop1+1:tot
			p2 += t[i]
		end
	
		if p1 == 0 || p1==pop1
			PiAf=0 
		else
			varA+=1
			Pi=pairwisedifference(pop1,p1)
			PiA+=Pi
			PiAf=Pi
		end
	
		if p2 == 0 || p2==pop2
			PiBf=0 
		else
			varB+=1
			Pi=pairwisedifference(pop2,p2)
			PiB+=Pi
			PiBf=Pi
		end
		
		PiT=((pop1*pop2)-(p2*p1+(pop2-p2)*(pop1-p1)))/(pop1*pop2)
		PiO+=PiT
		PiOf=PiT
		Fst += (((PiOf - ((PiAf+PiBf)/2))/PiOf)) #?????
	end
end
PiAf=(PiA/lsize)
PiBf=(PiB/lsize)
Dxy=(PiO/lsize)
Fst=(Fst/lsize)
println(varA)
Da=Dxy-(PiA+PiB)/2

##### Tajima's D #################
a1=0
for i in 1:(pop1-1)
a1+=1/i
end
println(PiA)


a2=0
for i in 1:(pop1-1)
a2+=1/i^2
end

b1=(pop1+1)/3*(pop1-1)
b2=(2*(pop1^2+pop1+3))/(9*pop1*(pop1-1))
c1=b1-1/a1
c2=b2-((pop1+2)/(a1*pop1))+a2/(a1^2)
e1=c1/a1
e2=c2/(a1^2+a2)

M=varA/a1 ##### S doit etre les siote variant dans pop 1, pas tout les site : a corriger ###
#d=(PiA/pop1)-M
d=PiA-M
println(d)
D=d/(sqrt(e1*varA+e2*varA*(varA-1)))

################################

println("PiA=",PiAf)
println("D=",D)
println("PiB=",PiBf)
println("Fst=",Fst)
println("Dxy=",Dxy)
println("Da=",Da)
close(f)
toc()
