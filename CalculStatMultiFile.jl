#!/usr/bin/julia
using Combinatorics

function pairwisedifference(n,s)
   egal=binomial(s,2) + binomial((n-s),2)
   tot=binomial(n,2)
   val=(tot-egal)/tot
    return val
end

function tajimasD(pop,Pi,var) #pop = sample size, Pi= nucleotide diversity, var=Number of segregating site in pop
	##### Tajima's D pop1#################
	a1=0
	for i in 1:(pop-1)
	a1+=1/i
	end
	
	
	a2=0
	for i in 1:(pop-1)
	a2+=1/i^2
	end
	
	b1=(pop+1)/3*(pop-1)
	b2=(2*(pop^2+pop+3))/(9*pop*(pop-1))
	c1=b1-1/a1
	c2=b2-((pop+2)/(a1*pop))+a2/(a1^2)
	e1=c1/a1
	e2=c2/(a1^2+a2)
	
	M=var/a1
	d=Pi-M
	D=d/(sqrt(e1*var+e2*var*(var-1)))
	return D
end
	
dictFile=Dict() #dict that will contain output file

cd(ARGS[1]) # grep the directiory contain split sequences
dirfiles=readdir()
for file in dirfiles # for each sequence interval
	println(file)
	if countlines(file) == 1 # if the file just contain header (empty), all value equal 0.0
		f = open("$file")
		head=readlines(f)[1]
		a=split(head,r"\t|\n")
		win=a[5]	
		if haskey(dictFile,win)
			continue
		else
			dictFile[win]=open("Interval.$win.txt","w")
			write(dictFile[win],"PiA","\t","PiB","\t","Tajima'D_pop1","\t","Tajima'D_pop2","\t","Fst","\t","Dxy","\t","Da","\n")
		end
		PiAf=0.0
		PiBf=0.0
		Dxy=0.0
		Fst=0.0
		Da=0.0
		D1=0.0
		D2=0.0
	else 
		f = open("$file")
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
				win=a[5]	# time interval corresponding to this sequence interval
				global sizeinter=a[7]-a[8] # size of the sequence interval
				if haskey(dictFile,win) # if this time interval ever have its outfile, do nothing
					continue
				else
					dictFile[win]=open("Interval.$win.txt","w") #else, open it in the directionnary value
					write(dictFile[win],"PiA","\t","PiB","\t","Tajima'D_pop1","\t","Tajima'D_pop2","\t","Fst","\t","Dxy","\t","Da","\n") # print header
				end
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
				p1=0 # Number of individual with a mutation in pop1
				for i in 1:pop1
					p1 += t[i]
				end
			
				p2=0 # Number of individual with a mutation in pop2
				for i in pop1+1:tot
					p2 += t[i]
				end
			
				if p1 == 0 || p1==pop1
					PiAf=0 ## if no or all individual have the mutation, piA = O
				else
					varA+=1 ## IF not, is a site with variation
					Pi=pairwisedifference(pop1,p1) # it 's a trick to calculate pi juste with the number individual with the mutation
					PiA+=Pi # We increment Pi for the whole sequence
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
				
				PiT=((pop1*pop2)-(p2*p1+(pop2-p2)*(pop1-p1)))/(pop1*pop2) # Tricks to calculate piO/Dxy between population 
				PiO+=PiT
				PiOf=PiT
				Fst += (((PiOf - ((PiAf+PiBf)/2))/PiOf)) 
			end
		end
	#	PiAf=(PiA/lsize) # We divide the sum of pi by the number of site
	#	PiBf=(PiB/lsize)
	#	Dxy=(PiO/lsize)
	#	Fst=(Fst/lsize)
	#	Da=Dxy-(PiAf+PiBf)/2
	############
	##
	## 	Check the good way to calculate mean of the statistique on the diferent sequence interval (to take into account that all interval have different size . A mean of means is not the same ,thing that an overall mean) Same thing for standard deviation. 
	##
	##############
		PiAf=PiA/sizeinter #divide by the number of site in the interval
		PiBf=PiB/sizeinter
		Dxy=PiO/sizeinter
		Fst=Fst/sizeinter
		Da=Dxy-(PiAf+PiBf)/2

	if PiA != 0.0 # to avoid division by 0.0
			D1=tajimasD(pop1,PiA,varA)
		else
			D1=0.0
		end
		if PiB != 0.0
			D2=tajimasD(pop2,PiB,varB)
		else
			D2=0.0
		end
		close(f)
	end
###### Tajima's D pop1#################
#a1=0
#for i in 1:(pop1-1)
#a1+=1/i
#end
#
#
#a2=0
#for i in 1:(pop1-1)
#a2+=1/i^2
#end
#
#b1=(pop1+1)/3*(pop1-1)
#b2=(2*(pop1^2+pop1+3))/(9*pop1*(pop1-1))
#c1=b1-1/a1
#c2=b2-((pop1+2)/(a1*pop1))+a2/(a1^2)
#e1=c1/a1
#e2=c2/(a1^2+a2)
#
#M=varA/a1
#d=PiA-M
#D1=d/(sqrt(e1*varA+e2*varA*(varA-1)))
#
#################################
###### Tajima's D pop2#################
#a1=0
#for i in 1:(pop2-1)
#a1+=1/i
#end
#
#
#a2=0
#for i in 1:(pop2-1)
#a2+=1/i^2
#end
#
#b1=(pop2+1)/3*(pop2-1)
#b2=(2*(pop2^2+pop2+3))/(9*pop2*(pop2-1))
#c1=b1-1/a1
#c2=b2-((pop2+2)/(a1*pop2))+a2/(a1^2)
#e1=c1/a1
#e2=c2/(a1^2+a2)
#
#M=varB/a1
#d=PiB-M
#D2=d/(sqrt(e1*varB+e2*varB*(varB-1)))
#
#################################
	#write(dictFile[win],"PiA","\t","PiB","\t","Tajima'D pop1","\t","Tajima'D pop2","\t","Fst","\t","Dxy","\t","Da","\n")
	write(dictFile[win],string(PiAf ,"\t",PiBf ,"\t",D1 ,"\t",D2 ,"\t",Fst ,"\t",Dxy ,"\t",Da,"\n"))

	println("PiA=",PiAf)
	println("PiB=",PiBf)
	println("Tajima'D pop1 =",D1)
	println("Tajima'D pop2 =",D2)
	println("Fst=",Fst)
	println("Dxy=",Dxy)
	println("Da=",Da)
end
for i in values(dictFile) # close all open file
	close(i)
end

dirfiles=readdir()
for file in dirfiles
	if startswith(file, "Interval.") # if the file are the one containing the statistique
		a=readdlm("$file")
		b=a[2:end,:] #avaid the header line
		c=mapslices(mean,b,1) # mean of column
		println(c)
	end
end
cd("..")
rm("window", recursive=true)
