#!/usr/bin/julia
### sisi ###
### sosos ###

using Combinatorics
using ArgParse
using ArgParse

function parse_commandline()
	s=ArgParseSettings()
	@add_arg_table s begin
		"--window","-w"
			help="size of the window"
			arg_type=Int
		"--increment","-s"
			help="size of the increment"
			arg_type=Int
			default=0
		"--file","-i"
			help="input file"
			required=true
		"--output","-o"
			help="output file"
			required=true

	end
	return parse_args(s)
end
parse_commandline()
S = 0 # segregating site

f = open("test.txt")
Fst=0
for l in eachline(f)
	S += 1
	if startswith(l,"#")
		a=split(l,r"\t|\n")
		global pop1=parse(Int,a[2])
		global pop2=parse(Int,a[3])
		global tot=pop1+pop2
		deb=pop1+1
		global lsize=parse(Float64,a[4])
		global compt=Array{Any}(binomial(tot,2),2)## Create an array of the size of the combination of n sample in group of 2
		compt[1:end,1]=0
		corres=Array{Any}(tot)
		corres[1:pop1]="A"
		corres[deb:tot]="B"
		corresAB=collect(combinations(corres,2))
		for i in 1:size(compt,1)
		 	if corresAB[i][1]==corresAB[i][2]
				compt[i,2]=corresAB[i][2]
			else
				compt[i,2]="O"
			end
		end

	else
		compta=Array{Any}(binomial(tot,2))## Create an array of the size of the combination of n sample in group of 2
		a=split(l, r"\t|\n|$")## split string to array
		t=include_string(a[3])## convert a string of the form [a,b] to an array
		#println(t)
		b=collect(combinations(t,2)) ############# Faire une simplification : si les pop1 premiere valeurs valle 0 ou valent 1, on rempli direct le tableau pour toute les combination AA, Sinon, on continue, et on increment un variable qui comple le nombre de site variabme intra espece. meme chose pour B. Verifier aussi tajima's D !
		difA=0 	
		difB=0
		difO=0
 		for i in 1:length(b)
 			if b[i][1]!=b[i][2]
				#println("pata")
 				compt[i,1] +=1
				compta[i] = 1
			else
				compta[i] = 0
 			end
		end
		o=0
		a=0
		b=0
		for i in 1:length(compta)
			if compt[i,2] == "O"
				difO += compta[i]
				o +=1
			end
			if compt[i,2] == "A"
				difA += compta[i]
				a +=1
			end
			if compt[i,2] == "B"
				difB += compta[i]
				b += 1
			end
		end
		PiA=difA/a
		PiB=difB/b
		PiO=difO/o
		Fst += (((PiO - ((PiA+PiB)/2))/PiO))/lsize#?????
	end
end

PiT=0
PiA=0
PiB=0
PiT=0
PiO=0
Pi0=0
Dxy=0
o=0
t=0
a=0
b=0


#for i in 1:size(compt,1)
#	PiT += 2*(1/tot)*(1/tot)*(compt[i,1]/lsize)
#end

for i in 1:size(compt,1)
	if compt[i,2] == "O"
		Dxy += compt[i,1]/lsize
#		PiT += compt[i,1]/lsize
#		t += 1	
		o +=1
	end
	if compt[i,2] == "A"
	#	PiA += compt[i,1]/lsize
		PiA += compt[i,1]
#		PiT += compt[i,1]/lsize
		a +=1
#		t += 1	
	end
	if compt[i,2] == "B"
		PiB += compt[i,1]
#		PiT += compt[i,1]/lsize
		b +=1
#		t += 1	
	end
end

PiAf=(PiA/lsize)/a
PiBf=(PiB/lsize)/b
#PiT=PiT/t
Dxy=Dxy/o
Da=Dxy-(PiAf+PiBf)/2

##### Tajima's D #################
a1=0
for i in 1:(pop1-1)
	a1+=1/i
end


a2=0
for i in 1:(pop1-1)
	a2+=1/i^2
end

b1=(pop1+1)/3*(pop1-1)
b2=(2*(pop1^2+pop1+3))/(9*pop1*(pop1-1))
c1=b1-1/a1
c2=b2-((pop1+2)/(a1*pop1))+a2/a1^2
e1=c1/a1
e2=c2/(a1^2+a^2)

M=S/a1 ##### S doit etre les siote variant dans pop 1, pas tout les site : a corriger ###
d=(PiA/pop1)-M
D=d/(sqrt(e1*S+e2*S*(S-1)))

################################

#println("PiT=",PiT)
println("PiA=",PiAf)
println("D=",D)
println("PiB=",PiBf)
#println("PiO=",PiO)
println("Fst=",Fst)
println("Dxy=",Dxy)
println("Da=",Da)
close(f)
