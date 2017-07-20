#!/usr/bin/julia

f = open(ARGS[1])
w = parse(Int,ARGS[2]) ## Number of interval

MaxAr = zeros(1) ## Array containing the divergence time
Window = zeros(1) 
for l in eachline(f)
	a=split(l, r"\t|\n|$") ## split string to array
	b=parse(Float64,a[2]) # grep de divergence time
	push!(MaxAr,b) # push the divergence time into the divergence array
end

max=maximum(MaxAr) # give the maximum divergence time
WinSize=max/w # determine the size of an time interval : Maximum divergence time divide by the number of interval ## possible amelioration different interval size depending on the divergence (short in recent time, large in ancient time)
close(f)

f = open(ARGS[1])
out = open(ARGS[3],"w") #File to write the window pos
#PosEnd="0.0" # the start of the window is 0, the begin of the sequence

for l in eachline(f) # just for the first line, thing are different
	write(out,"0.0","-") #Write the first position
	a=split(l, r"\(|, |\)|\t|\n|$") # split string into array . pos 2 contain the first position of the interval, pos 3 containt the end of the interval, and pos 5 contain the divergence time
	win=div(parse(Float64,a[5]),WinSize) # contain the time  interval of this sequence interval. It correspond to division (entiere) of the divergence time by the winSize.
	push!(Window,win ) # put the first interval in the array
	break
end



for l in eachline(f) # each line contain an interval of position in the sequence and the corresponding divergence time. Start at line 2, as the fist line was ever readed in the last loop
	global a=split(l, r"\(|, |\)|\t|\n|$") # split string into array . pos 2 contain the first position of the interval, pos 3 containt the end of the interval, and pos 5 contain the divergence time
	win=div(parse(Float64,a[5]),WinSize) # contain the time  interval of this sequence interval. It correspond to division (entiere) of the divergence time by the winSize.
	if win != Window[end] ## the divergence time of the window is stocked in an array. The last value of this array correspond to the divergence time of the previous line. If $win is different to the last value, that mean that in a the beginning of this sequence interval, we move from a time interval to an other one.
		write(out,a[2],"\t", string(Window[end]), "\n", a[2],"-") #Write 
		push!(Window,win ) # put the new interval in the array
	end
end		
write(out,a[3],"\t",string(Window[end])) # Write the last position and the last time interval
close(f)
close(out)


#dictOut=Dict() #Contain the reference file for each time interval 
win = open(ARGS[3])
lines = readlines(win)
seq = open(ARGS[4])

for ln in eachline(seq) # open just the first line of the seq file
	global c=ln
	break
end

#for i in 1:w
#	dictOut[i] = "win$i.txt"	
#	write(dictOut[i], a)	
#end


count=1 ## start a the first line of the window file
a=split(lines[count], r"\t|-|$")  
w=open("Inter$(a[1])-$(a[2]).txt", "w") # open the first file containing the variant for this window

mkdir("window")
cd("window") ## create and moove to a directory containt all the window

for l in eachline(seq)
	s=split(l, r"\t|\n|$") ## split string to array
	#println(s[2])
	if parse(Float64,s[2]) < parse(Float64,a[2])
		write(w,l)
	else
		println(s[2],"\t", a[2])
		while parse(Float64,s[2]) > parse(Float64,a[2])
			count+=1
			global a=split(lines[count], r"\t|-|$")
			global w=open("Inter$(a[1])-$(a[2]).txt","w")
		end
		println(s[2],"\t", a[2], "\n")
		write(w,l)
	end
end
cd("..")
#println(a[1])
#w=write("Inter$(a[1])-$(a[2]).txt")

#
#write("win5.txt","soso")
#
##for l in eachline(seq)
#	
#
#
#println(a[1])
#mkdir("window")
#cd("window")
#cd("..")
#rm("window", recursive=true)
#
