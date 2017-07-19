#!/usr/bin/julia

f = open(ARGS[1])
w = parse(Int,ARGS[2])

MaxAr = zeros(1) 
Window = zeros(1) 
for l in eachline(f)
	a=split(l, r"\t|\n|$") ## split string to array
	b=parse(Float64,a[2])
	push!(MaxAr,b)
end

max=maximum(MaxAr)
WinSize=max/w
close(f)

f = open(ARGS[1])
out = open(ARGS[3],"w")
PosEnd="0"
println(typeof(PosEnd))
for l in eachline(f)
	a=split(l, r"\(|, |\)|\t|\n|$")
	#println(a[5])
	win=div(parse(Float64,a[5]),WinSize)
	println(typeof(PosEnd))
#	println(win)
	if win != Window[end]
		push!(Window,win )	
#		print(PosEnd, "\n", a[2],"-")
		write(out,PosEnd,"\n", a[2],"-")
	else
		PosEnd=a[3]
	end
end		
write(out,PosEnd)
close(f)
close(out)
#println(Window)
