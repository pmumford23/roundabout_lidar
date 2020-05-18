# roundabout_lidar
# Software for horizontal accuracy in airborne LiDAR data
#
# Peter Mumford  2020
# ------------------------------------------------------------
# School of Aviation, University of New South Wales, Australia
# ------------------------------------------------------------
#


# threshold circ

using FileIO
using LasIO
using LinearAlgebra
using Printf
using Plots
using Statistics

include("point_struct.jl")
include("helpers.jl")

function thresh_circ(dirpath::String, refString::String, h::LasHeader, d::Array{Point}, eastCent::Float64, northCent::Float64, surveyRadius::Float64, aVz::Float64, threshold::UInt16)

	boundingRadius = surveyRadius + 3.5
	limitRadius = surveyRadius - 2.0
	if limitRadius < 1.0
		limitRadius = 1.0
	end
	# removed lindex from here


	# build array of points for circle fitting
	f = Array{Float64,1}(undef, 0)
	rVec = Array{Float32,1}(undef, 0)
	flop::Int = 0

	i::Int = 1
	# search forward
	while i < length(d)# && d[i].r > limitRadius
		# check scan line
		if d[i].n # new scan line
			flop = 0
		elseif flop == 0 && d[i].i > threshold
				push!(f,d[i].x)
				push!(f,d[i].y)
				push!(rVec,d[i].r)
				flop = 1
		end
		i = i + 1
	end

	flop = 0
	# serach backward
	i = length(d)
	while i > 0# && d[i].r > limitRadius
		# check scan line
		if d[i].n # new scan line
			flop = 0
		elseif flop == 0 && d[i].i > threshold
				push!(f,d[i].x)
				push!(f,d[i].y)
				push!(rVec,d[i].r)
				flop = 1
		end
		i = i - 1
	end

	if length(f) < 10
		println("not enough points to form circle")
		return 0.0, 0.0, 0.0, 0.0
	end

	# get the mean and standard dev of the
	# candidate circumference points
	aVr = mean(rVec)
	stdr = std(rVec)
	@printf("Radius Average:%.1f  STD:%.1f\n",aVr,stdr)

	# build index of outliers
	# update recvec2
	rVec2 = Array{Float32,1}(undef, 0)
	q = Array{Int,1}(undef, 0)
	c = 1
	for rr in rVec
		#println("Trim:",rr," offset:",abs(rr - aVr))
		if abs(rr - aVr) > (stdr * 1.6) || (rr - surveyRadius) > 2.5  #1.6
			push!(q,c)
			c = c+1
			push!(q,c)
			c = c+1
			#@printf("outlier:%.1f\n",rr)
		else
			c = c+2
		end
	end
	# remove outliers
	@printf("Removing %i outliers\n",length(q)/2)
	deleteat!(f,q)

	# check there are enough points left
	if length(f) < 8
		println("not enough points to form circle")
		return 0.0, 0.0, 0.0, 0.0
	end

	# f is a linear array of e1, n1, e2, n2, e3, n3 ....
	# reshape into 2 X n array (2 rows, n columns)
	#  | e1 e2 e3 ...|
	#  | n1 n2 n3 ...|
	# transpose into n X 2 array (n rows, 2 columns)
	#  | e1 n1 |
	#  | e2 n2 |
	#  | e3 n3 |
	#  ........
	# finally, convert into Array{Float64,2} type required by function
	G = convert(Array{Float64,2},transpose(reshape(f,(2,:))))
	#println(G)
	East, North, r = findCircle1(G)
	#println("findCircle 1 --------")
	@printf("Center East:%11.3f North:%13.3f   Radius:%4.1f  std:%.3f\n",East,North,r,stdr)
	@printf("delta  East:%11.3f North:%13.3f\n",eastCent-East,northCent-North)


	if isnan(r)
		println("failed")
	else
		println("Making LAS file")
		lasCircle(dirpath, refString,h,East,North,G,aVz)
	end

	plot(G[:,1]',G[:,2]', seriestype = :scatter,
				 color=RGB(0.2, 0.2, 0.2),
				 title = string(refString), legend = false, aspect_ratio = 1,
				 xaxis=("East"),yaxis=("North"))
	savefig(string(dirpath,"\\",refString,".png"))

	# work out standard deviation of radius of points
	# using calculated center point
	rVec2 = Array{Float32,1}(undef, 0)
	for i = 1:length(G[:,1])
		rr = sqrt((G[i,1]-East)^2 + (G[i,2]-North)^2)
		push!(rVec2,rr)
	end
	stdr = std(rVec2)

	return eastCent-East, northCent-North, r, stdr

end
