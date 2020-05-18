# roundabout_lidar
# Software for horizontal accuracy in airborne LiDAR data
#
# Peter Mumford  2020
# ------------------------------------------------------------
# School of Aviation, University of New South Wales, Australia
# ------------------------------------------------------------
#

# gradiant circumference

using FileIO
using LasIO
using LinearAlgebra
using Printf
#using Plots
using Statistics

include("point_struct.jl")
include("helpers.jl")

function grad_circ(dirpath::String, refString::String, h::LasHeader, d::Array{Point}, eastCent::Float64, northCent::Float64, surveyRadius::Float64, aVz::Float64)

	# work out expected contrast
	# surveyRadius is the provides the extent to work with
	# but if the surveyed center and LiDAR scan are offset
	# by more than about 1.5m there will be issues with
	# estimating the intensities.
	# get average intensity for points near the center
	innerIntensity::Float64 = 0.0
	c::Int32 = 0
	for dd in d
		# get the center 2/3 of roundabout inner
		# in some roundabouts the center is darker
		if dd.r < (surveyRadius * 0.6)
			c = c+1
			innerIntensity = innerIntensity + convert(Float64,dd.i)
		end
	end
	innerIntensity = innerIntensity/c

	outerIntensity::Float64 = 0.0
	c = 0
	for dd in d
		if dd.r > (surveyRadius + 1.5)
			c = c+1
			outerIntensity = outerIntensity + convert(Float64,dd.i)
		end
	end
	outerIntensity = outerIntensity/c

	# estimate threshold if thresh parrameter not set
	# values between 2400 and 5000 are sensible
	deltaI_thresh::Float64 = innerIntensity-outerIntensity

	if deltaI_thresh > 10000.0
		deltaI_thresh = deltaI_thresh/2.5
	elseif deltaI_thresh > 6000.0
		deltaI_thresh = deltaI_thresh/2.0
	elseif deltaI_thresh > 3000.0
		deltaI_thresh = deltaI_thresh - deltaI_thresh/10.0
	end

	@printf("Estimated threshold:%.0f\n",deltaI_thresh)

	# removed lindex from here

	# build array of points for circle fitting
	f = Array{Float64,1}(undef, 0)
	rVec = Array{Float32,1}(undef, 0)
	gMax::Float64 = 0
	gMin::Float64 = 0
	gMax_i::Int = 1
	gMin_i::Int = 1
	# find largest negative and positive gradients
	# in each scan line
	# cannot use argmax or argmin due to the line segments
	for i = 1:length(d)
		# check scan line
		if d[i].n # new scan line

			# check if the points are sensible
			if gMax > deltaI_thresh# && d[gMax_i].r >  surveyRadius - 1.5
				push!(f,(d[gMax_i].x+d[gMax_i+1].x)/2.0)
				push!(f,(d[gMax_i].y+d[gMax_i+1].y)/2.0)
				push!(rVec,(d[gMax_i].r+d[gMax_i+1].r)/2.0)

			end
			if gMin < -deltaI_thresh# && d[gMin_i].r > surveyRadius - 1.5
				push!(f,(d[gMin_i].x+d[gMin_i+1].x)/2.0)
				push!(f,(d[gMin_i].y+d[gMin_i+1].y)/2.0)
				push!(rVec,(d[gMin_i].r+d[gMin_i+1].r)/2.0)

			end

			gMax = 0
			gMin = 0
		else
			if d[i].g > gMax
				gMax = d[i].g
				gMax_i = i
			end
			if d[i].g < gMin
				gMin = d[i].g
				gMin_i = i
			end
		end
	end

	if length(f) < 10
		println("not enough points to form circle")
		return 0,0,0,0,0
	end

	# get the mean and standard dev of the
	# candidate circumference points
	aVr = mean(rVec)
	stdr = std(rVec)
	@printf("Radius Average:%.1f  STD:%.1f\n",aVr,stdr)

	# build index of outliers
	q = Array{Int,1}(undef, 0)
	c = 1
	for rr in rVec
		#println("Trim:",rr," offset:",abs(rr - aVr))
		if abs(rr - aVr) > (stdr * 2.0) || abs(rr - surveyRadius) > 2.0  #1.6
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

	@printf("Center East:%11.3f North:%13.3f   Radius:%4.1f  std:%.3f\n",East,North,r,stdr)
	@printf("delta  East:%11.3f North:%13.3f\n",eastCent-East,northCent-North)

	if isnan(r)
		println("Circle fails")
	else
		lasCircle(dirpath,refString,h,East,North,G,aVz)
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
