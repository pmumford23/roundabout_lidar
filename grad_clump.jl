# grad_clump

using FileIO
using LasIO
using LinearAlgebra
using Printf
using Plots
using Statistics

include("point_struct.jl")
include("helpers.jl")

function grad_clump(dirpath::String, refString::String, h::LasHeader, d::Array{Point}, eastCent::Float64, northCent::Float64, surveyRadius::Float64, aVz::Float64)

	c::Int32 = 0

	# work out expected contrast
	# surveyRadius is the provides the extent to work with
	# but if the surveyed center and LiDAR scan are offset
	# by more than about 1.5m there will be issues with
	# estimating the intensities.
	# get average intensity for points near the center
	innerIntensity::Float64 = 0.0

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

	@printf("Intensity: Outer:%.0f  Inner:%.0f  diff:%.0f\n",outerIntensity,innerIntensity,innerIntensity-outerIntensity)

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


	# find the maximum and minimum intensity gradient
	# along each scan line. Make index of these points
	f = Array{Int,1}(undef, 0)
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
			if gMax > deltaI_thresh && gMin < -deltaI_thresh
				# set d[i].d to indicate the point is inside
				for k = (gMax_i + 1):gMin_i
					d[k].d = true
				end
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

	# find average easting and northing
	# from all points found to be inside
	# the roundabout
	East::Float64 = 0.0
	North::Float64 = 0.0
	c = 0
	for dd in d
		if dd.d
			East = East + dd.x
			North = North + dd.y
			c = c + 1
		end
	end

	East = East / c
	North = North / c

	@printf("Center East:%11.3f North:%13.3f\n",East,North)
	@printf("delta  East:%11.3f North:%13.3f\n",eastCent-East,northCent-North)

	if isnan(East)
		println("failed")
	else
		lasClump(dirpath,refString,h,East,North,d,aVz)
	end

	eastCent-East, northCent-North, 0.0, 0.0
end
