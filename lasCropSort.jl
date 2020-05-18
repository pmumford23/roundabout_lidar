# roundabout_lidar
# Software for horizontal accuracy in airborne LiDAR data
#
# Peter Mumford  2020
# ------------------------------------------------------------
# School of Aviation, University of New South Wales, Australia
# ------------------------------------------------------------
#

# crop

using FileIO
using LasIO

include("point_struct.jl")

function lasCropSort(dirpath::String, fileName::String, refString::String, eastCent::Float64, northCent::Float64, surveyRadius::Float64, winsize::Int)

    println("open las file")
    h, p = load(fileName,mmap=true)

	boundingRadius = surveyRadius + 3.5

	# scale coordinates to LAS file values
	# scaledPoint::Int32 = round((eastCent - offset - h.x_offset)/h.x_scale)

	c::Int32 = 0
	r::Float64 = 0
	d = Array{Point,1}(undef, 0)
	x::Float64 = 0
	y::Float64 = 0
	z::Float64 = 0
	aVz::Float64 = 0

	for pp in p
	# scale LAS point to float
		# reconstitute values
		x = pp.x * h.x_scale + h.x_offset
		y = pp.y * h.y_scale + h.y_offset
        z = pp.z * h.z_scale + h.z_offset

	# calculate distance between point and center
		r = sqrt((x-eastCent)^2 + (y-northCent)^2)

		if (r < boundingRadius)
			c = c+1
			aVz = aVz + z
			newPoint = Point(pp.gps_time,x,y,z,pp.intensity,r,0,false,0.0,false,false,false)
			push!(d,newPoint)
		end
	end
	aVz = aVz / c

	@printf("Number of points:%i  Average z:%.1f\n",length(d),aVz)

	# remove high points
	k = Array{Int,1}(undef, 0)
	c = 0
	for dd in d
		c = c+1
		if dd.z > aVz + 0.5
			push!(k,c)
		end
	end
	println("Removing ",length(k)," high points")
	deleteat!(d,k)

	# sort into scan lines
	# sort by time
    sort!(d)

	dt::Float64 = 0 # delta time
	t1::Float64 = 0
	t2::Float64 = 0
	dtMax::Float64 = 0
	dtMin::Float64 = 10
	dtThresh::Float64 = 0

    # get min and max delta time between points
    for dd in d

		t1 = t2
		t2 = dd.t
		dt = t2 - t1
		#println("Delta:",dt)
		if t1 > 0 # don't use first entry
			if dt > dtMax
				dtMax = dt
			end
			if dt > 0 && dt < dtMin # filter multi-returns
				dtMin = dt
			end
		end
	end

	# threshold is set halfway between Min and Max
	dtThresh = (dtMax + dtMin) / 2.0

	#println("delta Max:",dtMax," Min:",dtMin," Thresh:",dtThresh)

	s::Int32 = 1 # scan line number

	# sort into flight lines ----
	# create index to line
	lindex = Array{Int,1}(undef,0)

	#lc = 0 # line counter
	c = 0 # index counter
	pcMin = 0
	t1 = t2 = 0
	for dd in d
		c = c+1
		t1 = t2
		t2 = dd.t
		dt = t2 - t1
		if t1 > 0
			if dt > dtThresh
				s = s+1
				dd.n = true # mark new scan line
				push!(lindex,c-1) # next index to end of line
			end
		end
		dd.s = s
		#lc = lc+1
	end
	push!(lindex,c)

	println("Number of scan lines:",s)

	# lindex is an array with the indexes to the last point in a scan line
	# mark the internal points in each scan line for use in calculating intensity gradients
	# do not use the bounding (w) points at the start and end of scan line
	#
	# window w=3      [          ][          ]
	#                  1   2   3   4   5   6   7   8   9
	# index                    ^
	# virtual point              ^    (average position of points 3 and 4)
	#
	# set window size 1, 2, 3, 4, 5, 6, ...
	if winsize > 0
		w = winsize
	else
		w = 4
	end

	c = 0
	for i in lindex
		if (i - c) > (w*2 + 4) # don't use short lines
			#for j = (c+4):(i-3)
			for j = (c+w+1):(i-w)
				d[j].u = true
			end
		end
		c = i
	end

	zb::Float64 = 0.0
	zf::Float64 = 0.0
	# process flight lines
	# form gradient over multiple points
	for i = 1:length(d)
		if d[i].u == true

			for j = 1:w
				zb = zb + convert(Float64,d[i-j+1].i)
				zf = zf + convert(Float64,d[i+j].i)
			end

			zb = zb/w
			zf = zf/w
			#diff = zf - zb
			#@printf(io,"zb: %7.2f  zf: %7.2f  diff:%7.2f\r\n",zb,zf,diff)
			d[i].g = zf - zb # gradient in delta intensity
			zb = 0
			zf = 0
		else
			d[i].g = 0
		end
		#println("G:",d[i].g)
	end

	# histagram of intensities
	# normalised by point radius (gives two nice peaks)
	# range is 10,000 to 40,000
	numberOfBins::Int32 = 30 #20
	binwidth::Float64 = 30000/numberOfBins
	hist = zeros(Float64,numberOfBins)
	for dd in d
		for bin = 1:numberOfBins
			if dd.i >= (bin*binwidth + 10000) && dd.i < ((bin+1)*binwidth + 10000)
				hist[bin] = hist[bin] + boundingRadius/dd.r
			end
		end
	end

	bin::Int = 2
	hist_d1 = zeros(Float64,numberOfBins)
	# traverse the hist curve
	# create smoothed gradient
	while bin < numberOfBins - 1
		hist_d1[bin] = (hist[bin+1] + hist[bin+2] - hist[bin] - hist[bin-1])/4.0
		bin = bin + 1
	end

	# find the first bin where the d1 goes from negative to positive
	# only search the middle section
	bin = numberOfBins / 5
	neg::Bool = false
	threshold::UInt16 = 0
	while bin < (numberOfBins - numberOfBins / 5)
			if hist_d1[bin] < -10 # make sure gradient has gone negative
				neg = true
			end
			# find first zero crossing
			if neg && hist_d1[bin] < 0  && (hist_d1[bin+1] > 0 || hist_d1[bin+2] > 0) && hist[bin+1] < 50
				#if (hist_d1[bin+2] + hist_d1[bin+3] - hist_d1[bin] - hist_d1[bin+1]) > 6
					bin = bin + 1
					threshold = convert(UInt16,round((bin+1)*binwidth + 10000))
					break
				#end
			end

		bin = bin + 1
	end

	if threshold == 0
		println("Threshold not detected !!!")
		plot(hist, w = 3, linecolor = :black, line = :solid, title = refString, legend = false)
		plot!(hist_d1, w = 3, linecolor = :black, line = :dash)
		savefig(string(dirpath,"\\",refString,".png"))
	else
		plot(hist, w = 3, linecolor = :black, line = :solid, title = refString, legend = false)
		plot!(hist_d1, w = 3, linecolor = :black, line = :dash)
		vline!([bin], w = 1, linecolor = :black, line = :solid)
		savefig(string(dirpath,"\\",refString,".png"))
	end

	println("Bin:",bin,"  Threshold:",threshold)

	# return las header, point array, average Z
	return h, d, aVz, threshold
end
