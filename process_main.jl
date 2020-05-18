# roundabout_lidar
# Software for horizontal accuracy in airborne LiDAR data
#
# Peter Mumford  2020
# ------------------------------------------------------------
# School of Aviation, University of New South Wales, Australia
# ------------------------------------------------------------
#

# process main

# Arguments are:
# dirpath:  this is the directory path for saving plots, LAS files and output text from the process.
# infile:  the name of a text file, with each line providing details of the LAS file, roundabout surveyed radius and center point.
# outfile  the name for the output file.

using FileIO
using LasIO
using Printf
using CSV

include("point_struct.jl")
include("lasCropSort.jl")
include("grad_clump.jl")
include("grad_circ.jl")
include("thresh_circ.jl")


function main(dirpath::String, infile::String, outfile::String)
# open files for process output
	io1 = open(string(dirpath,"\\grad_circ_",outfile),"w")
	io2 = open(string(dirpath,"\\grad_clump_",outfile),"w")
	io3 = open(string(dirpath,"\\thresh_circ_",outfile),"w")

	@printf(io1,"date,process,type,ref_name,winsize,surv_radius,calc_radius,delta_radius,std_radius,delta_E,delta_N\r\n")
	@printf(io2,"date,process,type,ref_name,winsize,surv_radius,delta_E,delta_N\r\n")
	@printf(io3,"date,process,type,ref_name,winsize,surv_radius,calc_radius,delta_radius,std_radius,delta_E,delta_N\r\n")

	D = CSV.read(infile)

	for row in eachrow(D)
		println("===== ",row.date,"    ",row.refname,"    ",row.radius," =====")
		println(row.filename)
		h, d, aVz, threshold = lasCropSort(dirpath, row.filename, string(row.date,"_",row.refname), row.easting, row.northing, row.radius, row.winsize)

		println("Number of points:",length(d))
		println("")
		println("   --- gradiant circ ----")
		East, North, radius, std_rad = grad_circ(dirpath, string(row.date,"_",row.refname,"_grad_circ"), h, d, row.easting, row.northing, row.radius, aVz)
		dr = row.radius - radius
		@printf(io1,"%i,1,grad_circu,%s,%i,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f\r\n",row.date,row.refname,row.winsize,row.radius,radius,dr,std_rad,East, North)
		println("")
		println("   --- gradiant clump ----")
		East, North, radius, std_rad = grad_clump(dirpath, string(row.date,"_",row.refname,"_grad_clump"), h, d, row.easting, row.northing, row.radius, aVz)
		dr = 0.0
		@printf(io2,"%i,2,grad_clump,%s,%i,%.2f,%.3f,%.3f\r\n",row.date,row.refname,row.winsize,row.radius,East, North)

		if threshold != 0
			println("")
			println("   --- threshold circ ----")
			East, North, radius, std_rad = thresh_circ(dirpath, string(row.date,"_",row.refname,"_thresh_circ"), h, d, row.easting, row.northing, row.radius, aVz, threshold)
			if radius == 0.0
				dr = 0.0
			else
				dr = row.radius - radius
			end
			@printf(io3,"%i,3,thresh_circu,%s,%i,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f\r\n",row.date,row.refname,row.winsize,row.radius,radius,dr,std_rad,East, North)
			println("")
		else
			@printf(io3,"%i,3,thresh_circu,%s,%i,%.2f,,,,,\r\n",row.date,row.refname,row.winsize,row.radius)
		end
		println("Complete")
		println("")
		println("")
	end
	close(io1)
	close(io2)
	close(io3)

	println("Ends")
end
