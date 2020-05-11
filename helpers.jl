using FileIO
using LasIO
using LinearAlgebra
using Printf

# first function to find circle from points
# must be [n * 2] matrix (n rows, 2 columns)
function findCircle1(A::Array{Float64})

    # form B
    B = [A[:,1].^2 + A[:,2].^2 A[:,1] A[:,2] ones(length(A[:,1]),1)]

    # evaluate singular value decomposition
    # F.U => orthogonal matrix
    # F.S => singular values
    # F.Vt => singular vector

    local F::SVD # need this to avoid UndefVarError
    try
        F = svd(B)
    catch
        println("SVD fails")
        return
    end

    # find index of minimum singular value
    # Smin = argmin(F.S)
    # get row of U that corresponds to the minimum S
    #u = F.Vt[Smin,:] # always 4
    u = F.Vt[4,:]

	# calculate radius
	m = (u[2]*u[2]+u[3]*u[3])/(4*u[1]*u[1]) - u[4]/u[1]
    r = sqrt(m)
    #@printf("Radius:%6.3f\n",r)
    # return the center (x, y) and radius
    -u[2]/(2*u[1]), -u[3]/(2*u[1]), r
end

# generate las file of perimeter points
# input is n X 2 matrix of eastings and northings
function lasCircle(dirpath::String, refString::String, h::LasHeader, eastCent::Float64, northCent::Float64, A::Array{Float64}, zf::Float64)

	i::UInt16 = 60000 # max is 0xFFFF = 65535
	b = Array{LasPoint,1}(undef,0)
	x::Int32 = 0
	y::Int32 = 0
	z::Int32 = round((zf - h.z_offset)/h.z_scale)
	class::UInt8 = 11 #

	(n,m) = size(A)

	for j = 1 : n
			x = round((A[j,1] - h.x_offset)/h.x_scale) # easting
			y = round((A[j,2] - h.y_offset)/h.y_scale) # northing
			t = LasPoint1(x,y,z,i,0,class,0,0,0,0)
	        push!(b,t)
	end

	# and finally the center point
	x = round((eastCent - h.x_offset)/h.x_scale)
	y = round((northCent - h.y_offset)/h.y_scale)
	t = LasPoint1(x,y,z,i,0,class,0,0,0,0)
	push!(b,t)

	update!(h,b)
	h.system_id = "circ"
	h.software_id = "PJM"

	save(string(dirpath,"\\",refString,".las"), h, b)
end

# generate las file of included points
function lasClump(dirpath::String, refString::String, h::LasHeader, eastCent::Float64, northCent::Float64, p::Array{Point}, zf::Float64)

	i::UInt16 = 60000 # max is 0xFFFF = 65535
	b = Array{LasPoint,1}(undef,0)
	x::Int32 = 0
	y::Int32 = 0
	z::Int32 = round((zf - h.z_offset)/h.z_scale)
	class::UInt8 = 12 #
	for pp in p
		if pp.d
			x = round((pp.x - h.x_offset)/h.x_scale)
			y = round((pp.y - h.y_offset)/h.y_scale)
			t = LasPoint1(x,y,z,i,0,class,0,0,0,0)
	        push!(b,t)
		end
	end

	# and finally the center point
	x = round((eastCent - h.x_offset)/h.x_scale)
	y = round((northCent - h.y_offset)/h.y_scale)
	t = LasPoint1(x,y,z,i,0,class,0,0,0,0)
	push!(b,t)

	update!(h,b)
	h.system_id = "clump"
	h.software_id = "PJM"

	save(string(dirpath,"\\",refString,".las"), h, b)
end
