# common stuff

import Base: isless
# overload function for sort
# time sort


mutable struct Point
    t::Float64 # GPS time
    x::Float64 # easting
    y::Float64 # northing
    z::Float64 # height
    i::UInt16 # lidar point intensity
	r::Float64 # radius from surveyed center (calculated)
    s::Int32 # scan line number
	n::Bool # start of new scan line
	g::Float64 # change of intensity along line at point
	u::Bool # use or not
	d::Bool # inside point
	e::Bool # edge point
end

isless(a::Point,b::Point) = isless(a.t,b.t)
