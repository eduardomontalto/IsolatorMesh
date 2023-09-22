#=
3D Tetrahedral Mesh Generator using DistMesh

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 3d mesh for a block section using tetrahedral elements based on DistMesh. The size of the elements can be nonuniform and the location of the nodes are arbitrary in general, except for predefined points.
=#

#----------------------------INPUT----------------------------------------
#=
REQUIRED
d:          depth of block
w:          width of block
h:          height of block
nz:         vector with number of elements in z on interior and boundary
            such that nz[1] < nz[2]

OPTIONAL
x0:         zero x coordinate of the block
y0:         zero y coordinate of the block
z0:         zero z coordinate of the block
xb:         vector with widths of boundary sections in x
yb:         vector with widths of boundary sections in y
half_x:     model only half of the bearing in x (defaults to false)
half_y:     model only half of the bearing in y (defaults to false)
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y,z) coordinates. Size: nn x 3
t:      array of connectivity. Size: nel x 4
=#

#----------------------------START----------------------------------------

function tet_block(d::Float64,
    w::Float64,
    h::Float64,
    nz::Vector{Int};
    x0::Float64=0., 
    y0::Float64=0.,
    z0::Float64=0.,
    xb::Float64=0.2*d,
    yb::Float64=0.2*w,
    half_x::Bool=false,
    half_y::Bool=false)::Tuple{Matrix{Float64},Matrix{Int}}

    function fd(p)
        if half_x & !half_y
            d1 = dblock(p,x0+d/2,x0+d,y0,y0+w,z0,z0+h)
        elseif half_y & !half_x
            d1 = dblock(p,x0,x0+d,y0+w/2,y0+w,z0,z0+h)
        elseif half_x & half_y
            d1 = dblock(p,x0+d/2,x0+d,y0+w/2,y0+w,z0,z0+h)
        else
            d1 = dblock(p,x0,x0+d,y0,y0+w,z0,z0+h)
        end
        return d1
    end
    h0 = h/nz[1]
    hr = nz[1]/nz[2]
    xb = [xb*!half_x,xb]
    yb = [yb*!half_y,yb]
    function fh(p)
        eh = ones(size(p,1))
        if nz[1] != nz[2]
            for i in axes(p,1)
                if (p[i,1] < x0 + xb[1]) & (y0 + yb[1] <= p[i,2] <= y0 + w - yb[2])
                    eh[i] = hr + (1-hr)/xb[1]*(p[i,1]-x0)
                elseif (p[i,1] > x0 + d - xb[2]) & (y0 + yb[1] <= p[i,2] <= y0 + w - yb[2])
                    eh[i] = 1 + (hr-1)/xb[2]*(p[i,1]-(x0+d-xb[2]))
                elseif (p[i,2] < y0 + yb[1]) & (x0 + xb[1] <= p[i,1] <= x0 + d - xb[2])
                    eh[i] = hr + (1-hr)/yb[1]*(p[i,2]-y0)
                elseif (p[i,2] > y0 + w - yb[2]) & (x0 + xb[1] <= p[i,1] <= x0 + d - xb[2])
                    eh[i] = 1 + (hr-1)/yb[2]*(p[i,2]-(y0+w-yb[2]))
                elseif (p[i,1] < x0 + xb[1]) & (p[i,2] < y0 + yb[1])
                    eh[i] = hr*(x0+xb[1]-p[i,1])/xb[1]*(y0+yb[1]-p[i,2])/yb[1]
                    eh[i] += hr*(p[i,1]-x0)/xb[1]*(y0+yb[1]-p[i,2])/yb[1]
                    eh[i] += 1.0*(p[i,1]-x0)/xb[1]*(p[i,2]-y0)/yb[1]
                    eh[i] += hr*(x0+xb[1]-p[i,1])/xb[1]*(p[i,2]-y0)/yb[1]
                elseif (p[i,1] < x0 + xb[1]) & (y0 + w - yb[2] < p[i,2])
                    eh[i] = hr*(x0+xb[1]-p[i,1])/xb[1]*(y0+w-p[i,2])/yb[2]
                    eh[i] += 1.0*(p[i,1]-x0)/xb[1]*(y0+w-p[i,2])/yb[2]
                    eh[i] += hr*(p[i,1]-x0)/xb[1]*(p[i,2]-y0-w+yb[2])/yb[2]
                    eh[i] += hr*(x0+xb[1]-p[i,1])/xb[1]*(p[i,2]-y0-w+yb[2])/yb[2]
                elseif (p[i,1] > x0 + d - xb[2]) & (p[i,2] < y0 + yb[1])
                    eh[i] = hr*(x0+d-p[i,1])/xb[2]*(y0+yb[1]-p[i,2])/yb[1]
                    eh[i] += hr*(p[i,1]-x0-d+xb[2])/xb[2]*(y0+yb[1]-p[i,2])/yb[1]
                    eh[i] += hr*(p[i,1]-x0-d+xb[2])/xb[2]*(p[i,2]-y0)/yb[1]
                    eh[i] += 1.0*(x0+d-p[i,1])/xb[2]*(p[i,2]-y0)/yb[1]
                elseif (p[i,1] > x0 + d - xb[2]) & (y0 + w - yb[2] < p[i,2])
                    eh[i] = 1.0*(x0+d-p[i,1])/xb[2]*(y0+w-p[i,2])/yb[2]
                    eh[i] += hr*(p[i,1]-x0-d+xb[2])/xb[2]*(y0+w-p[i,2])/yb[2]
                    eh[i] += hr*(p[i,1]-x0-d+xb[2])/xb[2]*(p[i,2]-y0-w+yb[2])/yb[2]
                    eh[i] += hr*(x0+d-p[i,1])/xb[2]*(p[i,2]-y0-w+yb[2])/yb[2]
                end
            end
        end
        return eh
    end

    if half_x & !half_y
        bbox = [x0+d/2 y0 z0;x0+d y0+w z0+h]
        pfix = [x0+d/2 y0+w/2 z0+h;x0+d/2 y0 z0;x0+d/2 y0 z0+h;x0+d/2 y0+w z0;x0+d/2 y0+w z0+h;x0+d y0 z0;x0+d y0 z0+h;x0+d y0+w z0;x0+d y0+w z0+h]
    elseif half_y & !half_x
        bbox = [x0 y0+w/2 z0;x0+d y0+w z0+h]
        pfix = [x0+d/2 y0+w/2 z0+h;x0 y0+w/2 z0;x0 y0+w/2 z0+h;x0 y0+w z0;x0 y0+w z0+h;x0+d y0+w/2 z0;x0+d y0+w/2 z0+h;x0+d y0+w z0;x0+d y0+w z0+h]
    elseif half_x & half_y
        bbox = [x0+d/2 y0+w/2 z0;x0+d y0+w z0+h]
        pfix = [x0+d/2 y0+w/2 z0;x0+d/2 y0+w/2 z0+h;x0+d/2 y0+w z0;x0+d/2 y0+w z0+h;x0+d y0+w/2 z0;x0+d y0+w/2 z0+h;x0+d y0+w z0;x0+d y0+w z0+h]
    else
        bbox = [x0 y0 z0;x0+d y0+w z0+h]
        pfix = [x0 y0 z0;x0 y0 z0+h;x0 y0+w z0;x0 y0+w z0+h;x0+d y0 z0;x0+d y0 z0+h;x0+d y0+w z0;x0+d y0+w z0+h]
    end
    p, t = distmesh3d(fd, fh, h0, bbox, pfix=pfix)
    t = clean_delaunay!(p,t)

    return p, t
end