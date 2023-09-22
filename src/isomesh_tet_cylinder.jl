#=
3D Tetrahedral Mesh Generator using DistMesh

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 3d mesh for a cylinder section using tetrahedral elements based on DistMesh. The size of the elements can be nonuniform and the location of the nodes are arbitrary in general, except for predefined points.
=#

#----------------------------INPUT----------------------------------------
#=
REQUIRED
R:      radius of cylinder
h:      height of block
nz:     vector with number of elements in z direction on interior and boundary
        such that nz[1] < nz[2]

OPTIONAL
x0:     zero x coordinate of the cylinder
y0:     zero y coordinate of the cylinder
z0:     zero z coordinate of the cylinder
rb:     width of boundary sections in r
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y,z) coordinates. Size: nn x 3
t:      array of connectivity. Size: nel x 4
=#

#----------------------------START----------------------------------------

function tet_cylinder(R::Float64,
    h::Float64,
    nz::Vector{Int};
    x0::Float64=0., 
    y0::Float64=0.,
    z0::Float64=0.,
    rb::Float64=0.2*R,
    half_x::Bool=false,
    half_y::Bool=false)::Tuple{Matrix{Float64},Matrix{Int}}
    
    function fd(p)
        d1 = dcylinder(p,x0,y0,z0,R,h)
        if half_x & !half_y
            d2 = dblock(p,x0,x0+1.5*R,y0-1.5*R,y0+1.5*R,z0-0.5*h,z0+1.5*h)
            d = dintersect(d1,d2)
        elseif half_y & !half_x
            d2 = dblock(p,x0-1.5*R,x0+1.5*R,y0,y0+1.5*R,z0-0.5*h,z0+1.5*h)
            d = dintersect(d1,d2)
        elseif half_x & half_y
            d2 = dblock(p,x0,x0+1.5*R,y0,y0+1.5*R,z0-0.5*h,z0+1.5*h)
            d = dintersect(d1,d2)
        else
            d = d1
        end
        return d
    end
    h0 = h/nz[1]
    hr = nz[1]/nz[2]
    function fh(p)
        eh = ones(size(p,1))
        r = sqrt.((p[:,1] .- x0).^2 + (p[:,2] .- y0).^2)
        if nz[1] != nz[2]
            for i in axes(p,1)
                if r[i] > R - rb
                    eh[i] = 1 + (hr-1)/rb*(r[i]-(R-rb))
                end
            end
        end
        return eh
    end
    
    if half_x & !half_y
        pfix = [x0 y0 z0+h;x0 y0-R z0;x0 y0+R z0;x0 y0-R z0+h;x0 y0+R z0+h]
    elseif half_y & !half_x
        pfix = [x0 y0 z0+h;x0-R y0 z0;x0+R y0 z0;x0-R y0 z0+h;x0+R y0 z0+h]
    elseif half_x & half_y
        pfix = [x0 y0 z0;x0 y0+R z0;x0+R y0 z0;x0 y0 z0+h;x0 y0+R z0+h;x0+R y0 z0+h]
    else
        pfix = [x0 y0 z0+h]
    end
    bbox = [x0-1.5*R y0-1.5*R z0-0.5*h;x0+1.5*R y0+1.5*R z0+1.5*h]
    p, t = distmesh3d(fd, fh, h0, bbox, pfix=pfix)
    t = clean_delaunay!(p,t)

    return p, t
end