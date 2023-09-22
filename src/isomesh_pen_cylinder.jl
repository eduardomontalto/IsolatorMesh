#=
3D Tetrahedral Mesh Generator using DistMesh

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 3d mesh for a cylinder section using pentahedral elements based on DistMesh. The size of the elements can be nonuniform and the location of the nodes are arbitrary in general, except for predefined points.
=#

#----------------------------INPUT----------------------------------------
#=
R:      radius of cylinder
h:      height of block
nz:     number of elements in z direction
er:     target element aspect ratio at center
x0:     zero x coordinate of the block
y0:     zero y coordinate of the block
z0:     zero z coordinate of the block
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y,z) coordinates. Size: nn x 3
t:      array of connectivity. Size: nel x 6
=#

#----------------------------START----------------------------------------

function pen_cylinder(R::Float64,
    h::Float64,
    nz::Int,
    er::Int;
    x0::Float64=0., 
    y0::Float64=0.,
    z0::Float64=0.,
    half_x::Bool=false,
    half_y::Bool=false)::Tuple{Matrix{Float64},Matrix{Int}}
    
    function fd(p)
        d1 = dcircle(p,x0,y0,R)
        if half_x & !half_y
            d2 = drectangle(p,x0,x0+1.5*R,y0-1.5*R,y0+1.5*R)
            d = dintersect(d1,d2)
        elseif half_y & !half_x
            d2 = drectangle(p,x0-1.5*R,x0+1.5*R,y0,y0+1.5*R)
            d = dintersect(d1,d2)
        elseif half_x & half_y
            d2 = drectangle(p,x0,x0+1.5*R,y0,y0+1.5*R)
            d = dintersect(d1,d2)
        else
            d = d1
        end
        return d
    end

    if er == 1
        h0 = h/nz
    else
        h0 = h/nz*(1+er)/2
    end

    function fh(p)
        r = sqrt.((p[:,1] .- x0).^2 + (p[:,2] .- y0).^2)
        eh = er .- (er-1)*r/R
        return eh
    end

    if half_x & !half_y
        pfix = [x0 y0;x0 y0-R;x0 y0+R]
    elseif half_y & !half_x
        pfix = [x0 y0;x0-R y0;x0+R y0]
    elseif half_x & half_y
        pfix = [x0 y0;x0 y0+R;x0+R y0]
    else
        pfix = [x0 y0]
    end
    bbox = [x0-1.5*R y0-1.5*R;x0+1.5*R y0+1.5*R]
    pcs, tcs = distmesh2d(fd, fh, h0, bbox, pfix=pfix)
    tcs = clean_delaunay!(pcs,tcs)
    tcs = check_direction!(pcs,tcs)
    ncs = size(pcs,1)
    p = [pcs z0*ones(ncs)]
    t = zeros(0,6)

    for i in 1:nz
        p = [p;[pcs (z0+i*h/nz)*ones(ncs)]]
        t = [t;[tcs.+(i-1)*ncs tcs.+i*ncs]]
    end

    return p, t
end