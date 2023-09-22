#=
2D Tri Mesh Generator using DistMesh

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 2d mesh for a rectangular area using tri elements based on DistMesh. The size of the elements can be nonuniform and the location of the nodes are arbitrary in general, except for predefined points.
=#

#----------------------------INPUT----------------------------------------
#=
w:      width of rectangular area
h:      height of rectangular area
ny:     vector with number of elements in y direction on left and right sides
x0:     zero x coordinate of the block
y0:     zero y coordinate of the block
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y) coordinates. Size: nn x 2
t:      array of connectivity. Size: nel x 4
=#

#----------------------------START----------------------------------------

function tri(w::Float64,
    h::Float64,
    ny::Vector{Int};
    x0::Float64=0.,
    y0::Float64=0.)::Tuple{Matrix{Float64},Matrix{Int}}

    fd(p) = drectangle(p,x0,x0+w,y0,y0+h)
    fh(p) = ny[2]/ny[1]*(w .+ x0 .- p[:,1])/w + (p[:,1] .- x0)/w
    hmax = h/minimum(ny)
    hmin = h/maximum(ny)
    if ny[1] == ny[2]
        h0 = (hmax+hmin)/2
    else
        h0 = hmax/√2
    end
    bbox = [x0-0.5*w y0-0.5*h;x0+1.5*w y0+1.5*h]
    pfix = [x0 y0;x0+w y0;x0 y0+h;x0+w y0+h]
    for i=1:2
        if isless.(ny,reverse(ny))[i] || isequal.(ny,reverse(ny))[i]
            pfix = [pfix;
            [(x0+w*iseven(i))*ones(ny[i]-1) h/ny[i]*collect(1:ny[i]-1)]]
        end
        if  isless.(reverse(ny),ny)[i] || isequal.(ny,reverse(ny))[i]
            pfix = [pfix;
            [pfix[i,1]+(isodd(i)-iseven(i))*hmin/√2,pfix[i,2]+hmin/√2]']
            pfix = [pfix;
            [pfix[i+2,1]+(isodd(i)-iseven(i))*hmin/√2,pfix[i+2,2]-hmin/√2]']
        end
    end
    dptol = 0.0001
    Fscale = 1.2
    p, t = distmesh2d(fd, fh, h0, bbox, pfix=pfix, dptol=dptol, Fscale=Fscale)

    # Delete nodes on the sides with fixed nodes
    pkeep = ones(Bool,size(p,1))
    for i in axes(p,1)
        if (isapprox(p[i,1],x0) && ny[1] <= ny[2]) || (isapprox(p[i,1],x0+w) && ny[2] <= ny[1])
            if p[i,2] ∉ pfix[pfix[:,1] .== p[i,1],2]
                pkeep[i] = false
            end
        end
    end
    p = p[pkeep,:]
    t = delaunay(p)
    t = clean_delaunay!(p,t)

    return p, t
end