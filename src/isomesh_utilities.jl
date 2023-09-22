#=
2D Mesh Generator Utilities

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= Utility functions to clean the mesh and plotting.
=#

#----------------------------START----------------------------------------

# function delaunay(p)
#     tri = pyimport("matplotlib.tri")
#     t = tri[:Triangulation](p[:,1], p[:,2])
#     return Int64.(t[:triangles] .+ 1)
# end

# Function to identify edges in 2D
function all_edges(t)
    etag = vcat(t[:,[1,2]], t[:,[2,3]], t[:,[3,1]])
    etag = hcat(sort(etag, dims=2), 1:3*size(t,1))
    etag = sortslices(etag, dims=1)
    dup = all(etag[2:end,1:2] - etag[1:end-1,1:2] .== 0, dims=2)[:]
    keep = .![false;dup]
    edges = etag[keep,1:2]
    emap = cumsum(keep)
    invpermute!(emap, etag[:,3])
    emap = reshape(emap,:,3)
    dup = [dup;false]
    dup = dup[keep]
    bndix = findall(.!dup)
    return edges, bndix, emap
end

# Function to identify boundary nodes in 2D
function boundary_nodes(t)
    edges, bndix, _ = all_edges(t)
    return unique(edges[bndix,:][:])
end;

# Function to remove spurious elements
function clean_delaunay!(p::Matrix{Float64},
    t::Matrix{Int};
    ϵ::Float64 = 1e-4)::Matrix{Int}
    dim = size(p,2)
    nel = size(t,1)
    if dim == 2
        A = zeros(nel)
        for i in 1:nel
            A[i] = tri_area(p[t[i,:],1], p[t[i,:],2])       
        end
        id = A .> ϵ
        t = [t[id,1] t[id,2] t[id,3]]
    elseif dim == 3
        V = zeros(nel)
        for i in 1:nel
            V[i] = tet_volume(p[t[i,:],1], p[t[i,:],2],p[t[i,:],3])
        end
        id = V .> ϵ
        t = [t[id,1] t[id,2] t[id,3] t[id,4]]
    end
    return t
end

# Function to calculate area of a triangle
function tri_area(x::Vector{Float64}, 
    y::Vector{Float64})::Float64
    p1 = [x[1]-x[2],y[1]-y[2],0]
    p2 = [x[1]-x[3],y[1]-y[3],0]
    An = cross(p1,p2)
    A = 1/2*norm(An,2)
    return A
end

# Function to calculate volume of a tetrahedron
function tet_volume(x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64})::Float64
    p1 = [x[1]-x[2],y[1]-y[2],z[1]-z[2]]
    p2 = [x[1]-x[3],y[1]-y[3],z[1]-z[3]]
    p3 = [x[1]-x[4],y[1]-y[4],z[1]-z[4]]
    V = abs(dot(cross(p1,p2),p3))/6
    return V
end

# Function to plot mesh in 2D
function mplot2d(p::Matrix{Float64},
    t::Matrix{Int};
    idtype::Bool=false)
    
    Cl = [[0 0.447 0.741],[0.85 0.325 0.098],[0.929 0.694 0.125],
    [0.494 0.184 0.556],[0.466 0.674 0.188],[0.301 0.745 0.933],
    [0.635 0.078 0.184],[0 0 1],[0 0.5 0],[1 0 0],[0 0.75 0.75],
    [0.75 0 0.75],[0.75 0.75 0],[0.25 0.25 0.25],[0.75 0.75 0.75]]
    nne = size(t,2)
    if nne == 2
        bars = unique(sort(t,dims=2),dims=1)
    elseif nne == 3
        bars = [t[:,[1,2]];t[:,[2,3]];t[:,[1,3]]]
        bars = unique(sort(bars,dims=2),dims=1)
    elseif nne == 4
        bars = [t[:,[1,2]];t[:,[2,3]];t[:,[3,4]];t[:,[1,4]]]
        bars = unique(sort(bars,dims=2),dims=1)
    end
    if idtype && nne == 2
        for i in axes(bars,1)
            plot(p[bars[i,:],1],p[bars[i,:],2],"-o",color=Cl[nne],linewidth=1.0,markersize=2.0)
        end
    elseif idtype
        for i in axes(bars,1)
            plot(p[bars[i,:],1],p[bars[i,:],2],"-",color=Cl[nne],linewidth=0.5)
        end
    elseif nne == 2
        for i in axes(bars,1)
            plot(p[bars[i,:],1],p[bars[i,:],2],"-k",linewidth=1.5)
        end
    else
        for i in axes(bars,1)
            plot(p[bars[i,:],1],p[bars[i,:],2],"-k",linewidth=0.5)
        end
    end
    axis("equal")
end

#Function to plot mesh in 3D
function mplot3d(p::Matrix{Float64},
    t::Matrix{Int};
    fdp::Function=x->ones(size(x,1)),
    idtype::Bool=false)

    xrange = (maximum(p[:,1]) - minimum(p[:,1]))/2
    yrange = (maximum(p[:,2]) - minimum(p[:,2]))/2
    zrange = (maximum(p[:,3]) - minimum(p[:,3]))/2
    range = max(xrange,yrange,zrange)
    xmid = (maximum(p[:,1]) + minimum(p[:,1]))/2
    ymid = (maximum(p[:,2]) + minimum(p[:,2]))/2
    zmid = (maximum(p[:,3]) + minimum(p[:,3]))/2
    nne = size(t,2)
    if nne == 4
        n = cross(p[t[1,2],:]-p[t[1,1],:],p[t[1,3],:]-p[t[1,1],:])
        if isapprox(dot(n,p[t[1,4],:]-p[t[1,1],:]),0.0)
            ftype = "quad"
        else
            ftype = "tri"
        end
    end
    if nne == 3
        faces = copy(t)
    elseif nne == 4 && ftype == "tri"
        faces = [t[:,[1,2,3]];t[:,[1,2,4]];t[:,[1,3,4]];t[:,[2,3,4]]]
    elseif nne == 4 && ftype == "quad"
        faces = [t[:,[1,2,3]];t[:,[3,4,1]]]
    elseif nne == 6
        faces = [t[:,[1,2,5,4]];t[:,[2,3,6,5]];t[:,[3,1,4,6]]]
        faces = [t[:,1:3];t[:,4:6];faces[:,[1,2,3]];faces[:,[3,4,1]]]
    elseif nne == 8
        faces = [t[:,1:4];t[:,5:8];t[:,[1,2,6,5]];t[:,[2,3,7,6]];t[:,[3,4,8,7]];t[:,[4,1,5,8]]]
        faces = [faces[:,[1,2,3]];faces[:,[3,4,1]]]
    end
    faces = unique(sort(faces,dims=2),dims=1)
    idfdp = (fdp(p[faces[:,1],:]).>-0.01) .& (fdp(p[faces[:,2],:]).>-0.01) .& (fdp(p[faces[:,3],:]).>-0.01)
    faces = faces[idfdp,:]
    if idtype
        if (nne == 3) || (nne == 4 && ftype == "quad")
            plot_trisurf(p[:,1],p[:,2],p[:,3],triangles=faces.-1,color="lightsalmon",alpha=1.0,antialiased=true,shade=false,edgecolors="k",linewidth=0.5)
        else
            plot_trisurf(p[:,1],p[:,2],p[:,3],triangles=faces.-1,color="lightskyblue",alpha=1.0,antialiased=true,shade=false,edgecolors="k",linewidth=0.25)
        end
    else
        plot_trisurf(p[:,1],p[:,2],p[:,3],triangles=faces.-1,color="lightgray",alpha=1.0,antialiased=true,shade=false,edgecolors="k",linewidth=0.25)
    end
    gca().set_xlim(xmid-range,xmid+range)
    gca().set_ylim(ymid-range,ymid+range)
    gca().set_zlim(zmid-range,zmid+range)
end