#=
DistMesh Mesh Generator

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This file includes all the necessary functions to create an unstructured triangular mesh in 2D or an unstructured tetrahedral mesh in 3D using DistMesh
=#

#----------------------------INPUT----------------------------------------
#=
fd:     distance function to define geometry of the body 
fh:     function to indicate relative size of elements
h0:     desired size of elements
bbox:   box to generate intial random points for Delaunay triangulation
        [x1 y1 z1;x2 y2 z2]
dptol:  tolerance side length change to stop interpolation
ttol:   tolerance side length change to retriangulate
Fscale: parameter >=1 to ensure that triangulation fills the geometry
deltat: fictitious time step for time integration of the change in length
geps:   tolerance for defining if element is inside the boundary
eps:    machine tolerance
deps:   delta for finite difference computation of normal vector to boundary
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y) coordinates. Size: nn x 2
t:      array of connectivity. Size: nel x 3
=#

#----------------------------START----------------------------------------

function distmesh2d(fd::Function, 
    fh::Function, 
    h0::Float64, 
    bbox::Matrix{Float64}; 
    pfix::Matrix{Float64}=zeros(0,2),
    dptol::Float64=0.001, 
    ttol::Float64=0.1, 
    Fscale::Float64=1.2, 
    deltat::Float64=0.2, 
    geps::Float64=0.001*h0, 
    eps::Float64=1e-16, 
    deps::Float64=sqrt(eps)*h0)::Tuple{Matrix{Float64},Matrix{Int}}

    # Create initial distribution in bounding box (equilateral triangles)
    x = collect(bbox[1,1]:h0:bbox[2,1])
    y = collect(bbox[1,2]:h0*sqrt(3)/2:bbox[2,2])
    x = x*ones(length(y))'
    y = y*ones(size(x,1))'
    x[2:2:end,:] .+= h0/2
    p = [x[:] y[:]]

    # Remove points outside the region, apply the rejection method
    p = p[isless.(fd(p),geps),:]
    r0 = 1.0./fh(p).^2
    p = [pfix; p[isless.(rand(size(p,1)),r0./maximum(r0)),:]]
    N = size(p,1)

    pold = Inf*ones(size(p))
    dp = 1
    t = AbstractArray{Int,2}
    bars = AbstractArray{Int,2}
    nc = 0

    while dp > dptol
        # Retriangulation by the Delaunay algorithm
        if maximum(sqrt.(sum((p-pold).^2,dims=2))/h0) > ttol
            pold = p
            mesh = delaunay(p, [:Qt,:Qbb,:Qc])
            t = mesh.simplices
            pmid = (p[t[:,1],:] + p[t[:,2],:] + p[t[:,3],:])/3
            t = t[isless.(fd(pmid),-geps),:]
            # Describe each bar by a unique pair of nodes
            bars = [t[:,[1,2]];t[:,[1,3]];t[:,[2,3]]]
            bars = unique(sort(bars,dims=2),dims=1)
        end

        # Move meshpoints based on bar lengths L and forces F
        barvec = p[bars[:,1],:] - p[bars[:,2],:]
        L = sqrt.(sum(barvec.^2,dims=2))
        hbars = fh((p[bars[:,1],:]+p[bars[:,2],:])/2)
        L0 = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2))
        F = max.(L0-L,0)
        Fvec = F./L.*barvec
        Ftot = sparse(bars[:,1],ones(size(F,1))*1,Fvec[:,1],N,2)
        Ftot += sparse(bars[:,1],ones(size(F,1))*2,Fvec[:,2],N,2)
        Ftot += sparse(bars[:,2],ones(size(F,1))*1,-Fvec[:,1],N,2)
        Ftot += sparse(bars[:,2],ones(size(F,1))*2,-Fvec[:,2],N,2)
        Ftot = Matrix(Ftot)
        Ftot[1:size(pfix,1),:] .= 0
        p += deltat*Ftot

        # Bring outside points back to boundary
        d = fd(p)
        ix = isless.(0,d)
        dgradx = (fd(p[ix,:] .+ deps*[1 0])-d[ix])/deps
        dgrady = (fd(p[ix,:] .+ deps*[0 1])-d[ix])/deps
        p[ix,:] = p[ix,:]-d[ix].*[dgradx dgrady]

        # Termination criteria
        dp = maximum(sqrt.(sum(deltat*Ftot[isless.(d,-geps),:].^2,dims=2))/h0)
        nc += 1
        if nc == 10000
            break
        end
    end

    return p, t
end

function distmesh3d(fd::Function, 
    fh::Function, 
    h0::Float64,
    bbox::Matrix{Float64};
    pfix::Matrix{Float64}=zeros(0,3), 
    dptol::Float64=0.001, 
    ttol::Float64=0.1, 
    Fscale::Float64=1.1, 
    deltat::Float64=0.2, 
    geps::Float64=0.1*h0, 
    eps::Float64=1e-16, 
    deps::Float64=sqrt(eps)*h0, 
    ncmax::Int=250)::Tuple{Matrix{Float64},Matrix{Int}}

    # Create initial distribution in bounding box (equilateral triangles)
    x = collect(bbox[1,1]:h0:bbox[2,1])
    y = collect(bbox[1,2]:h0:bbox[2,2])
    z = collect(bbox[1,3]:h0:bbox[2,3])
    pp = Vector{Array}(undef,3)
    pp[1] = [i for i in x, j in y, k in z]
    pp[2] = [j for i in x, j in y, k in z]
    pp[3] = [k for i in x, j in y, k in z]
    p = [vec(pp[1]) vec(pp[2]) vec(pp[3])]

    # Remove points outside the region, apply the rejection method
    p = p[isless.(fd(p),geps),:]
    r0 = fh(p)
    p = [pfix; p[isless.(rand(size(p,1)),minimum(r0)^3 ./(r0.^3)),:]]
    N = size(p,1)

    pold = Inf*ones(size(p))
    dp = 1
    t = AbstractArray{Int,2}
    bars = AbstractArray{Int,2}
    nc = 0

    while dp > dptol
        # Retriangulation by the Delaunay algorithm
        if maximum(sqrt.(sum((p-pold).^2,dims=2))/h0) > ttol
            pold = copy(p)
            mesh = delaunay(p, [:Qt,:Qbb,:Qc])
            t = mesh.simplices
            pmid = (p[t[:,1],:] + p[t[:,2],:] + p[t[:,3],:] + p[t[:,4],:])/4
            t = t[isless.(fd(pmid),-geps),:]
            # Describe each bar by a unique pair of nodes
            bars = [t[:,[1,2]];t[:,[1,3]];t[:,[1,4]];t[:,[2,3]];t[:,[2,4]];t[:,[3,4]]]
            bars = unique(sort(bars,dims=2),dims=1)
            nc += 1
        end

        # Limit number of retriangulations
        if nc == ncmax; break; end

        # Move meshpoints based on bar lengths L and forces F
        barvec = p[bars[:,1],:] - p[bars[:,2],:]
        L = sqrt.(sum(barvec.^2,dims=2))
        hbars = fh((p[bars[:,1],:]+p[bars[:,2],:])/2)
        L0 = hbars*Fscale*(sum(L.^3)/sum(hbars.^3))^(1/3)
        F = max.(L0-L,0)
        Fvec = F./L.*barvec
        Ftot = sparse(bars[:,1],ones(Int,size(F,1))*1,zeros(size(F,1)),N,3)
        for i in 1:3
            Ftot += sparse(bars[:,1],ones(Int,size(F,1))*i,Fvec[:,i],N,3)
            Ftot += sparse(bars[:,2],ones(Int,size(F,1))*i,-Fvec[:,i],N,3)
        end
        Ftot = Matrix(Ftot)
        Ftot[1:size(pfix,1),:] .= 0
        p += deltat*Ftot

        # Bring outside points back to boundary
        d = fd(p)
        ix = isless.(0,d)
        dgradx = (fd(p[ix,:] .+ deps*[1 0 0])-d[ix])/deps
        dgrady = (fd(p[ix,:] .+ deps*[0 1 0])-d[ix])/deps
        dgradz = (fd(p[ix,:] .+ deps*[0 0 1])-d[ix])/deps
        p[ix,:] = p[ix,:]-d[ix].*[dgradx dgrady dgradz]

        # Termination criteria
        dp = maximum(sqrt.(sum(deltat*Ftot[isless.(d,-geps),:].^2,dims=2))/h0)
    end

    return p, t
end

# Distance functions

function drectangle(p,x1,x2,y1,y2)
    d = -min.(-y1.+p[:,2],y2.-p[:,2],-x1.+p[:,1],x2.-p[:,1])
    return d 
end

function dcircle(p,xc,yc,r)
    d = @. sqrt((p[:,1]-xc)^2 + (p[:,2]-yc).^2) - r
    return d
end

function dunion(d1,d2)
    d = min.(d1,d2)
    return d
end

function ddiff(d1,d2)
    d = max.(d1,-d2)
    return d
end

function dintersect(d1,d2)
    d = max.(d1,d2)
    return d
end

function dcylinder(p,xc,yc,zc,R,h)
    d1 = dcircle(p[:,1:2], xc, yc, R)
    d2 = @. p[:,3] - zc - h
    d3 = @. zc - p[:,3]
    d4 = sqrt.(d1.^2 + d2.^2)
    d5 = sqrt.(d1.^2 + d3.^2)
    d = max.(d1,d2,d3)
    d[(d1.>0) .& (d2.>0)] = d4[(d1.>0) .& (d2.>0)]
    d[(d1.>0) .& (d3.>0)] = d5[(d1.>0) .& (d3.>0)]
    return d
end

function dblock(p,x1,x2,y1,y2,z1,z2)
    d = -min.(-z1.+p[:,3],z2.-p[:,3],-y1.+p[:,2],y2.-p[:,2],-x1.+p[:,1],x2.-p[:,1])
    return d
end

function huniform(p)
    h = ones(size(p,1))
    return h
end