#=
3D Mesh Generator for Isolators in MSC Marc

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 3d mesh for a rectangular or circular isolator to input in MSC Marc. Function uses an unstructured mesh with tetrahedral elements created with DistMesh.
=#

#----------------------------INPUT----------------------------------------
#=
REQUIRED
a:          depth of rectangular cross section
b:          width of rectangular cross section
R:          radius of circular cross section
Ri:         radius of internal hole for annular cross section
h:          total height of rubber layers
nl:         number of rubber layers
epl:        number of elements per rubber layer
cs_type:    cross section type: "circular", "rectangular" or "annular"     

OPTIONAL
half_x:     model only half of the bearing in x (defaults to false)
half_y:     model only half of the bearing in y (defaults to false)
xb:         depth of boundary sections with refined mesh (defaults to h/nl)
yb:         width of boundary sections with refined mesh (defaults to h/nl)
rb:         radial distance of boundary with refined mesh (defaults to h/nl)
hex_type:   hex element type (7 disp., 84 mixed) (defaults to 7)
quad_type:  quad element type (18 membrane, 147 rebar) (defaults to 18)
pen_type:   pen element type (136 disp.) (defaults to 136)
tet_type:   tet element type (134 disp., 157 mixed) (defaults to 157)
tri_type:   tri element type (158 membrane) (defaults to 158)
conID:      contact body ID of quad elements (defaults to 1)
mat_geoID:  geometry ID of rubber matrix elements (defaults to 1)
mat_matID:  material ID of rubber matrix elements (defaults to 1)
mat_ornID:  orientation ID of rubber matrix elements (defaults to 1)
rnf_geoID:  geometry ID of reinforcement elements (defaults to 2)
rnf_matID:  material ID of reinforcement elements (defaults to 2)
root:       root directory where output files will be saved (defaults to pwd())
filename:   name of output text files (defaults to nothing)
plotmesh:   plot generated mesh (defaults to false)
plotid:     identify element types in plot (defaults to true)
=#

#----------------------------OUTPUT---------------------------------------
#=
FILES
mesh.txt:   text file containing connectivity and coordinates
set.txt:    text file containing element sets
=#

#----------------------------START----------------------------------------

function mesh3d(;a::Float64=200.0,
    b::Float64=a,
    R::Float64=100.0,
    Ri::Float64=0.0,
    h::Float64=50.0,
    nl::Int=5,
    cs_type::String="circular", 
    epl::Int=3,
    eplb::Int=epl,
    er::Int=1,
    nne::Int=4,  
    half_x::Bool=false,
    half_y::Bool=false, 
    xb::Float64=h/nl,
    yb::Float64=h/nl,
    rb::Float64=h/nl,
    hex_type::Int=7,
    quad_type::Int=18,
    pen_type::Int=136,
    tet_type::Int=157, 
    tri_type::Int=158, 
    conID::Int=1, 
    mat_geoID::Int=1, 
    mat_matID::Int=1, 
    mat_ornID::Int=0, 
    rnf_geoID::Int=2, 
    rnf_matID::Int=2, 
    root=pwd(), 
    filename=nothing, 
    emesh_root=pwd(),
    emesh_filename=nothing,
    plotmesh::Bool=false, 
    plotid::Bool=true)

    # Create mesh
    if cs_type == "rectangular" && nne == 6
        nne = 8
    end
    if cs_type == "annular"
        nne = 8
    end
    hl = h/nl
    hr = hl/eplb
    if nne == 4
        if cs_type == "circular"
            p, tmat = tet_cylinder(R,hl,[epl,eplb],rb=rb,half_x=half_x,half_y=half_y)
        elseif cs_type == "rectangular"
            p, tmat = tet_block(a,b,hl,[epl,eplb],x0=-a/2,y0=-b/2,xb=xb,yb=yb,half_x=half_x,half_y=half_y)
        end
        p[isapprox.(p[:,3],0.0,atol=hr/5),3] .= 0.0
        p[isapprox.(p[:,3],hl,atol=hr/5),3] .= hl
        if half_x; p[isapprox.(p[:,1],0.0,atol=hr/5),1] .= 0.0; end
        if half_y; p[isapprox.(p[:,2],0.0,atol=hr/5),2] .= 0.0; end
        if cs_type == "rectangular"
            if !half_x; p[isapprox.(p[:,1],-a/2,atol=hr/5),1] .= -a/2; end
            p[isapprox.(p[:,1],a/2,atol=hr/5),1] .= a/2
            if !half_y; p[isapprox.(p[:,2],-b/2,atol=hr/5),2] .= -b/2; end
            p[isapprox.(p[:,2],b/2,atol=hr/5),2] .= b/2
        end
        pl1 = copy(p)
        tl1 = copy(tmat)
        for i in 2:nl
            if iseven(i)
                pl2 = pl1 .- [0. 0. hl]
                pl2 .*= [1. 1. -1.]
                pl2 .+= [0. 0. (i-1)*hl]
            else
                pl2 = pl1 .+ [0. 0. (i-1)*hl]
            end
            tl2 = copy(tl1)
            p, _, tl2 = merge_mesh!(p, pl2, tl2)
            tmat = [tmat;tl2]
        end
        tmat = check_direction(p,tmat)
    elseif nne == 6 && cs_type == "circular"
        p, tmat = pen_cylinder(R,h,epl*nl,er,half_x=half_x,half_y=half_y)
        if half_x; p[isapprox.(p[:,1],0.0,atol=hr/5),1] .= 0.0; end
        if half_y; p[isapprox.(p[:,2],0.0,atol=hr/5),2] .= 0.0; end
    elseif nne == 8 
        if cs_type == "circular" || cs_type == "annular"
            p, tmat = hex_cylinder(h,epl*nl,emesh_root,emesh_filename)
            if half_x; p[isapprox.(p[:,1],0.0,atol=hr/5),1] .= 0.0; end
            if half_y; p[isapprox.(p[:,2],0.0,atol=hr/5),2] .= 0.0; end
        elseif cs_type == "rectangular"
            p, tmat = hex_block(a,b,h,epl*nl,er,x0=-a/2,y0=-b/2,half_x=half_x,half_y=half_y)
        end
    end

    # Create connectivity of reinforcement elements
    ϵ = 1e-3
    if nne == 4
        trnf = zeros(Int,0,3)
        faces = [tmat[:,[1,2,3]];tmat[:,[1,2,4]];tmat[:,[1,3,4]];tmat[:,[2,3,4]]]
        faces = unique(sort(faces,dims=2),dims=1)
        for i in 1:nl-1
            idl = (isapprox.(p[faces[:,1],3],i*hl,atol=ϵ,rtol=√eps())) .&
                (isapprox.(p[faces[:,2],3],i*hl,atol=ϵ,rtol=√eps())) .&
                (isapprox.(p[faces[:,3],3],i*hl,atol=ϵ,rtol=√eps()))
            trnf = [trnf;faces[idl,:]]
        end
        trnf = check_direction(p,trnf)
    elseif nne == 6
        trnf = zeros(Int,0,3)
        faces = [tmat[:,1:3];tmat[:,4:6]]
        faces = unique(faces,dims=1)
        for i in 1:nl-1
            idl = (isapprox.(p[faces[:,1],3],i*hl,atol=ϵ,rtol=√eps())) .&
                (isapprox.(p[faces[:,2],3],i*hl,atol=ϵ,rtol=√eps())) .&
                (isapprox.(p[faces[:,3],3],i*hl,atol=ϵ,rtol=√eps()))
            trnf = [trnf;faces[idl,:]]
        end
    elseif nne == 8
        trnf = zeros(Int,0,4)
        faces = [tmat[:,1:4];tmat[:,5:8]]
        faces = unique(faces,dims=1)
        for i in 1:nl-1
            idl = (isapprox.(p[faces[:,1],3],i*hl,atol=ϵ,rtol=√eps())) .&
                (isapprox.(p[faces[:,2],3],i*hl,atol=ϵ,rtol=√eps())) .&
                (isapprox.(p[faces[:,3],3],i*hl,atol=ϵ,rtol=√eps())) .&
                (isapprox.(p[faces[:,4],3],i*hl,atol=ϵ,rtol=√eps()))
            trnf = [trnf;faces[idl,:]]
        end
    end

    # Add control node for top support
    p = [p;[0. 0. h]]

    # Create coordinate and connectivity arrays
    nn = size(p,1)
    # optscoor: ncel, nn, 0, !print
    optscoor = [3,nn,0,1]
    coor = [collect(1:nn) Printf.format.(Ref(Printf.Format("%.3f")),p)]
    nmat = size(tmat,1)
    nrnf = size(trnf,1)
    #optsconn: nel, 0, !print, 0, geoID, matID, ornID, conID, cadID
    optsmat = [nmat,0,1,0,mat_geoID,mat_matID,mat_ornID,conID,0]
    optsrnf = [nrnf,0,1,0,rnf_geoID,rnf_matID,0,0,0]
    if nne == 4
        cmat = [collect(1:nmat) tet_type*ones(Int,nmat) tmat]
        crnf = [collect(nmat+1:nmat+nrnf) tri_type*ones(Int,nrnf) trnf]
    elseif nne == 6
        cmat = [collect(1:nmat) pen_type*ones(Int,nmat) tmat]
        crnf = [collect(nmat+1:nmat+nrnf) tri_type*ones(Int,nrnf) trnf]  
    elseif nne == 8
        cmat = [collect(1:nmat) hex_type*ones(Int,nmat) tmat]
        crnf = [collect(nmat+1:nmat+nrnf) quad_type*ones(Int,nrnf) trnf]
    end

    # Reference node
    idrn = isapprox.(p[:,1],0.) .& isapprox.(p[:,2],0.) .& isapprox.(p[:,3],h)
    refn = collect(1:nn)[idrn][1]

    # Create sets of points in symmetry axes
    if half_x
        idx = collect(1:nn)[isapprox.(p[:,1],0.0)]
    end
    if half_y
        idy = collect(1:nn)[isapprox.(p[:,2],0.0)]
    end

    # Define file names
    if isnothing(filename)
        meshfile = joinpath(root,"mesh.txt")
        setfile = joinpath(root,"set.txt")
    else
        meshfile = joinpath(root,string("mesh_",filename,".txt"))
        setfile = joinpath(root,string("set_",filename,".txt"))
    end

    # Write coordinate and connectivity arrays to mesh.txt
    open(meshfile,"w") do io
        println(io,"connectivity")
        writedlm(io,optsmat',",")
        writedlm(io,cmat,",")
        println(io,"connectivity")
        writedlm(io,optsrnf',",")
        writedlm(io,crnf,",")
        println(io,"coordinates")
        writedlm(io,optscoor',",")
        writedlm(io,coor,",")
    end

    # Write sets
    open(setfile,"w") do io
        println(io,"define,element,set,set_mat")
        println(io,1," to ",nmat)
        println(io,"define,element,set,set_rnf")
        println(io,nmat+1," to ",nmat+nrnf)
        println(io,"define,node,set,set_topnode")
        println(io,nn)
        println(io,"define,node,set,set_refnode")
        println(io,refn)
        if half_x
            println(io,"define,node,set,set_nodesx0")
            nid = 1
            rid = 0
            for i in eachindex(idx)
                print(io,idx[nid]," ")
                nid += 1
                rid += 1
                if nid > length(idx)
                    break
                end
                if rid == 10
                    println(io,"c")
                    rid = 0
                end
            end
        end
        if half_y
            println(io,"define,node,set,set_nodesy0")
            nid = 1
            rid = 0
            for i in eachindex(idy)
                print(io,idy[nid]," ")
                nid += 1
                rid += 1
                if nid > length(idy)
                    break
                end
                if rid == 10
                    println(io,"c")
                    rid = 0
                end
            end
        end
    end

    # Plot mesh
    function fdp(q)
        if cs_type == "circular"
            d1 = dcylinder(q,0.,0.,0.,R,h)
            if half_x & !half_y
                d2 = dblock(q,0.,1.5*R,-1.5*R,1.5*R,-0.5*h,1.5*h)
                dp = dintersect(d1,d2)
            elseif half_y & !half_x
                d2 = dblock(q,-1.5*R,1.5*R,0.,1.5*R,-0.5*h,1.5*h)
                dp = dintersect(d1,d2)
            elseif half_x & half_y
                d2 = dblock(q,0.,1.5*R,0.,1.5*R,-0.5*h,1.5*h)
                dp = dintersect(d1,d2)
            else
                dp = d1
            end
        elseif cs_type == "annular"
            de = dcylinder(q,0.,0.,0.,R,h)
            di = dcylinder(q,0.,0.,0.,Ri,h)
            d1 = ddiff(de,di)
            if half_x & !half_y
                d2 = dblock(q,0.,1.5*R,-1.5*R,1.5*R,-0.5*h,1.5*h)
                dp = dintersect(d1,d2)
            elseif half_y & !half_x
                d2 = dblock(q,-1.5*R,1.5*R,0.,1.5*R,-0.5*h,1.5*h)
                dp = dintersect(d1,d2)
            elseif half_x & half_y
                d2 = dblock(q,0.,1.5*R,0.,1.5*R,-0.5*h,1.5*h)
                dp = dintersect(d1,d2)
            else
                dp = d1
            end
        elseif cs_type == "rectangular"
            if half_x & !half_y
                dp = dblock(q,0.,a/2,-b/2,b/2,0.,h)
            elseif half_y & !half_x
                dp = dblock(q,-a/2,a/2,0.,b/2,0.,h)
            elseif half_x & half_y
                dp = dblock(q,0.,a/2,0.,b/2,0.,h)
            else
                dp = dblock(q,-a/2,a/2,-b/2,b/2,0.,h)
            end
        end
        return dp
    end

    if plotmesh
        figure(figsize=(8,8))
        mplot3d(p,tmat,fdp=fdp,idtype=plotid);
        figure(figsize=(8,8))
        mplot3d(p,trnf,idtype=plotid);
    end
    return p, tmat, trnf
end