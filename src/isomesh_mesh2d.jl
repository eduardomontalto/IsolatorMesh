#=
2D Mesh Generator for Strip Isolators in MSC Marc

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 2d mesh for a rectangular isolator to input in MSC Marc. Central section uses a structured mesh with Hermite interpolation to have more elements towards the ends of the section. Boundary sections use an unstructured mesh with triangular elements created with DistMesh2D.
=#

#----------------------------INPUT----------------------------------------
#=
REQUIRED
a:          width of isolator
h:          total height of rubber layers
nl:         number of rubber layers
epl:        number of elements per rubber layer
ear:        average element aspect ratio (maximum is larger e.g. 5.5 for ear=4)

OPTIONAL
eplb:       number of elements per rubber layer on sides (defaults to epl)
ab:         width of boundary sections with tri elements (defaults to h/nl)
cen_vsize:  variable size of elements in central section (defaults to true)
boun_vsize: variable size of elements in boundary section (defaults to true)
boun_nne:   number of nodes per element in boundary elements (defaults to 3)
quad_type:  quad element type (11 disp., 80 mixed) (defaults to 11)
tri_type:   tri element type (6 disp., 155 mixed) (defaults to 80)
lin_type:   linear element type (9 truss, 165 2-node rebar) (defaults to 165)
cen_conID:  contact body ID of quad elements (defaults to 1)
boun_conID: contact body ID of tri elements (defaults to 2)
mat_geoID:  geometry ID of rubber matrix elements (defaults to 1)
mat_matID:  material ID of rubber matrix elements (defaults to 1)
mat_ornID:  orientation ID of rubber matrix elements (defaults to 1)
rnf_geoID:  geometry ID of reinforcement elements (defaults to 2)
rnf_matID:  material ID of reinforcement elements (defaults to 3)
root:       root directory where output files will be saved
filename:   name of output text files
plotmesh:   plot generated mesh
plotid:     identify element types in plot
=#

#----------------------------OUTPUT---------------------------------------
#=
FILES
mesh.txt:   text file containing connectivity and coordinates
set.txt:    text file containing element sets
=#

#----------------------------START----------------------------------------

function mesh2d(a::Float64,
    h::Float64,
    nl::Int,
    epl::Int,
    ear::Int;
    ab::Float64=h/nl, 
    cen_vsize::Bool=true, 
    boun_vsize::Bool=true, 
    boun_nne::Int=3, 
    quad_type::Int=11, 
    tri_type::Int=155, 
    lin_type::Int=165, 
    cen_conID::Int=1, 
    boun_conID::Int=2, 
    mat_geoID::Int=1, 
    mat_matID::Int=1, 
    mat_ornID::Int=1, 
    rnf_geoID::Int=2, 
    rnf_matID::Int=3, 
    root=pwd(), 
    filename=nothing, 
    plotmesh::Bool=false, 
    plotid::Bool=false)

    # Create mesh of central section
    hl = h/nl                                   
    pcen, tcen = quad(a-2*ab, h, epl*nl, ear, x0=-a/2+ab, varsize=cen_vsize)
    
    # Create mesh of boundary sections
    if boun_nne == 3
        
        eplb = (1+boun_vsize)*epl
        pbounl, tbounl = tri(ab,hl,[eplb, epl],x0=-a/2)
        ncrn = size(tbounl,1)
        pbounl1 = copy(pbounl)
        tbounl1 = copy(tbounl)
        for i in 2:nl
            if iseven(i)
                pbounl2 = pbounl1 .- [0. hl]
                pbounl2 .*= [1. -1.]
                pbounl2 .+= [0. (i-1)*hl]
            else
                pbounl2 = pbounl1 .+ [0. (i-1)*hl]
            end
            tbounl2 = copy(tbounl1)
            pbounl, _, tbounl2 = merge_mesh!(pbounl, pbounl2, tbounl2)
            tbounl = [tbounl;tbounl2]
        end
        
        tbounl = check_direction(pbounl,tbounl)
        pbounr = [-pbounl[:,1] pbounl[:,2]]
        tbounr = copy(tbounl)
        tbounr = check_direction(pbounr,tbounr)

    elseif boun_nne == 4                                 
        pbounl, tbounl = quad(ab, h, epl*nl, 1, x0=-a/2, varsize=false)
        ncrn = size(tbounl,1)÷nl
        pbounr, tbounr = quad(ab, h, epl*nl, 1, x0=a/2-ab, varsize=false)
    end

    # Merge meshes
    p, _, tbounl = merge_mesh!(pcen, pbounl, tbounl)
    p, _, tbounr = merge_mesh!(p, pbounr, tbounr)
    tboun = [tbounl;tbounr]

    # Create connectivity of linear elements
    ϵ = 1e-5
    tlin = zeros(Int,0,2)

    for i in 1:nl-1
        idl = isapprox.(p[sortperm(p[:,1]),2],i*hl,atol=ϵ,rtol=√eps())
        tlin2 = sortperm(p[:,1])[idl]
        tlin2 = [tlin2[1:end-1] tlin2[2:end]]
        tlin = [tlin; tlin2]
    end

    # Add control node for top support
    p = [p;[0. h]]

    # Create coordinate and connectivity arrays
    nn = size(p,1)
    # optscoor: ncel, nn, 0, !print
    optscoor = [2,nn,0,1]
    coor = [collect(1:size(p,1)) Printf.format.(Ref(Printf.Format("%.3f")),p)]
    ncen = size(tcen,1)
    nboun = size(tboun,1)
    nbounl = size(tbounl,1)
    nbounr = size(tbounr,1)
    nlin = size(tlin,1)
    #optsconn: nel, 0, !print, 0, geoID, matID, ornID, conID, cadID
    optscen = [ncen,0,1,0,mat_geoID,mat_matID,mat_ornID,cen_conID,0]
    optsbounl = [nbounl,0,1,0,mat_geoID,mat_matID,mat_ornID,boun_conID,0]
    optsbounr = [nbounr,0,1,0,mat_geoID,mat_matID,mat_ornID,boun_conID+1,0]
    optslin = [nlin,0,1,0,rnf_geoID,rnf_matID,0,0,0]
    ccen = [collect(1:ncen) quad_type*ones(Int,ncen) tcen]
    if boun_nne == 3
        cbounl = [collect(ncen+1:ncen+nbounl) tri_type*ones(Int,nbounl) tbounl]
        cbounr = [collect(ncen+nbounl+1:ncen+nboun) tri_type*ones(Int,nbounr) tbounr]
    elseif boun_nne == 4
        cbounl = [collect(ncen+1:ncen+nbounl) quad_type*ones(Int,nbounl) tbounl]
        cbounr = [collect(ncen+nbounl+1:ncen+nboun) quad_type*ones(Int,nbounr) tbounr]
    end
    clin = [collect(ncen+nboun+1:ncen+nboun+nlin) lin_type*ones(Int,nlin) tlin]

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
        writedlm(io,optscen',",")
        writedlm(io,ccen,",")
        println(io,"connectivity")
        writedlm(io,optsbounl',",")
        writedlm(io,cbounl,",")
        println(io,"connectivity")
        writedlm(io,optsbounr',",")
        writedlm(io,cbounr,",")
        println(io,"connectivity")
        writedlm(io,optslin',",")
        writedlm(io,clin,",")
        println(io,"coordinates")
        writedlm(io,optscoor',",")
        writedlm(io,coor,",")
    end

    # Write sets
    open(setfile,"w") do io
        println(io,"define,element,set,set_cen")
        println(io,1," to ",ncen)
        println(io,"define,element,set,set_boun")
        println(io,ncen+1," to ",ncen+nboun)
        println(io,"define,element,set,set_rnf")
        println(io,ncen+nboun+1," to ",ncen+nboun+nlin)
        println(io,"define,node,set,set_topnode")
        println(io,nn)
        println(io,"define,element,set,set_corner1")
        println(io,ncen+1," to ",ncen+ncrn)
        println(io,"define,element,set,set_corner2")
        println(io,ncen+(nboun÷2)-ncrn+1," to ",ncen+(nboun÷2))
        println(io,"define,element,set,set_corner3")
        println(io,ncen+(nboun÷2)+1," to ",ncen+(nboun÷2)+ncrn)
        println(io,"define,element,set,set_corner4")
        println(io,ncen+nboun-ncrn+1," to ",ncen+nboun)
    end

    # Plot mesh
    if plotmesh
        figure(figsize=(12,8))
        mplot2d(p,tcen,idtype=plotid)
        mplot2d(p,tboun,idtype=plotid)
        mplot2d(p,tlin,idtype=plotid);
    end
    return p, tcen, tboun, tlin
end