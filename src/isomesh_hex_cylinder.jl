#=
3D Hexahedral Mesh Generator Reading External Mesh

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 3d mesh for a cylinder section using hexahedral elements reading an external unstructured quad mesh. 
=#

#----------------------------INPUT----------------------------------------
#=
h:          height of cylinder
nz:         number of elements in z direction
root:       root of external mesh file
filename:   file name of external mesh
z0:         zero z coordinate of the block
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y,z) coordinates. Size: nn x 3
t:      array of connectivity. Size: nel x 6
=#

#----------------------------START----------------------------------------

function hex_cylinder(h::Float64,
    nz::Int,
    root::String,
    filename::String;
    z0::Float64=0.)::Tuple{Matrix{Float64},Matrix{Int}}

    pcs, tcs = readmesh(root,filename)
    ncs = size(pcs,1)
    p = [pcs z0*ones(ncs)]
    t = zeros(0,8)

    for i in 1:nz
        p = [p;[pcs (z0+i*h/nz)*ones(ncs)]]
        t = [t;[tcs.+(i-1)*ncs tcs.+i*ncs]]
    end

    return p, t
end