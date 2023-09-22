#=
3D Hex Mesh Generator

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 3d mesh for a block section using hexahedral elements based on Hermitian interpolation to get an aspect ratio close to 1 at the borders, and a given aspect ratio at the center of the mesh.
=#

#----------------------------INPUT----------------------------------------
#=
REQUIRED
d:      depth of block
w:      width of block
h:      height of block
nz:     number of elements in z direction
er:     target element aspect ratio at center

OPTIONAL
x0:     zero x coordinate of the block
y0:     zero y coordinate of the block
z0:     zero z coordinate of the block
half_x:     model only half of the bearing in x (defaults to false)
half_y:     model only half of the bearing in y (defaults to false)
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y,z) coordinates. Size: nn x 3
t:      array of connectivity. Size: nel x 8
=#

#----------------------------START----------------------------------------

function hex_block(d::Float64,
    w::Float64, 
    h::Float64,
    nz::Int,
    er::Int;
    x0::Float64=0., 
    y0::Float64=0.,
    z0::Float64=0.,
    half_x::Bool=false,
    half_y::Bool=false)::Tuple{Matrix{Float64},Matrix{Int}}

    xleft = 0.
    xright = 1.
    if (er != 1) & (!half_x)
        xξleft = 1/er
        xξright = 1/er
    elseif (er != 1) & (half_x)
        xξleft = 1
        xξright = 1/er
    else
        xξleft = 1
        xξright = 1
    end
    ybottom = 0.
    ytop = 1.
    if (er != 1) & (!half_y)
        yηbottom = 1/er
        yηtop = 1/er
    elseif (er != 1) & (half_y)
        yηbottom = 1
        yηtop = 1/er
    else
        yηbottom = 1
        yηtop = 1
    end
    H(β) = [(1 + 2*β)*(1 - β)^2,(1 - 2*(β - 1))*β^2,β*(1 - β)^2,(β - 1)*β^2]
    xint(ξ) = [xleft,xright,xξleft,xξright]'*H(ξ)
    if half_x; nx = round(Int,d/2/h*nz/er); else; nx = round(Int,d/h*nz/er); end
    nx += isodd(nx)*1
    ξ = 1/nx*collect(0:nx)
    x = (d - half_x*d/2)*xint.(ξ) .+ x0 .+ half_x*d/2
    yint(η) = [ybottom,ytop,yηbottom,yηtop]'*H(η)
    if half_y; ny = round(Int,w/2/h*nz/er); else; ny = round(Int,w/h*nz/er); end
    ny += isodd(ny)*1
    η = 1/ny*collect(0:ny)
    y = (w - half_y*w/2)*yint.(η) .+ y0 .+ half_y*w/2
    dz = h*1/nz
    z = dz*collect(0:nz) .+ z0
    pp = Vector{Array}(undef,3)
    pp[1] = [i for i in x, j in y, k in z]
    pp[2] = [j for i in x, j in y, k in z]
    pp[3] = [k for i in x, j in y, k in z]
    p = [vec(pp[1]) vec(pp[2]) vec(pp[3])]
    nnxy = (nx+1)*(ny+1)
    nel = nx*ny*nz
    t = Matrix{Int}(undef,nel,8)
    iel = 1
    n1 = 1
    for z in 1:nz
        for j in 1:ny
            for i in 1:nx
                ttemp = [n1,n1+1,n1+nx+2,n1+nx+1]
                t[iel,:] = [ttemp;ttemp.+nnxy]
                iel += 1
                n1 += 1
            end
            n1 += 1
        end
        n1 += nx + 1
    end
    return p, t
end