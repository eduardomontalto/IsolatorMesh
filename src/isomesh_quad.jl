#=
2D Quad Mesh Generator

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function generates a 2d mesh for a rectangular area using quad elements based on Hermitian interpolation to get an aspect ratio close to 1 at the borders, and a given aspect ratio at the center of the rectangular mesh.
=#

#----------------------------INPUT----------------------------------------
#=
w:      width of rectangular area
h:      height of rectangular area
ny:     number of elements in y direction
er:     target element aspect ratio (width/height) at center
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y) coordinates. Size: nn x 2
t:      array of connectivity. Size: nel x 4 
=#

#----------------------------START----------------------------------------

function quad(w::Float64, 
    h::Float64,
    ny::Int,
    er::Int;
    x0::Float64=0.,
    y0::Float64=0.,
    varsize::Bool=true)::Tuple{Matrix{Float64},Matrix{Int}}

    xleft = 0.
    xright = 1.
    if varsize
        xξleft = 1/er
        xξright = 1/er
    else
        xξleft = 1
        xξright = 1
    end
    H(ξ) = [(1 + 2*ξ)*(1 - ξ)^2,(1 - 2*(ξ - 1))*ξ^2,ξ*(1 - ξ)^2,(ξ - 1)*ξ^2]
    xint(ξ) = [xleft,xright,xξleft,xξright]'*H(ξ)
    nx = round(Int,w/h*ny/er) 
    ξ = 1/nx*collect(0:nx)
    x = w*xint.(ξ) .+ x0
    dy = h*1/ny
    y = dy*collect(0:ny) .+ y0
    x = x*ones(ny+1)'
    y = ones(nx+1)*y'
    nn = (nx+1)*(ny+1)
    p = [reshape(x,nn) reshape(y,nn)]
    nel = nx*ny
    t = Matrix{Int}(undef,nel,4)
    iel = 1
    n1 = 1
    for j=1:ny
        for i=1:nx
            t[iel,:] = [n1,n1+1,n1+nx+2,n1+nx+1] 
            iel += 1
            n1 += 1
        end
        n1 += 1
    end
    return p, t
end

function quad(w::Float64, 
    h::Float64,
    ny::Int,
    er::Float64;
    x0::Float64=0.,
    y0::Float64=0.,
    varsize::Bool=true)::Tuple{Matrix{Float64},Matrix{Int}}

    xleft = 0.
    xright = 1.
    if varsize
        xξleft = 1/er
        xξright = 1/er
    else
        xξleft = 1
        xξright = 1
    end
    H(ξ) = [(1 + 2*ξ)*(1 - ξ)^2,(1 - 2*(ξ - 1))*ξ^2,ξ*(1 - ξ)^2,(ξ - 1)*ξ^2]
    xint(ξ) = [xleft,xright,xξleft,xξright]'*H(ξ)
    nx = round(Int,w/h*ny/er) 
    ξ = 1/nx*collect(0:nx)
    x = w*xint.(ξ) .+ x0
    dy = h*1/ny
    y = dy*collect(0:ny) .+ y0
    x = x*ones(ny+1)'
    y = ones(nx+1)*y'
    nn = (nx+1)*(ny+1)
    p = [reshape(x,nn) reshape(y,nn)]
    nel = nx*ny
    t = Matrix{Int}(undef,nel,4)
    iel = 1
    n1 = 1
    for j=1:ny
        for i=1:nx
            t[iel,:] = [n1,n1+1,n1+nx+2,n1+nx+1] 
            iel += 1
            n1 += 1
        end
        n1 += 1
    end
    return p, t
end