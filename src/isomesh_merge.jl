#=
Mesh Merger

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function takes two meshes and merges them by deleting repeated nodes. Numbering of the first mesh is maintained and numbering of the second mesh is altered for consistency 
=#

#----------------------------INPUT----------------------------------------
#=
REQUIRED
p1:     array of (x,y,(z)) coordinates of main mesh
p2:     array of (x,y,(z)) coordinates of second mesh
t2:     array of connectivity of second mesh

OPTIONAL
ϵ:      tolerance to check for repeated nodes (defaults to 1e-12)
=#

#----------------------------OUTPUT---------------------------------------
#=
p1:     array of (x,y,(z)) coordinates of merged mesh. Size: nn x ndim
p2:     array of (x,y,(z)) coordinates of points unique to second mesh. 
        Size: nn x ndim
t2:     array of connectivity of second mesh. Size: nel x nen
=#

#----------------------------START----------------------------------------

function merge_mesh!(p1::Matrix{Float64},
    p2::Matrix{Float64},
    t2::Matrix{Int};
    ϵ::Float64=1e-5)::Tuple{Matrix{Float64},Matrix{Float64},Matrix{Int}}
    
    nn1 = size(p1,1)
    nn2 = size(p2,1)
    idp = BitArray(ones(nn2))
    t2temp = t2 .+ nn1

    for j in 1:nn2
        for i in 1:nn1
            if all(isapprox.(p2[j,:],p1[i,:],atol=ϵ,rtol=√eps()))
                t2temp[t2 .== j] .= copy(i)
                idp[j] = false
                t2temp[t2 .> j] .+= -1
                break
            end
        end
    end

    p2 = p2[idp,:]
    p1 = [p1;p2]
    t2 = copy(t2temp)
    return p1, p2, t2
end