#=
Mesh Connection Check

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function checks that the connectivity of a triangular mesh has been established in a right handed system.
=#

#----------------------------INPUT----------------------------------------
#=
p:      array of (x,y) coordinates of mesh
t:      array of connectivity of mesh

=#

#----------------------------OUTPUT---------------------------------------
#=
t:      array of connectivity of mesh
=#

#----------------------------START----------------------------------------

function check_direction!(p::Matrix{Float64},
    t::Matrix{Int})::Matrix{Int}

    nel = size(t,1)
    nne = size(t,2)
    if nne == 3
        s1 = p[t[:,2],:]-p[t[:,1],:]
        s2 = p[t[:,3],:]-p[t[:,1],:]
        if size(p,2) == 2
            s1 = [s1 zeros(nel)]
            s2 = [s2 zeros(nel)]
        end
        for i in 1:nel
            dir = cross(s1[i,:],s2[i,:])[3]
            if dir < 0.
                t[i,:] = [t[i,1] t[i,3] t[i,2]]
            end
        end
    elseif nne == 4
        s1 = p[t[:,2],:]-p[t[:,1],:]
        s2 = p[t[:,3],:]-p[t[:,1],:]
        s3 = p[t[:,4],:]-p[t[:,1],:]
        for i in 1:nel
            dir = dot(cross(s1[i,:],s2[i,:]),s3[i,:])
            if dir < 0.
                t[i,:] = [t[i,1] t[i,3] t[i,2] t[i,4]]
            end
        end
    end
    return t
end