#=
Reader for External Mesh

University of California, Berkeley
Department of Civil & Environmental Engineering
CE299 Independent Study
Instructor: Dimitrios Konstantinidis
Student: Eduardo Montalto
=#

#----------------------------DESCRIPTION----------------------------------
#= This function reads an external unstructured quad mesh. 
=#

#----------------------------INPUT----------------------------------------
#=
root:       root of external mesh file
filename:   file name of external mesh
=#

#----------------------------OUTPUT---------------------------------------
#=
p:      array of (x,y,z) coordinates. Size: nn x 3
t:      array of connectivity. Size: nel x 6
=#

#----------------------------START----------------------------------------

function readmesh(root::String,
    filename::String)::Tuple{Matrix{Float64},Matrix{Int}}
    
    p = zeros(0,2)
    t = zeros(Int,0,4)
    lines = readlines(joinpath(root,filename))
    lineid = 1
    connectivity = false
    coordinates = false
    while true
        lineid += 1
        if lineid > length(lines); break; end;
        filteredline = rstrip(lstrip(lines[lineid]))
        if filteredline == "connectivity"
            connectivity = true
            lineid += 2
        elseif filteredline == "coordinates"
            connectivity = false
            coordinates = true
            lineid += 2
        elseif filteredline == "loadcase            job1"
            break
        end
        if connectivity
            linevec = filter(x -> x !="",split(lines[lineid]," "))[3:6]
            t = [t; parse.(Int,linevec)']
        end
        if coordinates
            linevec = [lines[lineid][11:30],lines[lineid][31:50]]
            for j in 1:2
                if contains(linevec[j],'+')
                    linevec[j] = string(rsplit(linevec[j],"+",limit=2)[1],"E+",rsplit(linevec[j],"+",limit=2)[2])
                elseif contains(linevec[j],'-')
                    linevec[j] = string(rsplit(linevec[j],"-",limit=2)[1],"E-",rsplit(linevec[j],"-",limit=2)[2])
                end
            end
            p = [p; parse.(Float64,linevec)']
        end

    end

    return p, t
end