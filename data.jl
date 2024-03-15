function read_file(filename)

    filepath = string(@__DIR__, "/instances/", filename)

    J = Int64[]
    w = Int64[]
    E = Vector{Int64}[]
    
    line_number = 0
    items_amount = -1
    bin_size = -1
    
    open(filepath) do file
        for line in eachline(file)

            line_number += 1

            if line_number == 1
                m = match(r"^(\d+) (\d+)$", line)
                if m !== nothing
                    items_amount = parse(Int, m.captures[1])
                    bin_size = parse(Int, m.captures[2])
                end
                continue
            end

            l = parse.(Int, split(line, " "))

            j = l[1]
            w_j = l[2]
            conflicts = l[3:end]

            push!(J, j)
            push!(w, w_j)
            for i in conflicts
                push!(E, [j, i])
            end
        end
        sort!(E)
    end

    return J, w, E
end

function write_to_file(J, w, E, filepath)

    conflicts = Vector{Int64}[Int64[] for i in J]
    for e in E
        push!(conflicts[e[1]], e[2])
        push!(conflicts[e[2]], e[1])
    end

    open(filepath, "w") do file
        for (n, j) in enumerate(J)
            
            line = string()
            line = string()
        end
        write(file, line)
    end
end


