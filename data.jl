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
            filter!(x -> x < items_amount, conflicts)


            push!(J, j)
            push!(w, w_j)
            for i in conflicts
                push!(E, [j, i])
            end
        end
        sort!(E)
    end

    return J, w, E, bin_size
end

function write_to_file(J, w, E, W, filepath)

    conflicts = Vector{Int64}[Int64[] for i in J]
    for e in E
        if e[1] < e[2]
            push!(conflicts[e[1]], e[2])
        else
            push!(conflicts[e[2]], e[1])
        end
    end
    for c in conflicts
        sort!(c)
    end

    open(filepath, "w") do file
        
        s = string(length(J), " ", W, "\n")

        for j in J
            s *= string(j, " ", w[j])
            for i in conflicts[j] 
                s *= string(" ", i)
            end
            
            s *= "\n"
        end

        write(file, s)
    end
end

function check_solution_viability(sol, J, w, E, W)

    println("Bins amount: $(length(sol))")

    missing_item = Bool[true for j in J]

    binarized_E = BitVector[falses(length(J)) for j in J]
    for (i, j) in E
        binarized_E[i][j] = true
        binarized_E[j][i] = true
    end
    
    # register of bin conflicts
    conflicts = BitVector[falses(length(J)) for bin in sol]

    # binarized version of the bin
    bins_binarized = BitVector[falses(length(J)) for bin in sol]

    for (i, bin) in enumerate(sol)
        for j in bin
            if !missing_item[j]
                println("item $(j) was found more than once")
            end

            missing_item[j] = false
            bins_binarized[i][j] = true

            conflicts[i] = conflicts[i] .|| binarized_E[j]
        end

        if any(conflicts[i] .&& bins_binarized[i])
            println("bin $(i) has conflicting items")
        end

        bin_weight = sum(bins_binarized[i] .* w)
        if bin_weight > W
            println("bin $(i) is too heavy: $(bin_weight) > $(W)")
        end
    end

    if any(missing_item)
        for (j, m) in enumerate(missing_item)
            if m
                println("missing item $(j)")
            end
        end
    end
end

