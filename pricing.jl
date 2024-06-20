struct Label
    rcost::Float64 
    weight::Int64 # current weight
    last_item_added::Int # item added in this label
    prev_lab::Vector{Label} # last label
    next_conflics::BitVector # conflics that will be found by moving forward (includes weight conflicts)

    # last_arc::Vector{Int} # last arc in the route, in the (entry, exit) format
    # prev_lab::Vector{Label} # last label
end

"Returns true if l1 dominates l2"
function Base.isless(l1::Label, l2::Label)

    # l1 dominates l2 if:
    #   the possibilities set of l1 *at least contains* the possibilities set of l2
    #   l1 has smaller reduced cost

    if l1.weight <= l2.weight
        if l1.rcost <= l2.rcost
            return all(l1.next_conflics .<= l2.next_conflics)
        else
            return false
        end
    else
        return false
    end
end


# base data
len_J = 10
J = [i for i in 1:len_J]

# reduced costs
# rc = [i for i in J]
rc = rand(-10:10, 10)

positive_rcost = Bool[i > 0 for i in rc]

w = [i*10 for i in J]
W = 120

# binarized_E[i][j] = 1 if (i, j) ∈ E
# binarized_E = Vector{Int64}[Int64[0 for j in 1:len_J] for i in 1:len_J]
# binarized_E = BitVector[BitVector(undef, len_J) for i in J]
binarized_E = BitVector[falses(len_J) for i in J]


for (i, j) in translated_E
    binarized_E[i][j] = true
    binarized_E[j][i] = true
end

buckets = Vector{Vector{Label}}[Label[] for i in J]

to_extend = BinaryMinHeap{Label}()
for i in 1:len_J
    label = Label(c[i], w[i], i, Label[], deepcopy(binarized_E[i][i+1:end]))

    push!(buckets[i], label)

    if i < len_J
        push!(to_extend, label)
    end
end

trash = Dict{Label, Nothing}()

while !isempty(to_extend)
    curr_label = pop!(to_extend)

    # if the label is marked for deletion, delete and continue
    if haskey(trash, curr_label)
        delete!(trash, curr_label)
        continue
    end


    for i in curr_label.last_item_added+1:len_J

        # if the item has positive reduced cost, skip it
        if positive_rcost[i]
            continue
        end

        new_weight = curr_label.weight + w[i]
        if new_weight > W
            continue
        end

        new_label = Label(
            curr_label.rcost + rc[i], 
            new_weight, 
            i, 
            Label[curr_label], 
            curr_label.next_conflics[i-curr_label.last_item_added+1:end] .|| binarized_E[i+1:end],
        )



        # check for domination
        dominated = Dict{Label, Nothing}()
        clean_bucket = false
        for label in buckets[i]
            
            # is the new label dominated?
            if label < new_label    
                new_label_is_dominated = true
                break
            end

            # does the new label dominate any label?
            if new_label < label

                # add dominated to trash
                trash[label] = nothing

                # register dominated for later deletion
                dominated[label] = nothing
                clean_matrix = true
            end
        end

        # remove the labels dominated by the new label, if necessary
        if clean_matrix
            deleteat!(buckets[i], findall(x -> x in keys(dominated), buckets[i]))
        end

        if new_label_is_dominated
            continue
        else # if the new label isn't dominated
            push!(buckets[i], new_label)
            push!(to_extend, new_label)
        end 
    end
end

min_rcost = Inf
best_label = nothing
for bucket in buckets
    for label in bucket
        if label.rcost < min_rcost
            min_rcost = label.rcost
            best_label = label
        end
    end
end

new_bin = falses(len_J)
if min_rcost > 0 && min_rcost < Inf
    label = best_label
    
    done = false
    while !done
        if isempty(label.prev_lab) 
            done = true
        else
            new_bin[label.last_item_added] = true
            label = label.prev_lab[1]
        end
    end
end








mat = Matrix{Vector{Label}}(undef, Q+1, length(N))

# reset the matrix
for c in N
    for alpha in 0:Q 
        mat[alpha+1,c-1] = Label[]
    end
end

to_extend = BinaryMinHeap{Label}()
for c in N
    d = delivery(data, c)
    p = pickup(data, c)

    alpha = p
    beta = max(d, alpha)

    label = Label(costs[1,c], alpha, beta, Int[1, c], Label[])
    push!(to_extend, label)
    push!(mat[alpha+1, c-1], label)
end

trash = Dict{Label, Nothing}()

while !isempty(to_extend)
    curr_label = pop!(to_extend)

    if haskey(trash, curr_label)
        delete!(trash, curr_label)
        continue
    end
    
    for c in N
        if c == curr_label.last_arc[2]
            continue
        end

        alpha = curr_label.alpha + pickup(data, c)
        if alpha > Q
            continue
        end
        
        beta = max(curr_label.beta + delivery(data, c), alpha) 
        if beta > Q
            continue
        end

        rcost = curr_label.rcost + costs[curr_label.last_arc[2], c]
        
        # It will become an loop in the list of labels
        new_label = Label(rcost, alpha, beta, Int[curr_label.last_arc[2], c], Label[curr_label])

        new_label_is_dominated = false

        # check for domination
        dominated = Dict{Label, Nothing}()
        clean_matrix = false
        for label in mat[alpha+1, c-1]
            
            # is the new label dominated?
            if label < new_label    
                new_label_is_dominated = true
                break
            end

            # does the new label dominate any label?
            if new_label < label

                # add dominated to trash
                trash[label] = nothing

                # register dominated for later deletion
                dominated[label] = nothing
                clean_matrix = true
            end
        end

        # remove the labels dominated by the new label, if necessary
        if clean_matrix
            deleteat!(mat[alpha+1, c-1], findall(x -> x in keys(dominated), mat[alpha+1, c-1]))
        end

        if new_label_is_dominated
            continue
        else # if the new label isn't dominated
            push!(mat[alpha+1, c-1], new_label)
            push!(to_extend, new_label)
        end
    end
end

min_rcost = Inf
sols = Tuple{Float64, Vector{JuMP.VariableRef}, Vector{Float64}}[]
for c in N
    for alpha in 0:Q
        for label in mat[alpha+1, c-1]
            sol = Dict{Tuple{Int, Int}, Int}()

            # from the current vertex, go backwards until the depot is reached
            done = false
            curlabel = label
            while !done # stop at the begginning
                arc = tuple(curlabel.last_arc...)
                sol[arc] = get(sol, arc, 0) + 1                        
                
                if arc[1] == 1
                    done = true
                else
                    curlabel = curlabel.prev_lab[1] # it's a reference contained in an array
                end
            end

            # add final arc, finishing the route 
            arc = (c, 1)
            sol[arc] = get(sol, arc, 0) + 1

            # get route reduced cost, and compare with global minimum
            rcost = label.rcost + costs[c, 1]
            min_rcost = min(min_rcost, rcost)
            # @show rcost, sol
            
            # Create the solution (send only variables with non-zero values)
            solvars = JuMP.VariableRef[]
            solvals = Float64[]
            for (arc, val) in sol
                push!(solvars, x[spid, arc])
                push!(solvals, val)
            end
        
            push!(sols, (rcost, solvars, solvals))
        end
    end
end
@show min_rcost
sort!(sols, by = x -> x[1])
for sol in sols
    @show sol
    MOI.submit(cvrp, BD.PricingSolution(cbdata), sol...)
    break
end

