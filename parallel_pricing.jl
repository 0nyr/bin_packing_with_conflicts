using DataStructures

struct Label
    rcost::Float64 # less is better
    fcost::Vector{Float64} # reduced cost considering cut violations (less is better)
    weight::Int64 # current weight
    last_item_added::Int # item added in this label
    items::BitVector # items in the bin so far
    conflicts::BitVector # conflics that will be found by moving forward
    m::Vector{Int64} # amount of custumers involved in cut i for i in cuts: |S_i ∩ V(L)|
    sigma_ref::Vector{Float64} # necessary for dominance criteria...
end

# auxiliary stuff
SPINLOCK = Threads.SpinLock()

"Returns true if l1 dominates l2"
function Base.isless(l1::Label, l2::Label)

    # l1 dominates l2 if:
    #   the possibilities set of l1 *at least contains* the possibilities set of l2
    #   l1 has smaller reduced cost (considering cut related effects)

    if l1.weight <= l2.weight
        if all(l1.conflicts .<= l2.conflicts)

            sum_sigma_q_in_Q = 0
            if length(l1.m) != 0
                
                # get sigmas smaller than zero
                negative_sigma = l1.sigma_ref .< -1e-4
            
                if any(negative_sigma)
                    
                    for (q, is_negative) in enumerate(negative_sigma) # start by checking sigma, most likely positive 
                        if is_negative

                            # tau_l1 > tau_l2 ?
                            if l1.m[q] % 2 > l2.m[q] % 2 
                                sum_sigma_q_in_Q += l1.sigma_ref[q]     
                            end
                        end
                    end
                end
            end
            return l1.fcost[1] - sum_sigma_q_in_Q <= l2.fcost[1]
            # return l1.rcost - sum_sigma_q_in_Q <= l2.rcost
            # return l1.fcost[1] <= l2.fcost[1]
        else
            return false
        end
    else
        return false
    end

end

"Updates the amount of customers involved in each cut after the label visits customer i"
function update_m(label, i, cuts_binary_data)
    for (k, cut_k) in enumerate(cuts_binary_data)
        if cut_k[i]
            label.m[k] += 1
        end
    end
end

"Updates the label final cost (the reduced cost considering cut violations)"
function update_fcost(label::Label)
    label.fcost[1] = label.rcost - sum(label.sigma_ref.*floor.(label.m / 2))
end

"Dynamic programming (labelling) price"
function dp_price(J::Vector{Int64}, len_J::Int64, rc::Vector{Float64}, sigma::Vector{Float64}, positive_rcost::Vector{Bool}, w::Vector{Int64}, binarized_E::Vector{BitVector}, W::Int64, subset_row_cuts::Vector{Vector{Int64}}, cuts_binary_data::Vector{BitVector}, sr_k::Vector{Int64}; verbose=3, epsilon=1e-4)

    # fast_labelling = false

    # auxiliary data structure
    # sigma_multiplier = Float64[0.0 for i in subset_row_cuts]

    buckets = Vector{Label}[Label[] for i in J]

    to_extend = BinaryMinHeap{Label}()
    for i in 1:len_J
        if !positive_rcost[i]
            continue
        end

        # label = Label(rc[i], w[i], i, Label[], deepcopy(binarized_E[i][i+1:end]))
        # label = Label(1-rc[i], w[i], i, Label[], deepcopy(binarized_E[i]))
        label = Label(1-rc[i], Float64[0.0], w[i], i, falses(len_J), deepcopy(binarized_E[i]), Int64[0 for _ in subset_row_cuts], sigma)
        label.items[i] = true
        update_m(label, i, cuts_binary_data)
        update_fcost(label)

        push!(buckets[i], label)

        if i < len_J
            push!(to_extend, label)
        end
    end

    # was_extended = false
    trash = Dict{Label, Nothing}()
    # println("starting extensions")
    while !isempty(to_extend)

        # println(to_extend)

        curr_label = pop!(to_extend)

        # println(curr_label)

        # if the label is marked for deletion, delete and continue
        if haskey(trash, curr_label)
            delete!(trash, curr_label)
            continue
        end

        was_extended = false
        Threads.@threads for i in curr_label.last_item_added+1:len_J

            # if the item has negative reduced cost, skip it
            if !positive_rcost[i]
                continue
            end

            # if the item is in conflict with the label, skip it
            if curr_label.conflicts[i]
                continue
            end

            new_weight = curr_label.weight + w[i]
            if new_weight > W
                continue
            end

            # println("here")
            # new_next_conflicts =  curr_label.conflicts[i-curr_label.last_item_added+1:end] .|| binarized_E[i][i+1:end]
            new_next_conflicts = deepcopy(curr_label.conflicts)
            new_next_conflicts[i+1:end] =  curr_label.conflicts[i+1:end] .|| binarized_E[i][i+1:end]
            # println(new_next_conflicts)

            new_label = Label(
                curr_label.rcost - rc[i], 
                Float64[0.0],
                new_weight, 
                i, 
                # Label[curr_label], 
                deepcopy(curr_label.items), 
                new_next_conflicts,
                # curr_label.conflicts[i-curr_label.last_item_added+1:end] .|| binarized_E[i+1:end],
                deepcopy(curr_label.m),
                sigma,
            )

            new_label.items[i] = true
            update_m(new_label, i, cuts_binary_data)
            update_fcost(new_label)

            # check for domination
            dominated = Dict{Label, Nothing}()
            clean_matrix = false
            new_label_is_dominated = false
            for label in buckets[i]
                
                # is the new label dominated?
                if label < new_label    
                    new_label_is_dominated = true
                    break
                end

                # does the new label dominate any label?
                if new_label < label

                    # add dominated to trash
                    Threads.lock(SPINLOCK) do 
                        trash[label] = nothing
                    end

                    # register dominated for later deletion
                    dominated[label] = nothing
                    clean_matrix = true
                end
            end

            # remove the labels dominated by the new label, if necessary
            if clean_matrix
                Threads.lock(SPINLOCK) do 
                    deleteat!(buckets[i], findall(x -> x in keys(dominated), buckets[i]))
                end
            end

            if new_label_is_dominated
                continue
            else # if the new label isn't dominated
                Threads.lock(SPINLOCK) do 
                    push!(buckets[i], new_label)
                    push!(to_extend, new_label)
                end
                was_extended = true
                # println(new_label)
            end 
        end

        # if !was_extended && fast_labelling && 
        #     break
        # end
    end

    # if fast_labelling && !was_extended

    # Multiple labels!
    good_labels = vcat(Vector{Label}[filter((x) -> x.fcost[1] < -1e-4, bucket) for bucket in buckets]...)

    if !isempty(good_labels) && length(good_labels) > len_J/3
        avg_fcost = sum([i.fcost[1] for i in good_labels])/length(good_labels)
        filter!((x) -> x.fcost[1] < avg_fcost, good_labels)
    end

    # Single label!
    min_fcost = Inf
    best_label = nothing
    for label in good_labels
        if label.fcost[1] < min_fcost
            min_fcost = label.fcost[1]
            best_label = label
        end
    end

    if !isempty(good_labels)
        if best_label.weight > W
            error("label $(best_label) is too heavy: $(best_label.weight)")
        end
    
        # println("best new bin:")
        # println("sigma: $(sigma)")
        # println("m: $(best_label.m)")
        # println("k: $(sr_k)")
        # println("rcost: $(best_label.rcost)")
        # println("fcost: $(best_label.fcost[1])")
        
        println("found $(length(good_labels)) bins")
        
        # println("")    
    end


    # return min_fcost, best_label.items
    return min_fcost, Vector{Float64}[label.items for label in good_labels]
end


# # base data
# len_J = 10
# J = [i for i in 1:len_J]

# # reduced costs
# # rc = [i for i in J]
# rc = rand(-10:10, 10)

# positive_rcost = Bool[i > 0 for i in rc]

# w = [i*10 for i in J]
# W = 120

# translated_E = []
# for k in 1:rand(1:len_J*(len_J-1)/2)
#     i = rand(1:len_J-1)
#     j = rand(i+1:len_J)
#     push!(translated_E, [i,j])
# end
# sort!(translated_E)

# println(translated_E)
# println(rc)

# # binarized_E[i][j] = 1 if (i, j) ∈ E
# # binarized_E = Vector{Int64}[Int64[0 for j in 1:len_J] for i in 1:len_J]
# # binarized_E = BitVector[BitVector(undef, len_J) for i in J]
# binarized_E = BitVector[falses(len_J) for i in J]

# for i in J
#     binarized_E[i] .= false
# end

# for (i, j) in translated_E
#     binarized_E[i][j] = true
#     binarized_E[j][i] = true
# end

# min_rcost, new_bin = dp_price(J, len_J, rc, positive_rcost, w, binarized_E, translated_E, W)

# println(min_rcost)
# println(new_bin)




