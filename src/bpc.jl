using JuMP
using Gurobi
using LinearAlgebra
using Combinatorics

include("utils.jl")
include("parallel_pricing.jl")

const GUROBI_ENV = Gurobi.Env()

global LOG_IO = stdout

# struct representing a node in the BB tree
mutable struct Node
    id::Int64
    priority::Int64 
    J::Vector{Int64}
    E::Vector{Vector{Int64}}
    w::Vector{Int64}
    W::Int64
    S::Vector{Vector{Float64}} # lambdas
    mandatory_bags::Vector{Vector{Int64}} # mandatory q ∈ S
    mandatory_bag_amount::Int64
    forbidden_bags::Vector{Vector{Int64}} # forbidden q ∈ S
    item_address::Vector{Int64}
    bounds::Vector{Int64} # [lower_bound, upper_bound]
    solution::Vector{Vector{Int64}}
    bounds_status::Int64 # 0: not optimal, 1: locally optimized, 2: globally optimized 
    subset_row_cuts::Vector{Vector{Int64}} # cut_i for i in cuts | cut_i = [Si_1, Si_2 ... Si_n] 
    branch_history::Vector{Vector{Int64}} # (is_merge, i, j) for each branch
end


"node parameters getter"
get_node_parameters(node::Node) = node.J, node.E, node.w, node.W, node.S

"merges two items i and j such that i < j, merging conflicts, summing their weights"
function merge_items(i::Int64, j::Int64, J::Vector{Int64}, original_w::Vector{Int64}, item_address::Vector{Int64})
    
    # first, make sure i is the lesser value
    # i, j = sort([i,j])

    # println(LOG_IO, "merging $(i) and $(j)")
    # println(LOG_IO, "item_address: $(item_address)")

    # j is the old address
    old_address = j

    # now we move everyone previously at j to i
    for item_i in unmerge_bag_items([j], item_address)
        item_address[item_i] = i
    end
    
    # from the perspective of the current node, there will be one less item
    new_J = J[1:end-1]
    new_w = Int64[0 for j in new_J]

    # println(LOG_IO, "original_w: $(original_w)")
    # println(LOG_IO, "new_w: $(new_w)")
    # println(LOG_IO, "old_address: $(old_address)")
    
    if i == j
        error()
    end

    # "for j in original J"
    for k in eachindex(item_address)

        if item_address[k] == 0 # skip items in mandatory bags
            continue
        end

        # update addresses of items where address > old address of j
        if item_address[k] > old_address
            item_address[k] -= 1
        end

        # add weight to address
        # println(LOG_IO, "k: $(k), $(original_w[k]) at $(item_address[k])")
        new_w[item_address[k]] += original_w[k]
        # println(LOG_IO, "$(new_w)\n$(item_address)")
    end

    return new_J, new_w
end

"makes children with ryan and foster branching"
function make_child_node_with_rf_branch(
    node::Node, 
    j::Int64, 
    q::Vector{Float64},  
    original_w::Vector{Int64}, 
    nodes::Vector{Node}, 
    node_counter::Vector{Int64}, 
    bags_in_use::Vector{Int64}, 
    cuts_binary_data::Vector{BitVector}, 
    cuts_in_use
)
    w = node.w
    W = node.W

    # variable length representation
    # println(LOG_IO, "q: $(q)")
    items_in_q = Int64[i for (i, val) in enumerate(q) if val > 1e-4] 
    # println(LOG_IO, "j: $(j), items_in_q: $(items_in_q)")

    # get items that can be merged with j
    available_to_merge = Int64[i for i in items_in_q if i != j]
    # println(LOG_IO, "available_to_merge: $(available_to_merge)")
            
    # get largest item in bag
    _, i_index = findmax(x -> w[x], available_to_merge)
    i = available_to_merge[i_index]

    # first, make sure i is the lesser value
    i, j = sort([i,j])

    println(LOG_IO, "rf branching on items $(i) and $(j)")

    # pass important bags to child
    J_len = length(node.J)
    new_S = Vector{Float64}[Float64[0.0 for _1 in 1:J_len-1] for _2 in bags_in_use]
    for (new_index, old_index) in enumerate(bags_in_use)
        
        # pass the bags without j
        new_S[new_index] = vcat(node.S[old_index][1:j-1], node.S[old_index][j+1:end])

        # if j was in the old bag, make sure that it is also in the new bag
        if node.S[old_index][j] > .5
            node.S[new_index][i] = 1.0
        end
    end

    # translate cuts considering the merge of i and j
    new_sr_cuts = Vector{Int64}[]
    for n in cuts_in_use
        cut_binary = cuts_binary_data[n]

        if cut_binary[i] && cut_binary[j]

            # update positions and ignore the "new i" (j's new name)
            new_cut = Int64[r > j ? r-1 : r for r in node.subset_row_cuts[n] if r != j]
        elseif cut_binary[j] 
            
            # transform j into i and update positions
            new_cut = sort(Int64[r > j ? r-1 : r == j ? i : r for r in node.subset_row_cuts[n] if r != j])
        else

            # update positions
            new_cut = Int64[r > j ? r-1 : r for r in node.subset_row_cuts[n]]
        end

        if new_cut ∈ new_sr_cuts # does the cut already exists? Possible after merging
            continue
        else # add translated cut
            push!(new_sr_cuts, new_cut)
        end
    end

    # make child
    # pos_child = deepcopy(node)
    # node_counter[1] += 1
    pos_child = Node(
        node_counter[1]+1, # id
        # 1*node.priority,
        length(node.J)+length(node.E)-1,
        deepcopy(node.J),
        deepcopy(node.E),
        deepcopy(node.w),
        deepcopy(node.W),
        deepcopy(new_S), # S
        # Vector{Float64}[],
        deepcopy(node.mandatory_bags),
        deepcopy(node.mandatory_bag_amount),
        deepcopy(node.forbidden_bags),
        deepcopy(node.item_address),
        deepcopy(node.bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
        new_sr_cuts,
        vcat(node.branch_history, Vector{Int64}[[1, i, j]])
    )
    pos_child.bounds[2] = pos_child.mandatory_bag_amount + length(node.J) + 1 # remove prior upper bound

    # update graph
    J, E, w, W, S = get_node_parameters(pos_child)
    pos_child.J, pos_child.w = merge_items(i, j, J, original_w, pos_child.item_address)

    # filter bags that are now too heavy after the merge
    filter!((x) -> sum([pos_child.w[k] for (k, val) in enumerate(x) if val > .5]) < pos_child.W, pos_child.S)


    # filter bags that are now in conflict after the merge
    translated_pos_child_E = translate_edges(pos_child.E, pos_child.item_address)

    binarized_pos_child_E = BitVector[falses(length(pos_child.J)) for _ in pos_child.J]
    for (item_i, item_j) in translated_pos_child_E
        binarized_pos_child_E[item_i][item_j] = true
        binarized_pos_child_E[item_j][item_i] = true
    end

    filter!((x) -> !any(Bool[binarized_pos_child_E[i][k] for (k, val) in enumerate(x) if val > .5]), pos_child.S)
    unique!(pos_child.S)




    # Adding positive child to list
    push!(nodes, pos_child)
    println(LOG_IO, "added node $(pos_child.id) to list")            

    # pass important bags to child, while removing bags that violate the new conflict
    # new_S = Vector{Float64}[deepcopy(node.S[q]) for q in bags_in_use if node.S[q][i] < .5 || node.S[q][j] < .5]
    new_S = unique(Vector{Float64}[deepcopy(node.S[q]) for q in bags_in_use if node.S[q][i] < .5 || node.S[q][j] < .5])
    # new_S = Vector{Float64}[]

    # split branch
    # neg_child = deepcopy(node)
    # node_counter[1] += 1
    neg_child = Node(
        node_counter[1]+2, # id
        # 1*node.priority,
        length(node.J)+length(node.E)+1,
        # length(node.J)+length(node.E) - 100000,
        deepcopy(node.J),
        deepcopy(node.E),
        deepcopy(node.w),
        deepcopy(node.W),
        deepcopy(new_S), # S
        deepcopy(node.mandatory_bags),
        deepcopy(node.mandatory_bag_amount),
        deepcopy(node.forbidden_bags),
        deepcopy(node.item_address),
        deepcopy(node.bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
        deepcopy(Vector{Int64}[node.subset_row_cuts[n] for n in cuts_in_use]),
        vcat(node.branch_history, Vector{Int64}[[0, i, j]])
    )
    neg_child.bounds[2] = neg_child.mandatory_bag_amount + length(node.J) + 1 # remove prior upper bound
    
    # println(LOG_IO, "i: $(i), j: $(j), $(node.item_address)")
    # println(LOG_IO, node.w)

    # E stores conflicts in terms of the original graph
    i = unmerge_bag_items([i], node.item_address)[1]
    j = unmerge_bag_items([j], node.item_address)[1]
    
    # println(LOG_IO, "item_address: $(neg_child.item_address)")
    # println(LOG_IO, "E before: $(neg_child.E)")

    if sort([i, j]) ∈ neg_child.E
        error()
    end

    push!(neg_child.E, sort([i,j]))
    # println(LOG_IO, "E now: $(neg_child.E)")

    # Adding negative child to list
    push!(nodes, neg_child)
    println(LOG_IO, "added node $(neg_child.id) to list")
end

"removes items in q from J, w and E, updating addresses as necessary"
function remove_from_graph(q, q_on_original_G, J, E, original_w, item_address)

    items_amount = length(J)
    amount_to_remove = length(q)

    # println(LOG_IO, "amount_to_remove: $(amount_to_remove)")
    # println(LOG_IO, "q: $(q)")
    # println(LOG_IO, "q_on_original_G: $(q_on_original_G)")
    
    # println(LOG_IO, "item_address: $(item_address)")
    # println(LOG_IO, "w: $(w)")

    # 
    new_J = Int64[j for j in 1:items_amount-amount_to_remove]
    new_w = Int64[0 for j in new_J]

    # println(LOG_IO, "q_on_original_G: $(q_on_original_G)")
    # removing
    for i in q_on_original_G
        item_address[i] = 0
    end

    # updating addresses
    # from largest to smallest k ∈ q:
    #   move to the left all items which address' > k
    for k in amount_to_remove:-1:1
        println(LOG_IO, "removing q[$(k)] = $(q[k])")
        for (j, address) in enumerate(item_address)
            if address > q[k]
                item_address[j] -= 1 
            end
        end
    end
    # println(LOG_IO, "here: $(item_address)")

    # println(LOG_IO, "new item_address: $(item_address)")

    # remove edges containing removed items
    new_E = [e for e in E if !(e[1] ∈ q_on_original_G) && !(e[2] ∈ q_on_original_G)]

    # translate edges

    # update weights
    for (j, address) in enumerate(item_address)
        if address != 0
            new_w[address] += original_w[j]
        end
    end
    # println(LOG_IO, "new_w: $(new_w)")


    return new_J, new_E, new_w
end

"makes children with bag branching and adds them to the list"
function make_child_node_with_bag_branch(node::Node, q::Vector{Float64}, original_w::Vector{Int64}, nodes::Vector{Node}, node_counter::Vector{Int64})
    
    # println(LOG_IO, "q: $(q)")
    
    q = Int64[i for (i, val) in enumerate(q) if val > .5] # variable length representation
    q_on_original_G = unmerge_bag_items(q, node.item_address) # convert q to original G = (V, E), variable length
    
    println(LOG_IO, "branching on bag q: $(q_on_original_G)")
    # println(LOG_IO, "q: $(q)")
    # println(LOG_IO, "q_on_original_G: $(q_on_original_G)")


    # who lives at address j?
    # items_in_address = Vector{Int64}[Int64[] for j in J]
    # for (j, address) in enumerate(item_address)
    #     if address > 0
    #         push!(items_in_address[address], j)
    #     end
    # end

    # get positive child (variable to branch on >= 1) 
    # the items in mandatory bags (λ >= 1) are *removed* from the graph and only considered when computing bounds
    # pos_child = deepcopy(node)
    pos_child = Node(
        node_counter[1]+1, # id
        # 1*node.priority,
        length(node.J)-length(q),
        deepcopy(node.J),
        deepcopy(node.E),
        deepcopy(node.w),
        deepcopy(node.W),
        Vector{Float64}[], # S
        deepcopy(node.mandatory_bags),
        deepcopy(node.mandatory_bag_amount) + 1, # add the new mandatory bag
        deepcopy(node.forbidden_bags),
        deepcopy(node.item_address),
        deepcopy(node.bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
        Vector{Int64}[],
        vcat(node.branch_history, Vector{Int64}[[0, i, j]]),
    )
    pos_child.bounds[2] = pos_child.mandatory_bag_amount + length(node.J)-length(q) + 1 # remove prior upper bound

    J, E, w, W, S = get_node_parameters(pos_child)
    
    # add bag to mandatory set
    push!(pos_child.mandatory_bags, q_on_original_G)

    # remove the items in the mandatory bag from the graph
    pos_child.J, pos_child.E, pos_child.w = remove_from_graph(q, q_on_original_G, J, E, original_w, pos_child.item_address)


    # get negative child (variable to branch on <= 0)
    # the forbidden bag will be cut with a no good cut (non-robust!)
    # neg_child = deepcopy(node)
    neg_child = Node(
        node_counter[1]+2, # id
        # 3*node.priority,
        length(node.J),
        deepcopy(node.J),
        deepcopy(node.E),
        deepcopy(node.w),
        deepcopy(node.W),
        Vector{Float64}[], # S
        deepcopy(node.mandatory_bags),
        deepcopy(node.mandatory_bag_amount),
        deepcopy(node.forbidden_bags),
        deepcopy(node.item_address),
        deepcopy(node.bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
        Vector{Int64}[],
        vcat(node.branch_history, Vector{Int64}[[0, i, j]]),
    )
    neg_child.bounds[2] = node.mandatory_bag_amount + length(node.J) + 1 # remove prior upper bound

    push!(neg_child.forbidden_bags, q_on_original_G)


    # Adding positive child to list
    push!(nodes, pos_child)
    println(LOG_IO, "added node $(pos_child.id) to list")

    # Adding negative child to list
    push!(nodes, neg_child)
    println(LOG_IO, "added node $(neg_child.id) to list")
    
    # count new nodes
    node_counter[1] += 2

    # println(LOG_IO, "bag branching node $(node.id) into $(pos_child.id) and $(neg_child.id) done")
end


"adds new node to queue and node list"
function register_node(node, nodes, queue)

    # clear fields that need clearing
    # node.S = Vector{Float64}[]
    # node.solution = Vector{Int64}[]
    # node.mandatory_bag_amount = length(node.mandatory_bags)
    # node.bounds[2] = node.mandatory_bag_amount + length(node.J) + 1
    # node.interval_graph = false # will need to check again

    println(LOG_IO, "adding node to node list")

    # add new node to list
    push!(nodes, node)
    # node.id = length(nodes)

    println(LOG_IO, "added node $(node.id) to list")
    # println(LOG_IO, node)

    # add to queue at appropriate position
    added = false
    for (i, node_id) in enumerate(queue)
        if node.priority <= nodes[node_id].priority
            insert!(queue, i, node.id)
            added = true
        end
    end
    if !(added)
        push!(queue, node.id)
    end
    
end

"returns simplest lower bound ( ceil(sum(w)/W) )"
function get_simple_lower_bound(w, W)
    return Int64(ceil(sum(w)/W))
end

"L2 bound - Martello1990"
function get_l2_lower_bound(alpha, W, w)
    n1 = Int64[i for i in w if i > W-alpha]
    n2 = Int64[i for i in w if W-alpha >= i && i > W/2]
    n3 = Int64[i for i in w if W/2 >= i && i >= alpha]

    return Int64(length(n1) + length(n2) + max(0, ceil( (sum(n3) - (length(n2)*W - sum(n2) ) )/W  ) )) 
end

"apply First Fit Decreasing heuristic considering conflicts"
function first_fit_decreasing_with_conflicts(J, w, W, E, conflicts; verbose=true)
    bags = Vector{Int64}[Int64[0 for j in J] for i in J]
    # bags = Vector{Bool}[Bool[0 for j in J] for i in J]
    bags_sizes = Int64[0 for i in J]
    bags_slacks = Int64[W for i in J]
    bags_conflicts = Vector{Bool}[Bool[0 for j in J] for i in J]

    J = sort(J, rev=true, by=(x)->w[x])

    bags_amount = 1
    # for each item in decreasing order
    for item_i in J
        verbose && print(LOG_IO, "bagging item $(item_i), weight of $(w[item_i])")
        
        # try lowest-index bag
        for bag_number in 1:length(J)
            if bag_number > bags_amount
                bags_amount = bag_number
            end
        
            item_weight = w[item_i]
            bag_slack = bags_slacks[bag_number]
            bag_size = bags_sizes[bag_number]

            # if item fits and there is no conflict
            if bag_slack - item_weight >= 0 && !(bags_conflicts[bag_number][item_i]) 
                verbose && println(LOG_IO, " at bag $(bag_number)")

                # add item to bag
                bag_size += 1
                # bags[bag_number][bag_size] = item_i
                bags[bag_number][item_i] = 1

                # update bag slack
                bags_slacks[bag_number] -= item_weight
                bags_sizes[bag_number] = bag_size

                # track conflicts
                for conflict_j in conflicts[item_i]
                    verbose && println(LOG_IO, "\tadding conflict with item $(conflict_j)")
                    bags_conflicts[bag_number][conflict_j] = true 
                end
                break
            end
        end
    end



    bags = bags[1:bags_amount]

    return bags, bags_amount
end



"Runs pricing linear programming"
function price_lp(pi_bar, sigma_bar, w, W, J, E, S, forbidden_bags, sr_cuts; verbose=3, epsilon=1e-4)
    # price = Model(Gurobi.Optimizer)
    price = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_silent(price)
    
    # @variable(price, 1 >= x[1:length(J)] >= 0)
    @variable(price, x[1:length(J)], Bin)
    
    @constraint(price, sum([w[j]*x[j] for j ∈ J]) <= W, base_name="capacity")
    for e in E
        @constraint(price, x[e[1]] + x[e[2]] <= 1, base_name="e_($(e[1]), $(e[2]))")
    end
    
    # cut forbidden bags
    for (i, q) in enumerate(forbidden_bags)
        @constraint(price, sum([q[j]*x[j] + (1-q[j])*(1-x[j]) for j in J]) <= length(J) - 1, base_name="forbidden_bag_$(i)")
    end
    

    # Subset row cuts?
    if !isempty(sr_cuts)
        
        # Subset row cuts!
        @variable(price, z[1:length(sr_cuts)] >= 0, Int)
        for (i, cut) in enumerate(sr_cuts)
            @constraint(price, z[i] >= sum([0.5*x[j] for j in cut]) - 1 + 1e-3, base_name="cut_$(i)")
        end

        @objective(price, Min, 1- sum([pi_bar[j]*x[j] for j ∈ J]) - sum([sigma_i*z[i] for (i, sigma_i) in enumerate(sigma_bar)]))

    else
        
        # @objective(price, Min, sum([(1- pi_bar[j])*x[j] for j ∈ J]))
        @objective(price, Min, 1- sum([pi_bar[j]*x[j] for j ∈ J]))
    end

    
    # println(LOG_IO, pi_bar)
    verbose >=3 && println(LOG_IO, price)
    # if !(print_once[1])
    #     println(LOG_IO, price)
    #     print_once[1] = true
    # end
    # println(LOG_IO, price)
    optimize!(price)

    # is the price feasible?
    if termination_status(price) != OPTIMAL
        verbose >= 3 && println(LOG_IO, "price infeasible")
        return 1, nothing
    end
    
    p_obj = objective_value(price)
    verbose >=2 && println(LOG_IO, "̄c = $(p_obj)")

    new_x_bar = value.(price[:x])
    
    # println("best new bin:")
    # println("sigma: $(sigma_bar)")
    # println("m: $(Int64[ length(Int64[j for j in cut if new_x_bar[j] > 0.5]) for cut in sr_cuts ])")
    
    if !isempty(sr_cuts)    
        z_bar = value.(price[:z])
        # println("z_bar: $(z_bar)")
        # println("rcost: $(p_obj + sum(sigma_bar.*z_bar))")    
    else
        # println("rcost: $(p_obj)")
    end
    
    # println("fcost: $(p_obj)")
    
    if !isempty(sr_cuts)    
        for (i, _) in enumerate(sr_cuts)
            println(constraint_by_name(price, "cut_$(i)"))
        end
    end

    # println("")

        
    return p_obj, value.(price[:x])
end

"Searches for subset row cuts by MIP"
function cut_separation(J, lambda_bar, S; verbose=3, epsilon=1e-4)
        
    best_x = Float64[0.0 for j in J]
    n=3
    k=2
    
    cut_separator = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_silent(cut_separator)

    # main constraint
    @variable(cut_separator, x[1:length(J)], Bin)    
    @constraint(cut_separator, sum([x[j] for j ∈ J]) == n, base_name="main_constraint")
    
    # auxiliary constraints
    aux_constraints = ConstraintRef[]
    @variable(cut_separator, aux_w[1:length(lambda_bar)] >= 0, Int)
    for (p, l) in enumerate(lambda_bar)
        new_const = @constraint(cut_separator, sum([S[p][j]*x[j] for j ∈ J])/k + epsilon >= aux_w[p], base_name="aux_constraint_$(p)")
        push!(aux_constraints, new_const)
    end
    
    @objective(cut_separator, Max, sum([aux_w[p]*l  for (p, l) in enumerate(lambda_bar)]) - floor(n/k) )
    
    verbose >=4 && println(LOG_IO, cut_separator)
    optimize!(cut_separator)

    obj = objective_value(cut_separator)
    
    if obj > 0
        best_x = value.(cut_separator[:x])
    end
 
    return obj, k, Int64[i for (i, val) in enumerate(best_x) if val > .5]
end


get_triplets(J, w, W, binarized_E) = filter((x) -> w[x[1]] + w[x[2]] + w[x[3]] < W && !binarized_E[x[1]][x[2]] && !binarized_E[x[1]][x[3]] && !binarized_E[x[2]][x[3]], collect(combinations(J, 3)))

"Searches for subset row cuts by enumeration"
function cut_separation_enum(J, lambda_bar, S, triplets, triplets_tracker, max_per_check; verbose=3, epsilon=1e-4)

    cuts_index = Vector{Int64}[]
    violations = Float64[]
    
    for (i, triplet_i) in enumerate(triplets)

        # cut already added to master
        if triplets_tracker[i]
            continue
        end

        violation = sum(floor(sum(Float64[S[q][j] for j in triplet_i])/2)*l for (q, l) in enumerate(lambda_bar)) - 1

        if violation > epsilon
            push!(violations, violation)
            push!(cuts_index, Int64[i, length(violations)])
            triplets_tracker[i] = true
        end
    end

    if length(cuts_index) > max_per_check
        cuts_index = sort(cuts_index, by=(x) -> -violations[x[2]])[1:max_per_check]
    end

    return violations, cuts_index
end

"""
Column Generation algorithm (CGA)

The CGA searches for new columns (bins) to add based on the dual values of the Master Problem.

# Arguments
- `master`: The master problem model.
- `w`: Vector of item weights.
- `W`: Bin capacity.
- `J`: Set of items.
- `E`: Set of edges representing conflicts between items.
- `lambdas`: Vector of lambda variables representing the packing schemes.
- `S`: Set of packing schemes.
- `S_len`: Length of the set of packing schemes.
- `forbidden_bags`: Set of forbidden packing schemes.
- `subset_row_cuts`: Set of subset row cuts.
- `cuts_binary_data`: Binary representation of the subset row cuts.
- `demand_constraints_ref`: References to the demand constraints in the master problem.
- `cut_constraints_ref`: References to the cut constraints in the master problem.
- `binarized_E`: Binary representation of the conflict graph.
- `verbose`: Verbosity level for logging.
- `max_iter`: Maximum number of iterations for the CGA.
- `epsilon`: Tolerance for numerical stability.
- `using_dp`: Flag to indicate whether to use dynamic programming for pricing.

# Returns
- `m_obj`: Objective value of the master problem.
- `cga_lower_bound`: Lower bound obtained from the CGA.
- `S_len`: Updated length of the set of packing schemes.
"""
function cga(
    master, # master problem
    w, 
    W, 
    J, 
    E, 
    lambdas, 
    S, 
    S_len, 
    forbidden_bags, 
    subset_row_cuts, 
    cuts_binary_data, 
    demand_constraints_ref, 
    cut_constraints_ref, 
    binarized_E; 
    verbose=3, 
    max_iter=10e2, 
    epsilon=1e-4, 
    using_dp=true
)
    check_next_lambda = false

    len_J = length(J)

    m_obj = Inf

    # global print_once = [false]

    # run price, add new columns, check solution, repeat if necessary
    for iteration in 1:max_iter

        optimize!(master)
        # println(LOG_IO, "termination optimal: ", termination_status(master) == OPTIMAL)

        if termination_status(master) != OPTIMAL
            verbose >= 3 && println(LOG_IO, "master infeasible")
            break
        end

        # get values to build price
        m_obj = objective_value(master)
        # println(LOG_IO, "m_obj: $(m_obj)")

        pi_bar = dual.(demand_constraints_ref)
        if isempty(cut_constraints_ref) # trying to get from an empty array will raise error
            sigma_bar = Float64[]
        else
            sigma_bar = dual.(cut_constraints_ref)
        end

        verbose >= 1 && println(LOG_IO, "Z = $(m_obj)")
        

        # run price lp
        if using_dp
            positive_rcost = Bool[i > 0 for i in pi_bar]
            p_obj, new_bins = dp_price(J, len_J, pi_bar, sigma_bar, positive_rcost, w, binarized_E, W, subset_row_cuts, cuts_binary_data, verbose=verbose, epsilon=epsilon)
        else
            p_obj, q = price_lp(pi_bar, sigma_bar, w, W, J, E, S, forbidden_bags, subset_row_cuts, verbose=verbose, epsilon=epsilon)
            new_bins = Vector{Float64}[q]
        end

        # println("last 10 lambdas: $(value.(lambdas)[max(end-10, 1):end])")
        # println("largest lambda: $(max(value.(lambdas)...))")
        # println("p_obj: $(p_obj)")


        if p_obj < -epsilon

            # if using_dp # checking if dp is correct
            #     println(LOG_IO, "p_obj: $(p_obj), adding lambda: $([i for (i, j) in enumerate(q) if j > .5])")
            #     chk_obj, chk_bin = price_function(pi_bar, w, W, J, E, S, forbidden_bags, verbose=0, epsilon=epsilon)
            #     chk_bin = [i for (i, j) in enumerate(chk_bin) if j > .5]
            #     println(LOG_IO, "comparing to mip: $(chk_obj) $(chk_bin)")
            # end

            # pretty_q = Int64[n for (n, v) in enumerate(q) if v > .5]
            # if pretty_q == Int64[1, 5, 27, 29, 43]
            #     print_node_status(node)

            #     println("master: \n")
            #     println(master)

            #     error("here!")
            # end

            for q in new_bins

                # price
                verbose >= 4 && println(LOG_IO, "p_obj: $(p_obj), adding lambda $(S_len+1): $(Int64[n for (n, v) in enumerate(q) if v > .5])")


                if q == S[end]
                    error("repeated bin")
                    # println("repeated bin, stopping cga")
                elseif q ∈ S
                    error("AHA!")
                end

                # # check the value of lambdas added by the cga
                # if check_next_lambda
                #     if value(lambdas[end]) < epsilon
                #         error("AHA!")
                #     end
                # else
                #     check_next_lambda = true
                # end

                # if S[end] == q
                #     println("last bin: $(Int64[k for (k,v) in enumerate(label.items) if v > 0.5])")
                #     error()
                # end

                # add new packing scheme to list
                push!(S, q)

                S_len += 1

                # create new lambda
                push!(lambdas, @variable(master, lower_bound=0, base_name="λ_$(S_len)"))

                # set variable cost on master
                set_objective_function(master, objective_function(master) + lambdas[end])

                # set coefficient of new variable on demand constraints
                for i in J
                    if S[end][i] > 0.5
                        set_normalized_coefficient(demand_constraints_ref[i], lambdas[end], S[end][i])
                    end
                end

                # set coefficient of new variable on cut constraints
                for (i, cut) in enumerate(subset_row_cuts)
                    set_normalized_coefficient(cut_constraints_ref[i], lambdas[end], floor(sum([q[j] for j in cut])/2))
                end

            end

            # println("objective_function: $(objective_function(master))")
            # if !isempty(cut_constraints_ref)
            #     println("cut: $(subset_row_cuts[end])")
            #     println("cut constraint: $(cut_constraints_ref[end])")
            # end
            # println("\n\n")


            # show updated master
            verbose >= 3 && println(LOG_IO, master)

        else
            break
        end
    end

    if m_obj == Inf
        cga_lower_bound = Inf
    else
        cga_lower_bound = Int(ceil(m_obj - epsilon))
    end

    return m_obj, cga_lower_bound, S_len
end

"""Updates the current global bounds and best node if appropriate, and returns the bounds' status:

    0: node not done;
    1: node done;
    2: UB = LB;
"""
function update_bounds_status(node::Node, bounds, best_node, nodes; verbose=1)

    # update global lower bound
    bounds[1], _ = findmin(x -> x.bounds[1], vcat(node, nodes, best_node))


    # default status (continue processing node)
    status = 0
    
    if node.bounds[2] < bounds[2] # is there a new best solution?
        bounds[2] = node.bounds[2]
        best_node[1] = node

        # pretty_solution = get_pretty_solution(translate_solution(node), bounds[2])
        # println(LOG_IO, "node $(node.id): $(bounds[2]) -> $(pretty_solution)")

        if bounds[1] == bounds[2] # is the new solution guaranteed to be globally optimal?
            status = 2
        end
    end

    if node.bounds[1] == node.bounds[2] # is it locally optimal? 
        status = 1

    elseif node.bounds[1] >= bounds[2] # should the algorithm stop processing the node?
        status = 1
    end

    node.bounds_status = status

    # if status ∈ [1,2]
        # pretty_solution = get_pretty_solution(translate_solution(node), bounds[2])
        # println(LOG_IO, "node $(node.id): $(bounds), $(node.mandatory_bag_amount) -> $(pretty_solution)")
    # end
    # return status
end

function solve_bpc(
    J::Vector{Int64}, 
    E::Vector{Vector{Int64}}, 
    w::Vector{Int64}, 
    W::Int64; 
    time_limit::Int64=1800,
    verbose::Int64=1, 
    run_ffd::Bool=true, 
    epsilon::Float64=1e-4,
    max_iter::Int64=100,
    dp::Bool=true,
    )


    item_amount = length(J)

    # [LB, UB]
    bounds = Int64[1, item_amount+1]

    # where can item j be found? (for merging items)
    base_item_adress = Int64[j for j in J]

    # original w
    original_w = deepcopy(w)
    
    # initialize node list
    nodes = Node[Node(
        1, # id
        0, # priority
        J, 
        E, 
        w, 
        W, 
        Vector{Float64}[], # S
        Vector{Int64}[], # mandatory_bags
        0, # mandatory_bag_amount
        Vector{Int64}[], # forbidden_bags
        Int64[j for j in J], # item_address
        deepcopy(bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
        Vector{Int64}[],
        Vector{Int64}[],
    )]
    
    best_node = Node[nodes[1]]
    node_counter = Int64[1]

    not_first_node = false

    node = nodes[1]

    start_time = time()

    is_optimal = true


    # Start the tree
    while !(isempty(nodes))

        if time() - start_time >= time_limit
            println(LOG_IO, "out of time")
            is_optimal = bounds[1] == bounds[2]
            break
        end


        println(LOG_IO, "global bounds: $(bounds)")
        if not_first_node
            pretty_solution = get_pretty_solution(translate_solution(node), bounds[2])
            println(LOG_IO, "node $(node.id): $(bounds), $(node.mandatory_bag_amount) -> $(pretty_solution)")
        end
        println(LOG_IO, "node $(node.id): |J| = $(length(J))")
        # println(LOG_IO, "mandatory_bags: $(node.mandatory_bags)")
        # println(LOG_IO, "solution: $(Vector{Int64}[Int64[i for (i, v) in enumerate(bin) if v > 0.5 ] bin for bin in node.solution])")
        # println(LOG_IO, "solution: $(node.solution)")
        # println(LOG_IO, "item_address: $(node.item_address)")
        println(LOG_IO, "\n")


        # get next node
        _, next_node_position = findmin(x -> x.priority, nodes)
        node = splice!(nodes, next_node_position)

        # for all nodes except the first, check if there is a point in processing it (prune the tree)
        if not_first_node

            filter!(x -> x.bounds[1] < bounds[2], nodes)

            # update bounds status
            update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
            if node.bounds_status != 0 # should we continue processing this node?
                
                if node.bounds_status == 1 # no need to continue *this* node
                    # prune the tree
                    continue
                else # global optimal
                    break
                end
            
            end
        
        else
            not_first_node = true
        end

        J, E, w, W, S = get_node_parameters(node)
        verbose >=1 && println(LOG_IO, "node $(node.id)")


        # println(LOG_IO, "$(J)\n$(w)")
        # println(LOG_IO, "$(node.item_address)")

        # get translated edges
        translated_E = translate_edges(E, node.item_address)
        # println(LOG_IO, "translated_E: $(translated_E)")

        # get translated mandatory/forbidden bags
        forbidden_bags = Vector{Int64}[merge_bag_items(bag, node.item_address, J) for bag in node.forbidden_bags]
        # mandatory_bags = Vector{Int64}[merge_bag_items(bag, node.item_address, J) for bag in node.mandatory_bags]


        ## first try solving the node with heuristics and bounds' properties

        # get initial lower bound ( ⌈∑w/W⌉ ) 
        node.bounds[1] = get_simple_lower_bound(w, W) + node.mandatory_bag_amount
        verbose >= 1 && println(LOG_IO, "⌈∑w/W⌉ lower bound: $(node.bounds[1])")
        
        # naive solution (one item per bag)
        naive_solution = get_naive_solution(J)

        # remove forbidden bags from naive solution
        naive_solution_is_good = remove_forbidden_bags(naive_solution, forbidden_bags)

        if naive_solution_is_good
            
            node.bounds[2] = item_amount + node.mandatory_bag_amount
            node.solution = naive_solution
    
            verbose >= 1 && println(LOG_IO, "Naive upper bound: $(node.bounds[2])")
    
            # update bounds status
            update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
            if node.bounds_status != 0 # is it a global or local optimal?
            
                if node.bounds_status == 1 # no need to continue
                    # prune the tree
                    continue
                else # global optimal
                    break
                end
            
            end

        end

        # [conflicts of i for each i in J]
        conflicts = get_edges(J, translated_E)

        # FFD heuristic for initial solution and upper bound
        if run_ffd
            ffd_solution, ffd_upper_bound = first_fit_decreasing_with_conflicts(J, w, W, translated_E, conflicts, verbose=verbose>1)

            # remove forbidden bags from ffd solution
            ffd_solution_is_good = remove_forbidden_bags(ffd_solution, forbidden_bags)

            if ffd_solution_is_good

                verbose >= 1 && println(LOG_IO, "FFD heuristic upper bound: $(ffd_upper_bound + node.mandatory_bag_amount)")
                
                # if a better solution than one item per bag was found 
                if ffd_upper_bound + node.mandatory_bag_amount < node.bounds[2]
                    node.bounds[2] = ffd_upper_bound + node.mandatory_bag_amount
                    node.solution = ffd_solution
                    
                    # update bounds status
                    update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
                    if node.bounds_status != 0 # is it a global or local optimal?
                        
                        if node.bounds_status == 1 # no need to continue
                            # prune the tree
                            continue
                        else # global optimal
                            break
                        end
                    
                    end   
                
                end

                initial_solution = deepcopy(ffd_solution)
            end

        else
            ffd_solution_is_good = false 
        end
    
        if !(ffd_solution_is_good)
            if naive_solution_is_good
                initial_solution = deepcopy(naive_solution)
            else # no solution
                initial_solution = Vector{Int64}[]     
            end
        end
        
        # try to improve lower bound with martello L2 lower bound
        for alpha in get_not_greater_than_half_capacity(w, W)
            lower_bound = get_l2_lower_bound(alpha, W, w)
            if lower_bound + node.mandatory_bag_amount > node.bounds[1]

                node.bounds[1] = lower_bound + node.mandatory_bag_amount
                verbose >= 1 && println(LOG_IO, "L2 lower bound with α = $(alpha): $(node.bounds[1])")

                # update bounds status
                update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
                if node.bounds_status != 0 # is it a global or local optimal?
                    if node.bounds_status == 1 # no need to continue
                        # prune the tree
                        continue
                    else # global optimal
                        break
                    end
                end    
            end
        end
        

        ## no easy way out, build the LP and try to solve with integer pricing

        # store best solution
        best_solution = deepcopy(initial_solution)
        
        # build master
        # master = Model(Gurobi.Optimizer)
        master = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
        #  set_time_limit_sec(master, 600)
        set_silent(master)
        
        for q in best_solution
            if !(q ∈ S)
                push!(S, q)
            end
        end

        node.S = S
        S_len = length(S)
        # if S_len != length(unique(S))
        #     error("AHA")
        # end
            
        # create lambda variables from existing q ∈ S
        lambdas = VariableRef[]
        for (i, q) in enumerate(S) 
            var = @variable(master, lower_bound=0, base_name="λ_$(i)")
            push!(lambdas, var)
        end
        
        # demand constraints and artificial variables
        artificial_variables = VariableRef[]
        demand_constraints_ref = ConstraintRef[]
        for i in J            
            av = @variable(master, lower_bound=0, base_name="av_$(i)")
            con_ref = @constraint(master, sum([sum(S[q][i]*lambdas[q]) for q in 1:S_len]) + av >= 1, base_name="demand_$(i)")
            push!(demand_constraints_ref, con_ref)  
            push!(artificial_variables, av)
        end

        # subset_row_cuts (Jepsen, 2008)
        cut_artificial_variables = VariableRef[]
        cut_constraints_ref = ConstraintRef[]
        for (n, cut_n) in enumerate(node.subset_row_cuts)
            con_ref = @constraint(master, sum([floor(sum([S[p][i] for i in cut_n])/2)*l_p for (p, l_p) in enumerate(lambdas)]) <= 1, base_name="sr_cut_$(n)")
            push!(cut_constraints_ref, con_ref)
        end
    
        # objective function
        @objective(master, Min, sum(lambdas) + 10000*item_amount*sum(artificial_variables))

        # show initial master
        verbose >= 2 && println(LOG_IO, master)
    
        ## BCPA

        # apply cga
        
        # cuts auxiliary data for processing
        J_len = length(J)
        cuts_binary_data = BitVector[falses(J_len) for i in node.subset_row_cuts]
        for (i, cut) in enumerate(node.subset_row_cuts)
            for j in cut
                cuts_binary_data[i][j] = true
            end
        end

        # binarized E, for fast conflict checks
        binarized_E = BitVector[falses(J_len) for i in J]
        for (i, j) in translated_E
            binarized_E[i][j] = true
            binarized_E[j][i] = true
        end

        # cuts control
        max_cuts = 10*length(J)
        max_checks_per_node = 10
        max_cuts_per_check = 10
        checks_done_this_node = 0

        triplets = get_triplets(J, w, W, binarized_E)
        triplets_tracker = falses(length(triplets))

        
        lambda_bar = Float64[]
        z, cga_lb = Inf, Inf
        cga_lb_break = false
        continue_adding_cuts = true
        while continue_adding_cuts # cga and cut adding loop
        # for i in 1:max_checks_per_node # cga and cut adding loop
            
            z, cga_lb, S_len = cga(master, w, W, J, translated_E, lambdas, node.S, S_len, forbidden_bags, node.subset_row_cuts, cuts_binary_data, demand_constraints_ref, cut_constraints_ref, binarized_E, verbose=verbose, epsilon=epsilon, max_iter=max_iter, using_dp=true)
            if termination_status(master) != OPTIMAL
                println(LOG_IO, "node $(node.id) linear programming failed to optimize")
                break
            end
    
            if cga_lb + node.mandatory_bag_amount > node.bounds[1]
    
                verbose >= 1 && println(LOG_IO, "CGA lower bound: $(cga_lb + node.mandatory_bag_amount)")
    
                node.bounds[1] = cga_lb + node.mandatory_bag_amount
    
                # update bounds status
                update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
                if node.bounds_status != 0 # is it a global or local optimal?
                    if node.bounds_status == 1 # no need to continue
                        # prune the tree
                        cga_lb_break = true
                        break
                    end
                end  
            end

            if time() - start_time >= time_limit
                println(LOG_IO, "out of time")
                is_optimal = bounds[1] == bounds[2]
                break
            end                

            # if too many cuts already, skip the stop the cut adding loop
            if checks_done_this_node >= max_checks_per_node || length(node.subset_row_cuts) >= max_cuts
                continue_adding_cuts = false
                break
            end

            # get lambda values of the solution
            lambda_bar = value.(lambdas)

            # try to find subset row cuts
            # violation, k, cut_data = cut_separation(J, lambda_bar, S)
            violations, new_cuts_index = cut_separation_enum(J, lambda_bar, S, triplets, triplets_tracker, max_cuts_per_check)
            if isempty(new_cuts_index) # no cuts to add
                break
            else
                checks_done_this_node += 1
                for (t_index, v_index) in new_cuts_index

                    cut_data = deepcopy(triplets[t_index])
                    # println(LOG_IO, "adding cut with violation = $(violations[v_index]): $(cut_data)")
    
                    # add cut data to node
                    push!(node.subset_row_cuts, cut_data)
                    n = length(node.subset_row_cuts)
                    
                    # add cut to master
                    con_ref = @constraint(master, sum([floor(sum([S[p][i] for i in cut_data])/2)*l_p for (p, l_p) in enumerate(lambdas)]) <= 1, base_name="sr_cut_$(n)")
                    # println(con_ref)
    
                    # add to constraint to reference array
                    push!(cut_constraints_ref, con_ref)
    
                    # update auxiliary cut data
                    push!(cuts_binary_data, falses(J_len))
                    for r in cut_data
                        cuts_binary_data[end][r] = true
                    end
    
                    # println(master)
                end
            end
        end

        if cga_lb_break
            continue # prune the tree
        end

        # get lambda values of the solution
        lambda_bar = value.(lambdas)
        
        # get branching candidates
        bags_in_use = get_bags_in_use(lambda_bar, S, S_len, J; epsilon=epsilon)
        # most_fractional_bag = make_branching_analysis(bags_in_use, lambda_bar, S, S_len, conflicts, J, w, epsilon=1e-4)
        most_fractional_bag, most_fractional_item = find_ryan_foster_branch(bags_in_use, lambda_bar, S, w, epsilon=epsilon)

        # println(LOG_IO, "lambda_bar: $(lambda_bar)")
        # println(LOG_IO, "bags_in_use: $(bags_in_use)")
        # println(LOG_IO, "most_fractional_bag: $(most_fractional_bag)")
        # println(LOG_IO, "most_fractional_item: $(most_fractional_item)")

        # is there an integer bag to branch on? 
        # That is, 0 < λ_i < 1 | j ∈ {0,1} ∀ j ∈ λ_i
        # if most_fractional_bag != -1

        #     # get q to branch on
        #     q = S[most_fractional_bag]

        #     make_child_node_with_bag_branch(node, q, original_w, nodes, node_counter)

        # (Ryan and Foster branching)
        if most_fractional_item[1] != -1

            q = S[most_fractional_bag]
            j = most_fractional_item

            println(LOG_IO, "enum_q: $([i for i in enumerate(q)])")
            println(LOG_IO, "j: $(j)")

            # get cuts that are actually in use
            if isempty(cut_constraints_ref) # trying to get from an empty array will raise error
                sigma_bar = Float64[]
            else
                sigma_bar = dual.(cut_constraints_ref)
            end
            cuts_in_use = Int64[k for (k,v) in enumerate(sigma_bar) if v < -epsilon]
            # println("sigma_bar = $(sigma_bar)")
            # println("cuts_in_use = $(cuts_in_use)")

            make_child_node_with_rf_branch(node, j, q, original_w, nodes, node_counter, bags_in_use, cuts_binary_data, cuts_in_use)
            node_counter[1] += 2

        else # the solution is integer!

            verbose >= 1 && println(LOG_IO, "node $(node.id) finished with integer solution by CGA")

            # get solution values
            lambda_bar = value.(lambdas)
            x_bar, cga_ub = get_x(lambda_bar, S, S_len, J, epsilon=epsilon)
            
            
            # println("x_bar_unbin = $(Vector{Int64}[Int64[k for (k,v) in enumerate(bin) if v > 0.5] for bin in x_bar])")
            # println("x_bar_values = $(Vector{Vector{Float64}}[Vector{Float64}[Float64[k,v] for (k,v) in enumerate(bin) if v > 0] for bin in x_bar])")
        
            # treat current solution
            current_solution = round_up_solution(x_bar)
            current_solution, cga_ub = prune_excess_with_priority(current_solution, J, w, epsilon=epsilon)

            # was there an improvement?
            if cga_ub + node.mandatory_bag_amount < node.bounds[2]
        
                node.bounds[2] = cga_ub + node.mandatory_bag_amount
                best_solution = deepcopy(current_solution)
                node.solution = best_solution
                
                # update bounds status
                update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
                if node.bounds_status != 0 # should we stop processing?
                    if node.bounds_status == 1 # no need to continue
                        # prune the tree
                        continue
                    else # global optimal
                        break
                        # return best_solution, bounds[2]
                    end
                end  
            end
        end
    end

    node = best_node[1]

    print_node_status(node, original_w)

    verbose >= 1 && println(LOG_IO, "tree finished")
    
    solution = translate_solution(node, epsilon=epsilon)
    final_solution = get_pretty_solution(solution, bounds[2])

    println(LOG_IO, "bounds: $(bounds)")
    println(LOG_IO, "node.bounds: $(node.bounds)")
    println(LOG_IO, "node.bounds_status: $(node.bounds_status)")

    return final_solution, bounds[2], is_optimal
end

