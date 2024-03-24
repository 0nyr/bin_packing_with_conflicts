using JuMP
using Gurobi
using LinearAlgebra

include("utils.jl")

const GUROBI_ENV = Gurobi.Env()
    

# struct representing a node in the BB tree
mutable struct Node
    id::Int64
    priority::Int64 
    J::Vector{Int64}
    E::Vector{Vector{Int64}}
    w::Vector{Int64}
    W::Int64
    S::Vector{Vector{Float32}} # needs to be pruned after a Ryan-Foster merge!
    mandatory_bags::Vector{Vector{Int64}} # mandatory q ∈ S
    mandatory_bag_amount::Int64
    forbidden_bags::Vector{Vector{Int64}} # forbidden q ∈ S
    item_address::Vector{Int64}
    interval_graph::Bool # can DP-flow be used? 
    bounds::Vector{Int64} # [lower_bound, upper_bound]
    solution::Vector{Vector{Int64}}
    bounds_status::Int64 # 0: not optimal, 1: locally optimized, 2: globally optimized 
end


"node parameters getter"
get_node_parameters(node::Node) = node.J, node.E, node.w, node.W, node.S

"merges two items i and j, merging conflicts, summing their weights"
function merge_items(i::Int64, j::Int64, J::Vector{Int64}, original_w::Vector{Int64}, item_address::Vector{Int64})
    
    # first, make sure i is the lesser value
    i, j = sort([i,j])

    # println("merging $(i) and $(j)")
    # println("item_address: $(item_address)")

    # j is the old address
    old_address = j

    # now we move everyone previously at j to i
    for item_i in unmerge_bag_items([j], item_address)
        item_address[item_i] = i
    end
    
    # from the perspective of the current node, there will be one less item
    new_J = J[1:end-1]
    new_w = Int64[0 for j in new_J]

    # println("original_w: $(original_w)")
    # println("new_w: $(new_w)")
    # println("old_address: $(old_address)")
    
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
        # println("k: $(k), $(original_w[k]) at $(item_address[k])")
        new_w[item_address[k]] += original_w[k]
        # println("$(new_w)\n$(item_address)")
    end

    return new_J, new_w
end

"makes children with ryan and foster branching"
function make_child_node_with_rf_branch(node::Node, j::Int64, q::Vector{Float32},  original_w::Vector{Int64}, nodes::Vector{Node}, node_counter::Vector{Int64})
    
    w = node.w
    W = node.W

    # variable length representation
    # println("q: $(q)")
    items_in_q = Int64[i for (i, val) in enumerate(q) if val > 1e-4] 
    # println("j: $(j), items_in_q: $(items_in_q)")

    # get items that can be merged with j
    available_to_merge = Int64[i for i in items_in_q if i != j && w[i] + w[j] < W]
    # println("available_to_merge: $(available_to_merge)")
    
    
    # is there an item such that merging with j is feasible? (w_i + w_j < W) 
    if !(isempty(available_to_merge))
        
        # get largest item in bag, except the fractional item
        _, i_index = findmax(x -> w[x], available_to_merge)
        i = available_to_merge[i_index]

        # make child
        # pos_child = deepcopy(node)
        node_counter[1] += 1
        pos_child = Node(
            node_counter[1], # id
            # 1*node.priority,
            length(node.J)-1,
            deepcopy(node.J),
            deepcopy(node.E),
            deepcopy(node.w),
            deepcopy(node.W),
            Vector{Float32}[], # S
            deepcopy(node.mandatory_bags),
            deepcopy(node.mandatory_bag_amount),
            deepcopy(node.forbidden_bags),
            deepcopy(node.item_address),
            false, # interval_graph
            deepcopy(node.bounds), # node bounds
            Vector{Int64}[], # solution
            0, # bounds_status
        )
        pos_child.bounds[2] = pos_child.mandatory_bag_amount + length(node.J) + 1 # remove prior upper bound

        # update graph
        J, E, w, W, S = get_node_parameters(pos_child)
        pos_child.J, pos_child.w = merge_items(i, j, J, original_w, pos_child.item_address)
    
        # Adding positive child to list
        push!(nodes, pos_child)
        # println("added node $(pos_child.id) to list")            
        
    else # merging is infeasible, only create negative child, with heaviest i (higher impact on solution) 

        # get the heaviest i | i != j
        i_weight = -Inf
        heaviest_i = -1
        # println("items_in_q: $(items_in_q)")
        for i in items_in_q
            if node.w[i] > i_weight && i != j
                i_weight = node.w[i]
                heaviest_i = i 
            end
        end
        i = heaviest_i

    end

    println("rf branching on items $(i) and $(j)")

    # split branch
    # neg_child = deepcopy(node)
    node_counter[1] += 1
    neg_child = Node(
        node_counter[1], # id
        # 1*node.priority,
        length(node.J),
        deepcopy(node.J),
        deepcopy(node.E),
        deepcopy(node.w),
        deepcopy(node.W),
        Vector{Float32}[], # S
        deepcopy(node.mandatory_bags),
        deepcopy(node.mandatory_bag_amount),
        deepcopy(node.forbidden_bags),
        deepcopy(node.item_address),
        false, # interval_graph
        deepcopy(node.bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
    )
    neg_child.bounds[2] = neg_child.mandatory_bag_amount + length(node.J) + 1 # remove prior upper bound
    
    println("i: $(i), j: $(j), $(node.item_address)")
    println(node.w)

    # E stores conflicts in terms of the original graph
    i = unmerge_bag_items([i], node.item_address)[1]
    j = unmerge_bag_items([j], node.item_address)[1]
    
    println("item_address: $(neg_child.item_address)")
    println("E before: $(neg_child.E)")

    if sort([i, j]) ∈ neg_child.E
        error()
    end

    push!(neg_child.E, sort([i,j]))
    println("E now: $(neg_child.E)")

    # Adding negative child to list
    push!(nodes, neg_child)
    # println("added node $(neg_child.id) to list")
end

"removes items in q from J, w and E, updating addresses as necessary"
function remove_from_graph(q, q_on_original_G, J, E, original_w, item_address)

    items_amount = length(J)
    amount_to_remove = length(q)

    # println("amount_to_remove: $(amount_to_remove)")
    # println("q: $(q)")
    # println("q_on_original_G: $(q_on_original_G)")
    
    # println("item_address: $(item_address)")
    # println("w: $(w)")

    # 
    new_J = Int64[j for j in 1:items_amount-amount_to_remove]
    new_w = Int64[0 for j in new_J]

    # println("q_on_original_G: $(q_on_original_G)")
    # removing
    for i in q_on_original_G
        item_address[i] = 0
    end

    # updating addresses
    # from largest to smallest k ∈ q:
    #   move to the left all items which address' > k
    for k in amount_to_remove:-1:1
        println("removing q[$(k)] = $(q[k])")
        for (j, address) in enumerate(item_address)
            if address > q[k]
                item_address[j] -= 1 
            end
        end
    end
    # println("here: $(item_address)")

    # println("new item_address: $(item_address)")

    # remove edges containing removed items
    new_E = [e for e in E if !(e[1] ∈ q_on_original_G) && !(e[2] ∈ q_on_original_G)]

    # translate edges

    # update weights
    for (j, address) in enumerate(item_address)
        if address != 0
            new_w[address] += original_w[j]
        end
    end
    # println("new_w: $(new_w)")


    return new_J, new_E, new_w
end

"makes children with bag branching and adds them to the list"
function make_child_node_with_bag_branch(node::Node, q::Vector{Float32}, original_w::Vector{Int64}, nodes::Vector{Node}, node_counter::Vector{Int64})
    
    # println("q: $(q)")
    
    q = Int64[i for (i, val) in enumerate(q) if val > .5] # variable length representation
    q_on_original_G = unmerge_bag_items(q, node.item_address) # convert q to original G = (V, E), variable length
    
    println("branching on bag q: $(q_on_original_G)")
    # println("q: $(q)")
    # println("q_on_original_G: $(q_on_original_G)")


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
        Vector{Float32}[], # S
        deepcopy(node.mandatory_bags),
        deepcopy(node.mandatory_bag_amount) + 1, # add the new mandatory bag
        deepcopy(node.forbidden_bags),
        deepcopy(node.item_address),
        false, # interval_graph
        deepcopy(node.bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
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
        Vector{Float32}[], # S
        deepcopy(node.mandatory_bags),
        deepcopy(node.mandatory_bag_amount),
        deepcopy(node.forbidden_bags),
        deepcopy(node.item_address),
        false, # interval_graph
        deepcopy(node.bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
    )
    neg_child.bounds[2] = node.mandatory_bag_amount + length(node.J) + 1 # remove prior upper bound

    push!(neg_child.forbidden_bags, q_on_original_G)


    # Adding positive child to list
    push!(nodes, pos_child)
    # println("added node $(pos_child.id) to list")

    # Adding negative child to list
    push!(nodes, neg_child)
    # println("added node $(neg_child.id) to list")
    
    # count new nodes
    node_counter[1] += 2

    # println("bag branching node $(node.id) into $(pos_child.id) and $(neg_child.id) done")
end

"adds new node to queue and node list"
function register_node(node, nodes, queue)

    # clear fields that need clearing
    # node.S = Vector{Float32}[]
    # node.solution = Vector{Int64}[]
    # node.mandatory_bag_amount = length(node.mandatory_bags)
    # node.bounds[2] = node.mandatory_bag_amount + length(node.J) + 1
    # node.interval_graph = false # will need to check again

    println("adding node to node list")

    # add new node to list
    push!(nodes, node)
    # node.id = length(nodes)

    println("added node $(node.id) to list")
    # println(node)

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
        verbose && print("bagging item $(item_i), weight of $(w[item_i])")
        
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
                verbose && println(" at bag $(bag_number)")

                # add item to bag
                bag_size += 1
                # bags[bag_number][bag_size] = item_i
                bags[bag_number][item_i] = 1

                # update bag slack
                bags_slacks[bag_number] -= item_weight
                bags_sizes[bag_number] = bag_size

                # track conflicts
                for conflict_j in conflicts[item_i]
                    verbose && println("\tadding conflict with item $(conflict_j)")
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
function price_lp(pi_bar, w, W, J, E, S, forbidden_bags; verbose=3, epsilon=1e-4)
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
    
    # @objective(price, Min, sum([(1- pi_bar[j])*x[j] for j ∈ J]))
    @objective(price, Min, 1- sum([pi_bar[j]*x[j] for j ∈ J]))
    # set_silent(price)

    # println(pi_bar)
    verbose >=3 && println(price)
    # if !(print_once[1])
    #     println(price)
    #     print_once[1] = true
    # end
    # println(price)
    optimize!(price)

    # is the price feasible?
    if termination_status(price) != OPTIMAL
        verbose >= 3 && println("price infeasible")
        return 1, nothing
    end
    
    p_obj = objective_value(price)
    verbose >=2 && println("̄c = $(p_obj)")
        
    return p_obj, value.(price[:x])
end

"Runs relaxed pricing linear programming, but constrains the new lambda to be integer"
function rounded_relaxed_price_lp(pi_bar, w, W, J, E, S, forbidden_bags; verbose=3, epsilon=1e-4)
    # price = Model(Gurobi.Optimizer)
    price = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_silent(price)
    @variable(price, 1 >= x[1:length(J)] >= 0)
    # @variable(price, x[1:length(J)], Bin)
    @constraint(price, sum([w[j]*x[j] for j ∈ J]) <= W, base_name="capacity")
    for e in E
        @constraint(price, x[e[1]] + x[e[2]] <= 1, base_name="e_($(e[1]), $(e[2]))")
    end

    # cut forbidden bags
    for (i, q) in enumerate(forbidden_bags)
        @constraint(price, sum([q[j]*x[j] + (1-q[j])*(1-x[j]) for j in J]) <= length(J) - 1, base_name="forbidden_bag_$(i)")
    end
    
    # @objective(price, Min, sum([(1- pi_bar[j])*x[j] for j ∈ J]))
    @objective(price, Min, 1- sum([pi_bar[j]*x[j] for j ∈ J]))
    # set_silent(price)

    # println(pi_bar)
    verbose >=3 && println(price)
    optimize!(price)

    # is the price feasible?
    if termination_status(price) != OPTIMAL
        verbose >= 3 && println("price infeasible")
        return 1, nothing
    end
    
    p_obj = objective_value(price)
    verbose >=2 && println("̄c = $(p_obj)")
        
    # remove fractional part from new q
    q = floor_vector(value.(price[:x]), epsilon=epsilon)
    
    # after removing the fractional item, is the new q still useful?
    if !(sum(q) > 2 - epsilon) || q ∈ S 
        verbose >= 1 && println("integer pricing empty")
        return 1, q
    end

    return p_obj, q
end

function cga(master, price_function, w, W, J, E, lambdas, S, S_len, forbidden_bags; verbose=3, max_iter=10e2, epsilon=1e-4)
    
    m_obj = Inf

    # global print_once = [false]

    # run price, add new columns, check solution, repeat if necessary
    iteration = 1
    for iteration in 1:max_iter

        optimize!(master)
        # println("termination optimal: ", termination_status(master) == OPTIMAL)

        if termination_status(master) != OPTIMAL
            verbose >= 3 && println("master infeasible")
            break
        end


        # get values to build price
        m_obj, demand_constraints, pi_bar = get_master_data_for_pricing(master, J, verbose=verbose)
        
        # run price lp
        p_obj, q = price_function(pi_bar, w, W, J, E, S, forbidden_bags, verbose=verbose, epsilon=epsilon)

        if p_obj < -epsilon

            # price
            verbose >= 3 && println("adding lambda: $(q)")

            # add new packing scheme to list
            push!(S, q)
            
            S_len += 1

            # create new lambda
            push!(lambdas, @variable(master, lower_bound=0, base_name="λ_$(S_len)"))

            # set variable cost on master
            set_objective_function(master, objective_function(master) + lambdas[end])

            # set coefficient of new variable on demand constraints
            for i in J
                if S[end][i] > 0
                    set_normalized_coefficient(demand_constraints[i], lambdas[end], S[end][i])
                end
            end

            # show updated master
            verbose >= 3 && println(master)

        else
            break
        end
    end

    if iteration == max_iter && verbose >= 3
        println("CGA reached max iterations before exiting (check price objective value)")
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
        # println("node $(node.id): $(bounds[2]) -> $(pretty_solution)")

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
        # println("node $(node.id): $(bounds), $(node.mandatory_bag_amount) -> $(pretty_solution)")
    # end
    # return status
end

function solve_bpc(
    J::Vector{Int64}, 
    E::Vector{Vector{Int64}}, 
    w::Vector{Int64}, 
    W::Int64; 
    verbose::Int64=1, 
    run_ffd::Bool=true, 
    epsilon::Float64=1e-4,
    max_iter::Int64=100,
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
        Vector{Float32}[], # S
        Vector{Int64}[], # mandatory_bags
        0, # mandatory_bag_amount
        Vector{Int64}[], # forbidden_bags
        Int64[j for j in J], # item_address
        false, # interval_graph
        deepcopy(bounds), # node bounds
        Vector{Int64}[], # solution
        0, # bounds_status
    )]
    
    best_node = Node[nodes[1]]
    node_counter = Int64[1]

    not_first_node = false

    node = nodes[1]

    # Start the tree
    while !(isempty(nodes))

        pretty_solution = get_pretty_solution(translate_solution(node), bounds[2])
        println("node $(node.id): $(bounds), $(node.mandatory_bag_amount) -> $(pretty_solution)")
        println("node $(node.id): |J| = $(length(J))")
        println("mandatory_bags: $(node.mandatory_bags)")
        println("solution: $(node.solution)")
        println("item_address: $(node.item_address)")
        
        println("\n")


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
        verbose >=1 && println("node $(node.id)")


        # println("$(J)\n$(w)")
        # println("$(node.item_address)")

        # get translated edges
        translated_E = translate_edges(E, node.item_address)
        # println("translated_E: $(translated_E)")

        # get translated mandatory/forbidden bags
        forbidden_bags = Vector{Int64}[merge_bag_items(bag, node.item_address, J) for bag in node.forbidden_bags]
        # mandatory_bags = Vector{Int64}[merge_bag_items(bag, node.item_address, J) for bag in node.mandatory_bags]


        ## first try solving the node with heuristics and bounds' properties

        # get initial lower bound ( ⌈∑w/W⌉ ) 
        node.bounds[1] = get_simple_lower_bound(w, W) + node.mandatory_bag_amount
        verbose >= 1 && println("⌈∑w/W⌉ lower bound: $(node.bounds[1])")
        
        # naive solution (one item per bag)
        naive_solution = get_naive_solution(J)

        # remove forbidden bags from naive solution
        naive_solution_is_good = remove_forbidden_bags(naive_solution, forbidden_bags)

        if naive_solution_is_good
            
            node.bounds[2] = item_amount + node.mandatory_bag_amount
            node.solution = naive_solution
    
            verbose >= 1 && println("Naive upper bound: $(node.bounds[2])")
    
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

                verbose >= 1 && println("FFD heuristic upper bound: $(ffd_upper_bound + node.mandatory_bag_amount)")
                
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
        end
    
        if run_ffd && !(ffd_solution_is_good) && naive_solution_is_good
            initial_solution = deepcopy(naive_solution)
        else # no solution
            initial_solution = Vector{Int64}[]
        end
        
        # try to improve lower bound with martello L2 lower bound
        for alpha in get_not_greater_than_half_capacity(w, W)
            lower_bound = get_l2_lower_bound(alpha, W, w)
            if lower_bound + node.mandatory_bag_amount > node.bounds[1]

                node.bounds[1] = lower_bound + node.mandatory_bag_amount
                verbose >= 1 && println("L2 lower bound with α = $(alpha): $(node.bounds[1])")

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
        set_silent(master)
        
        # add the naive solution as lambda variables (can serve as artificial variables)
        S = Vector{Float32}[q for q in naive_solution]

        # if FFD was ran pass the relevant bags to S
        if run_ffd
            for q in ffd_solution
                if sum(q) > 1 # pass the relevant bags to S
                    push!(S, q)
                end
            end
        end
        node.S = S
        S_len = length(S)
    
        # create lambda variables from existing q ∈ S
        lambdas = VariableRef[]
        for (i, q) in enumerate(S) 
            var = @variable(master, lower_bound=0, base_name="λ_$(i)")
            push!(lambdas, var)
        end
        
        # demand constraints and artificial variables
        artificial_variables = VariableRef[]
        for i in J
            au = @variable(master, lower_bound=0, base_name="a_u_$(i)")
            al = @variable(master, lower_bound=0, base_name="a_l_$(i)")

            @constraint(master, sum([sum(S[q][i]*lambdas[q]) for q in 1:S_len]) + au - al == 1, base_name="demand_$(i)")

            push!(artificial_variables, au, al)
        end
    
        # objective function
        @objective(master, Min, sum(lambdas) + 1000*item_amount*sum(artificial_variables))
    
        # show initial master
        verbose >= 2 && println(master)
    
            
        # run column generation with specialized pricing
        # if node.interval_graph...

        # run column generation with integer pricing
        m_obj, cga_ub, S_len = cga(master, rounded_relaxed_price_lp, w, W, J, translated_E, lambdas, S, S_len, forbidden_bags, verbose=verbose, epsilon=epsilon, max_iter=max_iter)
        if termination_status(master) == OPTIMAL
            
            # get solution values
            lambda_bar = value.(lambdas)
            x_bar, cga_ub = get_x(lambda_bar, S, S_len, J, epsilon=epsilon)
        
            # treat current solution
            current_solution = round_up_solution(x_bar)
            current_solution, cga_ub = prune_excess_with_priority(current_solution, J, w, epsilon=epsilon)
            
            verbose >= 1 && println("Integer CGA upper bound: $(cga_ub) + $(node.mandatory_bag_amount) mandatory bags")
        
            # was there an improvement from the heuristic?
            if cga_ub + node.mandatory_bag_amount < node.bounds[2]
        
                node.bounds[2] = cga_ub + node.mandatory_bag_amount
                best_solution = deepcopy(current_solution)
                node.solution = best_solution
                
                # update bounds status
                update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
                if node.bounds_status != 0 # is it a global or local optimal?
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



        ## BCPA

        # apply cga
        z, cga_lb, S_len = cga(master, price_lp, w, W, J, translated_E, lambdas, node.S, S_len, forbidden_bags, verbose=verbose, epsilon=epsilon, max_iter=max_iter)
        if termination_status(master) != OPTIMAL
            println("node $(node.id) linear programming failed to optimize")
            break
        end

        # # is there already a better or equal solution?
        # if cga_lb + node.mandatory_bag_amount >= bounds[2]
        #     continue # close node
        # end

        if cga_lb + node.mandatory_bag_amount > node.bounds[1]

            verbose >= 1 && println("CGA lower bound: $(cga_lb) + $(node.mandatory_bag_amount) mandatory bags")

            node.bounds[1] = cga_lb + node.mandatory_bag_amount

            # update bounds status
            update_bounds_status(node, bounds, best_node, nodes, verbose=verbose)
            if node.bounds_status != 0 # is it a global or local optimal?
                if node.bounds_status == 1 # no need to continue
                    # prune the tree
                    continue
                end
            end  
        end


        # get lambda values of the solution
        lambda_bar = value.(lambdas)
        
        # get branching candidates
        bags_in_use, lambdas_in_use = get_bags_in_use(lambda_bar, S, S_len, J; epsilon=epsilon)
        most_fractional_bag = make_branching_analysis(bags_in_use, lambdas_in_use, lambda_bar, S, S_len, conflicts, J, w, epsilon=1e-4)

        # println("lambda_bar: $(lambda_bar)")
        # println("bags_in_use: $(bags_in_use)")
        # println("most_fractional_bag: $(most_fractional_bag)")
        # println("most_fractional_item: $(most_fractional_item)")

        # is there an integer bag to branch on? 
        # That is, 0 < λ_i < 1 | j ∈ {0,1} ∀ j ∈ λ_i
        if most_fractional_bag != -1

            # get q to branch on
            q = S[most_fractional_bag]

            make_child_node_with_bag_branch(node, q, original_w, nodes, node_counter)

        # if not, is there an item to branch on? (Ryan and Foster branching)
        # elseif most_fractional_item[1] != -1

        #     q = S[most_fractional_item[1]]
        #     j = most_fractional_item[2]

        #     println("q: $(q)")

        #     make_child_node_with_rf_branch(node, j, q, original_w, nodes, node_counter)

        else # the solution is integer!

            # get solution values
            lambda_bar = value.(lambdas)
            x_bar, cga_ub = get_x(lambda_bar, S, S_len, J, epsilon=epsilon)
        
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

    verbose >= 1 && println("tree finished")
    
    solution = translate_solution(best_node[1], epsilon=epsilon)
    final_solution = get_pretty_solution(solution, bounds[2])

    println("bounds: $(bounds)")
    println("node.bounds: $(node.bounds)")
    println("node.bounds_status: $(node.bounds_status)")

    return final_solution, bounds[2]
end