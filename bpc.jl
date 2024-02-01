using JuMP
using GLPK
using LinearAlgebra

include("utils.jl")

# struct representing a node in the BB tree
mutable struct Node
    id::Int64
    J::Array{Int64}
    E::Array{Array{Int64}}
    w::Array{Int64}
    W::Int64
    S::Array{Float32} # needs to be pruned after a Ryan-Foster merge!
    branches::Array{Int64} # q ∈ S: -1 if no bound, 0 if q is cut, 1 if mandatory
    # lambdas::Array{VariableRef}
    item_address::Array{Array{Int64}}
    is_interval_graph::Bool # can DP-flow be used? 
    bounds::Array{Int64} # [lower_bound, upper_bound]
    bounds_status::Int64
    solution::Array{Array{Int64}}
end

"node parameters getter"
get_node_parameters(node:Node) = node.J, node.E, node.w, node.W, node.S, node.branches

"merges two items i and j, merging conflicts, summing their weights"
function merge_items(i, j, J, w, item_address)
    
    # first, make sure i is the lesser value
    i, j = sort([i,j])

    old_address = item_address[j]

    new_J = J[1:end-1]
    new_w = Int64[0 for j in J]

    # update address of j
    item_address[j] = item_address[i]
    
    # "for j in original J"
    for k in eachindex(item_address)

        # update addresses of items where address > old address of j
        if item_address[k] > old_address
            item_address[k] -= 1
        end

        # add weight to address
        new_w[item_address[k]] += w[k]
    end

    return new_J, new_w
end

"makes children with ryan and foster branching"
function make_child_node_with_rf_branch(node::Node)

    child = deepcopy(node)
    J, E, w, W, S, bounds = get_node_parameters(child)

    # get largest item in bag and the fractional item to bound on


    weight_i, i = heaviest_in_address([1], item_address, w)
    weight_j, j = heaviest_in_address([2], item_address, w)

    # merge branch
    new_J, new_w = merge_items(i, j, J, w, item_address)

    # split branch
    push!(E, sort([i,j]))


    return child
end

"makes children with bag branching"
function make_child_node_with_bag_branch(node::Node)
    
    child = deepcopy(node)
    J, E, w, W, S, bounds = get_node_parameters(child)




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

"apply First Fit Decreasing heuristic"
function first_fit_decreasing(J, w, W, E; verbose=true)
    bags = Array{Int64}[Int64[0 for j in J] for i in J]
    # bags = Array{Bool}[Bool[0 for j in J] for i in J]
    bags_sizes = Int64[0 for i in J]
    bags_slacks = Int64[W for i in J]
    bags_conflicts = Array{Bool}[Bool[0 for j in J] for i in J]

    J = sort(J, rev=true, by=(x)->w[x])

    # [conflicts of i for each i in J]
    conflicts = get_edges(J, E)

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
function price_lp(pi_bar, w, W, J, E, S; verbose=3, epsilon=1e-4)
    price = Model(GLPK.Optimizer)
    @variable(price, 1 >= x[1:length(J)] >= 0)
    @constraint(price, sum([w[j]*x[j] for j ∈ J]) <= W, base_name="capacity")
    for e in E
        @constraint(price, x[e[1]] + x[e[2]] <= 1, base_name="e_($(e[1]), $(e[2]))")
    end
    
    # @objective(price, Min, sum([(1- pi_bar[j])*x[j] for j ∈ J]))
    @objective(price, Min, 1- sum([pi_bar[j]*x[j] for j ∈ J]))
    set_silent(price)

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
        
    return p_obj, value.(price[:x])
end

"Runs pricing linear programming, but constrains new lambda to be integer"
function int_price_lp(pi_bar, w, W, J, E, S; verbose=3, epsilon=1e-4)
    price = Model(GLPK.Optimizer)
    @variable(price, 1 >= x[1:length(J)] >= 0)
    @constraint(price, sum([w[j]*x[j] for j ∈ J]) <= W, base_name="capacity")
    for e in E
        @constraint(price, x[e[1]] + x[e[2]] <= 1, base_name="e_($(e[1]), $(e[2]))")
    end
    
    # @objective(price, Min, sum([(1- pi_bar[j])*x[j] for j ∈ J]))
    @objective(price, Min, 1- sum([pi_bar[j]*x[j] for j ∈ J]))
    set_silent(price)

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

function cga(master, price_function, w, W, J, E, lambdas, S, S_len; verbose=3, max_iter=10e2, epsilon=1e-4)
    
    m_obj = Inf

    # run price, add new columns, check solution, repeat if necessary
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
        p_obj, q = price_function(pi_bar, w, W, J, E, S, verbose=verbose, epsilon=epsilon)

        if p_obj < -epsilon

            # price
            verbose >= 1 && println("adding lambda: $(q)")

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
            verbose >= 2 && println(master)

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
    2: global done;
"""
function update_bounds_status(node, bounds, best_node; verbose=1)

    status = 0
    
    if node.bounds[2] < bounds[2] # is there a new best solution?
        status = 1
        bounds[2] = node.bounds[2]
        best_node = node.id

        if bounds[1] == bounds[2] # is the new solution guaranteed to be globally optimal?
            status = 2
        end

    elseif node.bounds[1] == node.bounds[2] # is it locally optimal? 
        status = 1

    elseif node.bounds[1] >= bounds[2] # should the algorithm stop processing the node?
        status = 1
    end

    return status
end

function solve_bpc(
    J::Array{Int64}, 
    E::Array{Int64, Int64}, 
    w::Array{Int64}, 
    W::Int64; 
    verbose::Int64=1, 
    run_ffd::Bool=true, 
    epsilon::Float64=1e-4,
    )
    
    item_amount = length(J)

    # [LB, UB]
    bounds = [1, item_amount+1]

    # where can item j be found? (for merging items)
    base_item_adress = Int64[j for j in J]
    
    # initialize node list
    nodes = Node[Node(1, J, E, w, W, Float32[], Int64[-1 for q in S], base_item_adress, false, deepcopy(bounds), 0, Array{Int64}[])]
    queue = Int64[1]
    
    best_node = Int64[1]

    # Start the tree
    while !(isempty(queue))

        # get next node
        next_node_id = splice!(queue, 1)
        node = nodes[next_node_id]
        J, E, w, W, S, bounds = get_node_parameters(node)
        verbose >=1 && println("node $(node.id)")

        ## first try solving the node with heuristics and bounds' properties

        # naive solution (one item per bag)
        naive_solution = get_naive_solution(J)
        node.bounds[2] = item_amount

        verbose >= 1 && println("Naive upper bound: $(node.bounds[2])")

        # get initial lower bound ( ⌈∑w/W⌉ ) 
        node.bounds[1] = get_simple_lower_bound(w, W)
        verbose >= 1 && println("⌈∑w/W⌉ lower bound: $(node_lb)")

        # update bounds status
        optimality = update_bounds_status(node, bounds, verbose=verbose)
        if optimality != 0 # is it a global or local optimal?
            if optimality == 1 # local optimal
                # prune the tree
                continue
            else # global optimal
                return naive_solution, bounds[2]
            end
        end

        # FFD heuristic for initial solution and upper bound
        if run_ffd
            ffd_solution, ffd_upper_bound = first_fit_decreasing(J, w, W, E, verbose=verbose>1)
            verbose >= 1 && println("FFD heuristic upper bound: $(ffd_upper_bound)")
            
            # if a better solution than one item per bag was found 
            if ffd_upper_bound < node.bounds[2]
                node.bounds[2] = ffd_upper_bound
                
                # update bounds status
                optimality = update_bounds_status(node, bounds, verbose=verbose)
                if optimality != 0 # is it a global or local optimal?
                    if optimality == 1 # local optimal
                        # prune the tree
                        continue
                    else # global optimal
                        return ffd_solution, bounds[2]
                    end
                end   
            end
    
            initial_solution = deepcopy(ffd_solution)
        else
            initial_solution = deepcopy(naive_solution)
        end
        
        # try to improve lower bound with martello L2 lower bound
        for alpha in get_not_greater_than_half_capacity(w, W)
            lower_bound = get_l2_lower_bound(alpha, W, w)
            if lower_bound > node.bounds[1]

                node.bounds[1] = lower_bound
                verbose >= 1 && println("L2 lower bound with α = $(alpha): $(node_lb)")

                # update bounds status
                optimality = update_bounds_status(node, bounds, verbose=verbose)
                if optimality != 0 # is it a global or local optimal?
                    if optimality == 1 # local optimal
                        # prune the tree
                        continue
                    else # global optimal
                        return initial_solution, bounds[2]
                    end
                end    
            end
        end
        

        ## no easy way out, build the LP and try to solve with integer pricing

        # store best solution
        best_solution = deepcopy(initial_solution)
        
        # build master
        master = Model(GLPK.Optimizer)
        set_silent(master)
        
        # add the naive solution as lambda variables (will serve as artificial variables)
        S = Array{Float32}[q for q in naive_solution]

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
        i = 1
        lambdas = VariableRef[]
        for q in S 
            var = @variable(master, lower_bound=0, base_name="λ_$(i)")
            push!(lambdas, var)
            i+=1
        end
    
        # demand constraints
        for i in J
            @constraint(master, sum([sum(S[q][i]*lambdas[q]) for q in 1:S_len]) >= 1, base_name="demand_$(i)")
        end
    
        # objective function
        @objective(master, Min, sum(lambdas))
    
        # show initial master
        verbose >= 2 && println(master)
    
        # run column generation with specialized pricing

        # run column generation with integer pricing
        m_obj, cga_ub, S_len = cga(master, int_price_lp, w, W, J, E, lambdas, S, S_len, verbose=verbose, epsilon=epsilon, max_iter=1e2)
    
        # get solution values
        lambda_bar = value.(lambdas)
        x_bar, cga_ub = get_x(lambda_bar, S, S_len, epsilon=epsilon)
    
        # treat current solution
        current_solution = round_up_solution(x_bar)
        current_solution, cga_ub = prune_excess_with_priority(current_solution, J, w, epsilon=epsilon)
        
        verbose >= 1 && println("Integer CGA upper bound: $(cga_ub)")
    
        # was there an improvement from the heuristic?
        if cga_ub < node.bounds[2]
    
            node.bounds[2] = cga_ub
            best_solution = deepcopy(current_solution)
            
            # update bounds status
            optimality = update_bounds_status(node, bounds, verbose=verbose)
            if optimality != 0 # is it a global or local optimal?
                if optimality == 1 # local optimal
                    # prune the tree
                    continue
                else # global optimal
                    return best_solution, bounds[2]
                end
            end  
        end


        ## BCPA

        # apply cga
        z, cga_lb, node.S_len = cga(node.master, price_lp, w, W, J, E, lambdas, node.S, node.S_len)

        # is there already a better or equal solution?
        if cga_lb >= bounds[2]
            continue
        end

        if cga_lb > node.bounds[1]
            verbose >= 1 && println("CGA lower bound: $(cga_lb)")

            node.bounds[1] = cga_lb

            # update bounds status
            optimality = update_bounds_status(node, bounds, verbose=verbose)
            if optimality != 0 # is it a global or local optimal?
                if optimality == 1 # local optimal
                    # prune the tree
                    continue
                end
            end  
        end

        # is there potential for a better solution? 
        if cga_lb < node.bounds[2]

            # get lambda values of the solution
            node.lambda_bar = value.(node.lambdas)
            
            # get where to bound
            bound_on = most_fractional_on_vector(node.lambda_bar)

            if bound_on == -1 # λ is integer, check x
                
                x_bar = get_x(node.lambda_bar, S, S_len; epsilon=epsilon)


                # node is finished (integer solution found)
                
            else # 


            end

        else # node is finished (no point in continuing)



        end




    end
    
    verbose >= 1 && println("tree finished")



    

    lambda_bar = value.(lambdas)
    println([i for i in lambda_bar if i > epsilon])

    solution = Array{Float32}[]
    for i in 1:length(lambda_bar)
        if lambda_bar[i] > epsilon
            println("lambda $(i): $(lambda_bar[i]) -> $(S[i])")
            push!(solution, lambda_bar[i]*S[i])
        end
    end

    
    solution = get_pretty_solution(solution, 5, J)

    return solution, bounds[2]


    return false, false
end