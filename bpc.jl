using JuMP
using GLPK
using LinearAlgebra

function get_edges(J, E)    
    edges = Array{Int64}[Int64[] for i in J]
    for i in J
        for edge in E
            if i == edge[1]
                push!(edges[i], edge[2])
            elseif i == edge[2]
                push!(edges[i], edge[1])
            end
        end
    end
    return edges
end

"returns simplest lower bound ( ceil(sum(w)/W) )"
function get_simple_lower_bound(w, W)
    return Int64(ceil(sum(w)/W))
end


function get_not_greater_than_half_capacity(w, W)
    return Int64[w_i for w_i in sort(w, rev=true) if w_i <= W/2]
end

"L2 bound - Martello1990"
function get_l2_lower_bound(alpha, W, w)
    n1 = Int64[i for i in w if i > W-alpha]
    n2 = Int64[i for i in w if W-alpha >= i && i > W/2]
    n3 = Int64[i for i in w if W/2 >= i && i >= alpha]

    return Int64(length(n1) + length(n2) + max(0, ceil( (sum(n3) - (length(n2)*W - sum(n2) ) )/W  ) )) 
end

function get_bag_weight(bag, w)
    return sum([w[i] for i in bag])
end

function get_slack(bag, W, w)
    return W - get_bag_weight(bag, w)
end

function get_compatible_items(bag, J, bag_conflicts, w, slack)
    return Int64[i for i in J if i ∉ bag_conflicts && slack - w[i] >= 0]
end

"transforms solution structure from binary, same length arrays to integer, variable length arrays"
function get_pretty_solution(bags, bags_amount, J)
    # clean_bags = [[bags[i][j] for j in 1:bags_sizes[i]] for i in 1:bags_amount]
    return Array{Int64}[ Int64[j for j in 1:length(J) if bags[i][j] > 0] for i in 1:bags_amount ]
end

"Utility to find most fractional item on a vector"
function most_fractional_on_bag(solution; epsilon=1e-4)
    bound_on = -1
    closest = 1
    for j in 1:length(solution)
        diff = solution[j] - floor(solution[j])
        if diff > epsilon && diff < 1-epsilon
            d = abs(diff - 0.5) 
            if d < closest
                closest = d
                bound_on = j  
            end
        end
    end
    return bound_on
end

"Utility to find most fractional x in a solution"
function most_fractional_on_solution(solution; epsilon=1e-4)

end

"from lambda, returns x"
function get_x(lambda_bar, S, S_len; epsilon=1e-4)
    bags = Array{Float32}[Float32[0.0 for j in J] for i in J]
    # bags_conflicts = Array{Bool}[Bool[0 for j in J] for i in J]
    # bags_sizes = Int64[0 for i in J]

    # get bags selected for use
    bag_amount = 0
    for q in 1:S_len
        if lambda_bar[q] > epsilon
            bag_amount+=1
            bags[bag_amount] = lambda_bar[q]*S[q]
        end
    end
    
    return bags[1:bag_amount], bag_amount
end

floor_vector(q; epsilon=1e-4) = floor.(q .+ epsilon)
ceil_vector(q; epsilon=1e-4) = ceil.(q .- epsilon)

"rounds up the solution bags, converting to integer"
round_up_solution(solution; epsilon=1e-4) = Array{Int64}[Int64.(ceil_vector(bag, epsilon=epsilon)) for bag in solution]


"Prunes the excess items, prioritizing the most heavy bags"
function prune_excess_with_priority(solution, J, w; epsilon=1e-4)
    
    item_count = sum(solution)
    excess = Int64[]

    # find excess items
    for i in J
        if item_count[i] > 1 + epsilon
            push!(excess, i)
        end
    end

    if isempty(excess) # no items to prune
        return solution, length(solution)
    end
    
    # get location of all excesses
    bag_i = 0
    item_locations = Array{Int64}[[] for j in J]
    for bag in solution
        bag_i += 1

        # if an item in excess is in the bag, register the bag number
        for j in excess
            if bag[j] > epsilon 
                push!(item_locations[j], bag_i)
            end
        end
    end

    # remove item from bags, prioritizing the most heavy bags
    bags_weights = Int64[sum([w[j] for j in bag]) for bag in solution]
    for j in excess

        # sort relevant bags by most empty first
        most_empty_first = sort(item_locations[j], by=(x)->bags_weights[x])
        
        for i in most_empty_first[2:end] # remove excess amount
            solution[i][j] = 0
        end
    end
    
    solution = solution[[i for i in 1:length(solution) if sum(solution[i]) > 0]]

    return solution, length(solution)
end

function next_fit_decreasing(J, w, W, E)
    bags = Array{Int64}[]
    J = sort(J, rev=true, by=(x)->w[x])
    J0 = Int64[i for i in J]

    conflicts = get_edges(J, E)

    for iter_i in 1:10e5
        bag = Int64[]
        push!(bag, splice!(J, 1))
        # println("new bag: $(bag)")
        

        bag_weight = get_bag_weight(bag, w)
        slack = W - bag_weight

        bag_conflicts = [i for i in conflicts[bag[1]]]

        cJ = get_compatible_items(bag, J, bag_conflicts, w, slack)
        # println("compatible items: $(cJ)")



        i = 1
        while i <= length(cJ)
            if slack - w[cJ[i]] >= 0 # add all items that fit in descending order
                push!(bag, splice!(cJ, i))

                bag_weight = get_bag_weight(bag, w)
                slack = W - bag_weight

                bag_conflicts = vcat(bag_conflicts, conflicts[bag[end]]) 
        
                cJ = get_compatible_items(bag, J, bag_conflicts, w, slack)
                i=1
            else
                i+=1
            end
        end

        push!(bags, bag)
        filter!(i->i ∉ bag, J)

        # println("new closed bag: $(bag)")
        # println("remaining items: $(J)")
        # println("remaining items' weights: $([w[i] for i in J])")

        if isempty(J)
            println("Done")
            break
        end
    end
    println("Final bags: $([i for i in bags])")
    println("Final items weights: $([[w[j] for j in i] for i in bags])")
    return bags
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

    # println("Final bags: $([i[] for i in bags])")
    # println("Final items weights: $([[w[j] for j in i] for i in bags])")
end

"Naive solution, one item per bag"
function get_naive_solution(J)
    naive_solution = Array{Int64}[Int64[0 for j in J] for i in J]
    for i in J # one item per bag
        naive_solution[i][i] = 1
    end
    return naive_solution
end

get_demand_constraints(model, J) = [constraint_by_name(model, "demand_$(i)") for i in J]
reduced_cost(x, pi_bar, J) = 1 - sum([pi_bar[j]*x[j] for j ∈ J])

"Utility function for retrieving master data necessary for the pricing step"
function get_master_data_for_pricing(master, J; verbose=2)
    m_obj = objective_value(master)
    verbose >= 2 && println("Z = $(m_obj)")

    demand_constraints = get_demand_constraints(master, J)
    pi_bar = dual.(demand_constraints)

    return m_obj, demand_constraints, pi_bar
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

    cga_lower_bound = Int(ceil(m_obj - epsilon))

    return m_obj, cga_lower_bound, S_len
end

function solve_bpc(J, E, w, W; verbose=1, run_ffd=true, epsilon=1e-4)

    # naive solution (one item per bag)
    naive_solution = get_naive_solution(J)
    UB = length(J)
    verbose >= 1 && println("Naive upper bound: $(UB)")
    
    # FFD heuristic for initial solution and upper bound
    if run_ffd
        ffd_solution, ffd_upper_bound = first_fit_decreasing(J, w, W, E, verbose=verbose>1)
        verbose >= 1 && println("FFD heuristic upper bound: $(ffd_upper_bound)")
        
        # if a better solution than one item per bag was found 
        if ffd_upper_bound < UB
            UB = ffd_upper_bound
        end

        initial_solution = deepcopy(ffd_solution)
    else
        initial_solution = deepcopy(naive_solution)
    end
    
    # get initial lower bound ( ⌈∑w/W⌉ ) 
    LB = get_simple_lower_bound(w, W)
    verbose >= 1 && println("⌈∑w/W⌉ lower bound: $(LB)")
    
    # if an optimal solution was found
    LB == UB && return initial_solution, UB
    
    # try to improve lower bound with martello L2 lower bound
    for alpha in get_not_greater_than_half_capacity(w, W)
        lower_bound = get_l2_lower_bound(alpha, W, w)
        if lower_bound > LB
            LB = lower_bound
            verbose >= 1 && println("L2 lower bound with α = $(alpha): $(LB)")
            
            # if an optimal solution was found
            LB == UB && return initial_solution, UB
        end
    end
    
    ## no easy way out, build the LP
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
    
    verbose >= 1 && println("Integer CGA: $(cga_ub)")

    # was there an improvement from the heuristic?
    if cga_ub < UB

        UB = cga_ub
        best_solution = deepcopy(current_solution)

        if LB == UB
            return best_solution, UB
        end
    end
    
    # Truly, no easy way out, do BCPA

    





    

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

    return solution, UB


    return false, false
end