using JuMP
using GLPK
using LinearAlgebra

"return set of items found in a set of addresses"
function get_items_in_address(addresses, item_address)
    
    items = Int64[]
    for (j, address) in enumerate(item_address)
        if address ∈ addresses
            push!(items, j)
        end
    end
    return items
end

"returns item and weight of heaviest item in a set of addresses"
heaviest_in_address(addresses, item_address, w) = findmax(x -> w[x], get_items_in_address(addresses, item_address))

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


function get_not_greater_than_half_capacity(w, W)
    return Int64[w_i for w_i in sort(w, rev=true) if w_i <= W/2]
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
function most_fractional_on_vector(solution; epsilon=1e-4)
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

"Naive solution, one item per bag"
function get_naive_solution(J)
    naive_solution = Array{Int64}[Int64[0 for j in J] for i in J]
    for i in J # one item per bag
        naive_solution[i][i] = 1
    end
    return naive_solution
end

"translate edges for new address"
function translate_edges(original_E, item_address)
    return unique(Array{Int64}[sort([item_address[e[1]], item_address[e[2]]]) for e in original_E])
end

"translates a solution, unmerging items"
function translate_solution(node, solution; epsilon=1e-4)
    translated_solution = Array{Int64}[Int64[0 for j in node.item_address] for i in 1:node.bounds[2]] 

    address_items = Array{Int64}[Int64[0 for j in node.item_address] for i in node.J]

    for (j, address) in enumerate(node.item_address) 
        address_items[address][j] = 1
    end


    for (i, bag) in enumerate(solution)
        for (address, is_here) in enumerate(bag)

            if is_here > epsilon
                translated_solution[i] += address_items[address]
            end
        end
    end

    return translated_solution
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
