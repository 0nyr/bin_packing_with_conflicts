using JuMP
using Gurobi
using LinearAlgebra
using Logging

global LOG_IO = stdout

"Returns set of items found in a set of addresses"
function unmerge_bag_items(addresses, item_address)
    
    items = Int64[]
    for (j, address) in enumerate(item_address)
        if address ∈ addresses
            push!(items, j)
        end
    end
    return items
end

"Translate a bag to account for mergings (returns binarized)"
function merge_bag_items(bag, item_address, J)
    merged = Int64[0 for j in J]
    for j in bag
        merged[item_address[j]] = 1
    end

    return merged
end


"returns item and weight of heaviest item in a set of addresses"
heaviest_in_address(addresses, item_address, w) = findmax(x -> w[x], unmerge_bag_items(addresses, item_address))

"returns item and weight of lightest item in a set of addresses"
lightest_in_address(addresses, item_address, w) = findmin(x -> w[x], unmerge_bag_items(addresses, item_address))

"returns an array where if i and j are in conflict, j ∈ array[i] and i ∈ array[j]"
function get_edges(J, E)    
    edges = Vector{Int64}[Int64[] for i in J]
    for i in J
        for edge in E
            if i == edge[1]
                push!(edges[i], edge[2])
            elseif i == edge[2]
                push!(edges[i], edge[1])
            end
        end
    end

    # for e in E
    #     push!(edges[e[1]], e[2])
    #     push!(edges[e[2]], e[1])
    # end
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

"Utility to find most fractional item on a vector"
function most_fractional_on_vector(v; epsilon=1e-4)
    bound_on = -1
    closest = 1
    for (n, j) in enumerate(v)
        diff = j - floor(j)
        if diff > epsilon && diff < 1-epsilon
            d = abs(diff - 0.5) 
            if d < closest
                closest = d
                bound_on = n  
            end
        end
    end
    return bound_on
end

"Utility to find most fractional bag and most fractional item in a solution"
function find_most_fractional_bag_and_item(bags_in_use, lambda_bar, S, S_len, conflicts, J; epsilon=1e-4)
    most_fractional_bag = -1
    most_fractional_item = [-1, -1]
    bag_closest = 1
    item_closest = 1

    for q in bags_in_use
        # println(LOG_IO, "checking integrality of $(q)")

        # check lambda integrality
        d = lambda_bar[q]
        diff = lambda_bar[q] - floor(lambda_bar[q])
        
        if diff > epsilon && diff < 1-epsilon

            d = abs(diff - 0.5) 
            if d < bag_closest
                bag_closest = d
                most_fractional_bag = q  
            end
        else # bag is integer
            continue
        end
        # println(LOG_IO, "λq: $(lambda_bar[q]), diff: $(diff), most_fractional_bag: $(most_fractional_bag)")


        # check items integrality
        for (j, x_j) in enumerate(S[q])

            d = lambda_bar[q]*x_j
            # d = x_j
            diff = d - floor(d)
            
            if diff > epsilon && diff < 1-epsilon

                d = abs(diff - 0.5) 
                if d < item_closest
                    item_closest = d
                    most_fractional_item = [q, j]  
                end
            end
            # println(LOG_IO, "x$(j): $(x_j), diff: $(diff), is_bag_integer: $(is_bag_integer)")
        end

    end

    return most_fractional_bag, most_fractional_item
end

"Analyze branching possibilities"
function make_branching_analysis(bags_in_use::Vector{Int64}, lambda_in_use::Vector{Float32}, lambda_bar::Vector{Float64}, S::Vector{Vector{Float32}}, S_len::Int64, conflicts::Vector{Vector{Int64}}, J::Vector{Int64}, w::Vector{Int64}; epsilon::Float64=1e-4)
    
    
    # find most fractional item (by weight) and most fractional bag 
    most_fractional_bag = -1
    most_fractional_item = [-1, -1]
    bag_closest = 1
    item_closest = 1

    # bags in use that are fractional
    # fractional_bags = Int64[]
    for q in bags_in_use

        # if !(isapprox(lambda_bar[q], 0, atol=epsilon)) || !(isapprox(lambda_bar[q], 1, atol=epsilon))
        #     push!(fractional_bags, q)
        # end

        d = lambda_bar[q]
        diff = lambda_bar[q] - floor(lambda_bar[q])
        
        if diff > epsilon && diff < 1-epsilon
    
            d = abs(diff - 0.5) 
            if d < bag_closest
                bag_closest = d
                most_fractional_bag = q  
            end
        end
    end

    # # fractional bag weights
    # fractional_bags_weights = Float64[sum(S[q].*w) for q in fractional_bags]

    # # just read the code...
    # frac_bags_frac_weights = fractional_bags_weights.*Float64[lambda_bar[q] for q in fractional_bags]

    # # amount of *items* in each bag
    # bag_item_amount = Float32[sum(S[q]) for q in fractional_bags]

    # # conflicts of each bag
    # bag_conflicts = Vector{Int64}[Int64[0 for j in J] for i in bags_in_use]

    # # amount of *edges* in each bag
    # bag_edges_amount = Int64[0 for i in fractional_bags]
    # for (bag_i, q) in enumerate(fractional_bags)
    #     for (j, val) in enumerate(q)
    #         if val > epsilon
    #             bag_edges_amount[bag_i] += length(conflicts[j])
    #             # for item_k in conflicts[j]
    #             #     bag_conflicts[bag_i][item_k] = 1
    #             # end
    #         end
    #     end
    # end

    # # amount of *conflicts* in each bag
    # bag_conflicts_amount = Int64[sum(i) for i in bag_conflicts]

    return most_fractional_bag
end

"returns {i | λ_i > 0 ∀ i}"
function get_bags_in_use(lambda_bar, S, S_len, J; epsilon=1e-4)
    # bags = Vector{Float32}[Float32[0.0 for j in J] for i in J]
    bags_in_use = Int64[]
    lambdas_in_use = Float32[]

    # get bags selected for use
    for q in 1:S_len
        if lambda_bar[q] > epsilon 
            push!(bags_in_use, q)
            push!(lambdas_in_use, q)
        end
    end
    
    return bags_in_use, lambdas_in_use
end

"from lambda, returns x"
function get_x(lambda_bar, S, S_len, J; epsilon=1e-4)
    bags = Vector{Float32}[Float32[0.0 for j in J] for i in J]

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
round_up_solution(solution; epsilon=1e-4) = Vector{Int64}[Int64.(ceil_vector(bag, epsilon=epsilon)) for bag in solution]

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
    item_locations = Vector{Int64}[[] for j in J]
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
    bags_weights = Int64[sum([w[j] for (j, value) in enumerate(bag) if value > .5]) for bag in solution]
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
    naive_solution = Vector{Int64}[Int64[0 for j in J] for i in J]
    for i in J # one item per bag
        naive_solution[i][i] = 1
    end
    return naive_solution
end

"remove forbidden bags from a solution. Returns true if bags were removed"
function remove_forbidden_bags(solution::Vector{Vector{Int64}}, forbidden_bags::Vector{Vector{Int64}})

    original_length = length(solution)

    filter!(x -> !(x ∈ forbidden_bags), solution)

    return original_length == length(solution)
end

"translate edges for new address. removes duplicates and edges involving items in mandatory bags"
function translate_edges(original_E, item_address)
    return unique(Vector{Int64}[
        sort([item_address[e[1]], item_address[e[2]]]) for e in original_E 
        if item_address[e[1]] != 0 && item_address[e[1]] != 0 && item_address[e[1]] != item_address[e[2]] 
    ]) 
end
 

"translates a solution, unmerging items and adding mandatory bags"
function translate_solution(node; epsilon=1e-4)

    # println(LOG_IO, "mandatory_bags: $(node.mandatory_bags)")

    translated_solution = Vector{Int64}[Int64[0 for j in node.item_address] for i in 1:node.bounds[2]] 

    # items at address i, for each i
    address_items = Vector{Int64}[Int64[0 for j in node.item_address] for i in node.J]

    for (j, address) in enumerate(node.item_address)
        if address != 0
            address_items[address][j] = 1
        end 
    end


    for (i, bag) in enumerate(node.solution)
        for (address, is_here) in enumerate(bag)

            if is_here > epsilon
                translated_solution[i] += address_items[address]
            end
        end
    end

    # convert mandatory_bags to constant length, binary arrays
    mandatory_bags_binary = Vector{Int64}[Int64[0 for j in node.item_address] for i in 1:node.mandatory_bag_amount]
    for (i, bag_items) in enumerate(node.mandatory_bags)
        for j in bag_items
            mandatory_bags_binary[i][j] = 1     
        end
    end

    # add the mandatory_bags
    # translated_solution = vcat(translated_solution, mandatory_bags_binary)
    for (i, j) in enumerate(node.bounds[2]-node.mandatory_bag_amount+1:node.bounds[2])
        translated_solution[j] = mandatory_bags_binary[i]
    end

    # println(LOG_IO, "translated_solution: $(translated_solution)")

    return translated_solution
end

"transforms solution structure from binary, same length arrays to integer, variable length arrays"
function get_pretty_solution(bags, bags_amount; epsilon=1e-4)
    return Vector{Int64}[ Int64[j for j in 1:length(bags[i]) if bags[i][j] > 0] for i in 1:bags_amount ]
end

get_demand_constraints(model, J) = [constraint_by_name(model, "demand_$(i)") for i in J]
reduced_cost(x, pi_bar, J) = 1 - sum([pi_bar[j]*x[j] for j ∈ J])

"Utility function for retrieving master data necessary for the pricing step"
function get_master_data_for_pricing(master, J; verbose=2)
    m_obj = objective_value(master)
    verbose >= 2 && println(LOG_IO, "Z = $(m_obj)")

    demand_constraints = get_demand_constraints(master, J)
    pi_bar = dual.(demand_constraints)

    return m_obj, demand_constraints, pi_bar
end
