using JuMP
using GLPK
using LinearAlgebra


"Solves with Revised Simplex Algorithm"
function rsa(A, b, c, integer=false)
    m = Model(GLPK.Optimizer)

    if integer
        @variable(m, x[1:4] >= 0, Int)
    else
        @variable(m, x[1:4] >= 0)
    end
    @objective(m, Min, dot(c,x))
    @constraint(m, c1, sum(A[1:1, :]*x) == b[1])
    @constraint(m, c2, A[2:end, :]*x .<= b[2:end])
    print(m)
    
    optimize!(m)
    
    solution = value.(x)
    obj = objective_value(m) 
    println("objective value = $(obj)")
    return solution, obj, m
end

"Solves with Branch-and-Bound Algorithm"
function bba(model, verbose=1)
    upper_bound = Inf
    best_solution = nothing
    best_model = nothing
    m, m_ref = copy_model(model)
    L = [m]
    t = 0
    tracker = ["S$(t)"]


    while !(isempty(L))
        m = splice!(L, 1)
        tag = splice!(tracker, 1)

        println("next in tree: $(tag)")
        if verbose >= 2
            print(m)
        end
        
        set_optimizer(m, GLPK.Optimizer)
        optimize!(m)

        if termination_status(m) == OPTIMAL
            z = objective_value(m)
            if verbose >= 1
                println("optimized. z = $(z). UB = $(upper_bound)")
            end
        else
            if verbose >= 1
                println("infeasible")
            end
            z = Inf
        end

        if z < upper_bound

            # check if x is integer
            # if not, assign most fractional to bound_on
            solution = value.(m[:x])
            bound_on = most_fractional_x(solution)
            
            if verbose >= 1
                println("solution = $(solution)")
            end

            # if integer
            if (bound_on == -1)
                if verbose >= 1
                    println("solution is integer")
                end
                upper_bound = z
                best_solution = solution
                best_model, _ = copy_model(m)

                if verbose >= 1
                    println("best solution: $(best_solution)")
                    println("UB = $(upper_bound)")
                end

                # Pruning the tree
                cleaning = true
                i=1
                j=0
                while cleaning
                    # println(i, "\t", L)
                    if i == length(L)
                        cleaning = false
                    end
                    if termination_status(L[i]) == OPTIMAL && upper_bound < objective_value(L[i])
                        splice!(L, i)
                        j+=1
                    else
                        i += 1
                    end
                end
                if verbose >= 1
                    println("pruned $(j) nodes")
                end
            else
                if verbose >= 1
                    println("solution is not integer")
                end

                # Create subproblems and add bounds
                # x >= ceil(x*)
                new_m, new_m_ref = copy_model(m)
                @constraint(new_m, new_m[:x][bound_on] >= ceil(solution[bound_on]) )
                push!(L, new_m)
                t+=1
                push!(tracker, "S$(t)")
                println("added subproblem $(tracker[end]) to tree")
                if verbose >= 2
                    print(new_m)
                end
                

               
                # x <= floor(x*)
                new_m, new_m_ref = copy_model(m)
                @constraint(new_m, new_m[:x][bound_on] <= floor(solution[bound_on]) )
                push!(L, new_m)
                t+=1
                push!(tracker, "S$(t)")
                println("added subproblem $(tracker[end]) to tree")
                if verbose >= 2
                    print(new_m)
                end
        
            end
        end
    end
    if verbose >= 1
        println("done. final model:")
        print(best_model)
    end
    return best_solution, upper_bound
end

"Solves with Branch-and-Price Algorithm"
function bpa(A, P, bm, bp, c, master=nothing, price=nothing, verbose=1, cga_verbose=1)
    upper_bound = Inf
    best_solution = nothing
    best_model = nothing

    A_nrows = size(A)[1]
    P_nrows = size(P)[1]
    x_amount = length(c)

    solved = Dict()

    L = Any[[A, P, bm, bp, c, master, price, nothing]]
    t = 0
    tracker = ["S$(t)"]

    while !(isempty(L))
        A, P, bm, bp, c, master, price, S = splice!(L, 1)
        tag = splice!(tracker, 1)

        solved[tag] = [A, P, bm, bp, c, master, price, S]

        if verbose >= 1
            println("next in tree: $(tag)")
        end
        # print(master)

        solution, m_obj, S, master, price = cga(
            A, 
            P, 
            bm, 
            bp, 
            c; 
            S=S,
            master=master, 
            price=price, 
            verbose=cga_verbose, 
            solved=solved)

        # set_optimizer(master, GLPK.Optimizer)
        # optimize!(master)

        if termination_status(master) == OPTIMAL
            z = objective_value(master)
            if verbose >= 1
                println("optimized. z = $(z). UB = $(upper_bound)")
            end
        else     
            if verbose >= 1
                println("infeasible")
            end
            z = Inf
        end

        if z < upper_bound

            # check if x is integer
            # if not, assign most fractional to bound_on
            # solution = value.(master[:x])
            bound_on = most_fractional_x(solution)
                
            if verbose >= 1
                println("solution = $(solution)")
            end

            # if integer
            if (bound_on == -1)
                
                if verbose >= 1
                    println("solution is integer")
                end
                upper_bound = z
                best_solution = solution
                best_model = copy(master)
                best_model_data = [A, P, bm, bp, c, best_model, copy(price), S]

                        
                if verbose >= 1
                    println("best solution: $(best_solution)")
                    println("UB = $(upper_bound)")
                end

                # Pruning the tree
                cleaning = !(isempty(L))
                i=1
                j=0
                while cleaning
                    # println(i, "\t", L)
                    if i == length(L)
                        cleaning = false
                    end
                    if termination_status(L[i][6]) == OPTIMAL && upper_bound < objective_value(L[i][6])
                        splice!(L, i)
                        j+=1
                    else
                        i += 1
                    end
                end
                        
                if verbose >= 1
                    println("pruned $(j) nodes")
                end
            else
                        
                if verbose >= 1
                    println("solution is not integer")
                    print("branching: ")
                end

                # lambda = [variable_by_name(master, "λ_$(i)") for i in 1:length(S)]
                # println(value.(lambda))
                # println(S)

                
                v_ceil = ceil(solution[bound_on])
                v_floor = floor(solution[bound_on])

                if verbose >= 1
                    println(" x[$(bound_on)] ∉ ($(v_floor), $(v_ceil))")
                end

                new_row = Matrix(UniformScaling(0), 1, x_amount)   
                new_row[bound_on] = 1
                new_A = vcat(A, new_row)

                # Create subproblems and add bounds
                # x >= ceil(x*)
                new_m = copy(master)

                m_c = get_Ax_constraints(new_m, size(new_A)[1])
                new_con_i = length(m_c)

                # retrieve lambda references 
                lambda = [variable_by_name(new_m, "λ_$(i)") for i in 1:length(S)]

                # @constraint(new_m, new_m[:x][bound_on] >= ceil(solution[bound_on]) )
                new_cons = @constraint(new_m, sum( [ S[q][bound_on]*lambda[q] for q in 1:length(S)] ) >= v_ceil, base_name="m_c[$(new_con_i)]")
                new_bm = vcat(bm, v_ceil)
                
                push!(L, [new_A, P, new_bm, bp, c, new_m, price, S])
                
                t+=1
                push!(tracker, "S$(t)")             
                if verbose >= 1
                    println("added subproblem $(tracker[end]) to tree:")
                    print(new_m)
                end

               
                # x <= floor(x*)
                new_m = copy(master)
                
                # retrieve lambda references 
                lambda = [variable_by_name(new_m, "λ_$(i)") for i in 1:length(S)]

                # @constraint(new_m, new_m[:x][bound_on] <= floor(solution[bound_on]) )
                new_cons = @constraint(new_m, sum( [ S[q][bound_on]*lambda[q] for q in 1:length(S)] ) <= v_floor, base_name="m_c[$(new_con_i)]")
                new_bm = vcat(bm, v_floor)
                
                push!(L, [new_A, P, new_bm, bp, c, new_m, price, S])
                
                t+=1
                push!(tracker, "S$(t)")             
                if verbose >= 1
                    println("added subproblem $(tracker[end]) to tree:")
                    print(new_m)
                end



            end
        end
    end
    println("done. final model:")
    print(best_model)
    return best_solution, upper_bound
end

"Utility to find most fractional item on a vector"
function most_fractional_x(solution, epsilon=1e-4)
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

function get_Ax_constraints(model, A_nrows)
    return [constraint_by_name(model, "m_c[$(i)]") for i in 1:A_nrows]
end

"Solves with Column Generation Algorithm"
function cga(A, P, bm, bp, c; master=nothing, price=nothing, S=nothing, M=9e6, verbose=1, solved=nothing)

    A_nrows = size(A)[1]
    P_nrows = size(P)[1]
    x_amount = length(c)

    if isnothing(S)
        S = []
    else
        S = copy(S)
    end
    Aq = [A*q for q in S]
    cq = [dot(c,q) for q in S] 

    ### build the master
    if isnothing(master)
        master = Model(GLPK.Optimizer)

        # artificial variable
        @variable(master, a[1:A_nrows] >= 0)

        # constraints 1 and convexity with artificial variables only
        @constraint(master, [i in 1:A_nrows], a[i] == bm[i], base_name="m_c")
        @constraint(master, convexity, 0 <= 1) # fix later

        @objective(master, Min, sum(M*a))

        lambda = [] # JuMP variable reference 

        println("master built")
    else
        master = copy(master)
        set_optimizer(master, GLPK.Optimizer)
        # retrieving references
        # m_c = master[:m_c]
        convexity = master[:convexity]
        a = master[:a]
    
        lambda = [variable_by_name(master, "λ_$(i)") for i in 1:length(S)] # JuMP variable reference 
        # println(lambda)
        println("master copied")
    end

    # if !(isnothing(solved))
    #     for s in solved
    #         println(s[1])
    #         println(s[2][6])
    #     end
    # end
    m_c = get_Ax_constraints(master, A_nrows)
    
    if verbose >= 1
        println(master)
        println("S = $(S)")
    end
    optimize!(master)
    m_obj = objective_value(master)
    
    pi_bar = dual.(m_c)
    pi_bar = reshape(pi_bar, 1, A_nrows)
    nu_bar = dual.(convexity)

    ### build the price
    # if isnothing(price)
    #     price = Model(GLPK.Optimizer)
    #     set_silent(price)

    #     @variable(price, x[1:x_amount] >= 0)

    #     @constraint(
    #         price, 
    #         p_c[i in 1:P_nrows],  
    #         sum([P[i, j]*x[j] for j in 1:x_amount]) <= bp[i] 
    #     )
    # else
    #     price = copy(price)
    #     set_optimizer(price, GLPK.Optimizer)

    #     # retrieving references
    #     p_c = price[:p_c]
    #     x = price[:x]
    # end

    cga_done = false
    # k = 1
    # while !(cga_done)
    for k in 1:20
        
        if verbose >= 1
            println("k = $(k)")        
            println("m_obj = $(m_obj)\n")
        end

        ### build the price
        price = Model(GLPK.Optimizer)
        set_silent(price)

        @variable(price, x[1:x_amount] >= 0)
        @constraint(
            price, 
            p_c[i in 1:P_nrows],  
            sum([P[i, j]*x[j] for j in 1:x_amount]) <= bp[i] 
        )

        m_c = get_Ax_constraints(master, A_nrows)
        pi_bar = dual.(m_c)
        pi_bar = reshape(pi_bar, 1, A_nrows)
        nu_bar = dual.(convexity)
        
        # println(c)
        # println(pi_bar)
        # println(A)
        # println(x)
        @objective(price, Min, sum((c-pi_bar*A)*x) - nu_bar )
                
        if verbose >= 1
            println("current price:")
            println(price)
        end
        optimize!(price)

        # println(termination_status(price))

        p_obj = objective_value(price)
        if verbose >= 1
            println("p_obj = $(p_obj)")
        end

        if p_obj < -1e-5
            # extract and save relevant values
            q = value.(price[:x])

            if verbose >= 1
                println("adding x = $(q)\n")
            end

            push!(S, q)
            push!(Aq, A*q)
            push!(cq, dot(c,q))

            ### generating the new column
            
            # new variable number
            # var_n = length(S)

            # create variable and add to master
            new_lambda = @variable(master, lower_bound=0, base_name="λ_$(length(S))")
            push!(lambda, new_lambda)
            # println(lambda)

            # set variable cost on master
            set_objective_function(master, objective_function(master) + cq[end]*new_lambda)

            # set coefficient of new variable on master's original constraints
            for i in 1:A_nrows
                set_normalized_coefficient(m_c[i], new_lambda, Aq[end][i])
            end

            # set coefficient of new variable on convexity constraint
            set_normalized_coefficient(convexity, new_lambda, 1)

                
            if verbose >= 1
                println("updated master:")
                println(master)
            end
            optimize!(master)

            m_obj = objective_value(master)

            # pi_bar = dual.([m_c[i] for i in 1:A_nrows])
            # pi_bar = dual.(m_c)
            # nu_bar = dual.(convexity)

            # pi_bar = dual.(master[:m_c])
            m_c = get_Ax_constraints(master, A_nrows)
            pi_bar = dual.(m_c)
            pi_bar = reshape(pi_bar, 1, A_nrows)
            nu_bar = dual.(convexity)

        else
            cga_done = true
        end
        if cga_done
            break
        end
    end

    if !(cga_done) 
        println("CGA failed to converge")
    end
    
    lambda_bar = value.(lambda)
    solution = sum([lambda_bar[j]*S[j] for j in 1:length(S)])

    return solution, m_obj, S, master, price

end

# "Solves with Column Generation Algorithm"
# function cga_old(A, b, c, sl=1, M=9e6)

#     constraints_amount = size(A)[1]
#     x_amount = length(c)
#     S = []
#     Aq = []
#     cq = [] 
#     lambda = [] # JuMP variable reference 

#     ### BUILDING THE MASTER

#     master = Model(GLPK.Optimizer)

#     # artificial variable
#     @variable(master, a[1:sl] >= 0)

#     # constraints 1 and convexity with artificial variables only
#     @constraint(master, m_c[i in 1:sl], a[i] == b[i])
#     @constraint(master, convexity, 0 <= 1)

#     @objective(master, Min, sum(M*a))

#     println(master)
#     optimize!(master)
#     m_obj = objective_value(master)

#     pi_bar = dual.([m_c[i] for i in 1:sl])
#     nu_bar = dual.(convexity)

#     ### BUILDING THE PRICE

#     price = Model(GLPK.Optimizer)
#     set_silent(price)

#     @variable(price, x[1:x_amount] >= 0)

#     @constraint(
#         price, 
#         p_c[i in sl+1:constraints_amount],  
#         sum([A[i, j]*x[j] for j in 1:x_amount]) <= b[i] 
#     )

#     for k in 1:1e3
#         println("k = $(k)")
#         println("m_obj = $(m_obj)")

#         @objective(price, Min, sum((c-pi_bar*A[1:sl, :])*x) - nu_bar )

#         println(price)
#         optimize!(price)

#         p_obj = objective_value(price)

#         if p_obj < 0
#             # extract and save relevant values
#             x_bar = value.(x)
#             Ax_bar = A*x_bar
#             cx_bar = dot(c,x_bar)
            
#             println("adding x = $(x_bar)\n")

#             push!(S, x_bar)
#             push!(Aq, Ax_bar)
#             push!(cq, cx_bar)

#             # new variable number
#             var_n = length(S)

#             # create variable and add to master
#             new_lambda = @variable(master, lower_bound=0, base_name="λ_$(length(S))")
#             push!(lambda, new_lambda)

#             # set variable cost on master
#             set_objective_function(master, objective_function(master) + dot(c,x_bar)*new_lambda )

#             # set coefficient of new variable on master's original constraints  
#             for i in 1:sl
#                 set_normalized_coefficient(m_c[i], new_lambda, Ax_bar[i])
#             end

#             # set coefficient of new variable on convexity constraint 
#             set_normalized_coefficient(convexity, new_lambda, 1)

#             println(master)
#             optimize!(master)

#             m_obj = objective_value(master)

#             pi_bar = dual.([m_c[i] for i in 1:sl])
#             nu_bar = dual.(convexity)

#         else
#             break
#         end
#     end

#     println("done")

#     lambda_bar = value.(lambda)
#     solution = sum([lambda_bar[j]*S[j] for j in 1:length(S)])

#     return solution, m_obj, master, price
# end


# c = [4 9 7 3]

# A = [
#     2 3 -1 3;
#     1 4 0 0;
#     2 3 0 0;
#     0 0 4 3;
#     0 0 2 5
# ]

# b = [13, 5, 7, 20, 25]

# sl = 1

# M=9e6

# sol, obj_v, model = rsa(A, b, -c, false)
# println(sol)

# sol2, obj_v2, m = rsa(A, b, -c, true)
# println(sol2)

# best_solution, upper_bound = bba(model)