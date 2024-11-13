using CategoricalArrays
using Distributions
using Plots
# --------------------------------------------------------------------------- #

# Fonctions pour résoudre le problème avec des heuristiques de construction

# Greedy Metaheuristic
function greedy_construction(C, A)
    m, n = size(A)
    T = sum(A, dims=1)[:]  

    prio = C ./ T
    prio[T .== 0] .= 0 
    sorted_indices = sortperm(prio, rev=true)
    
    x0 = zeros(Int, n)
    current_constraints = zeros(Float64, m)
    
    for j in sorted_indices
        if all(current_constraints .+ A[:, j] .<= 1)
            x0[j] = 1
            current_constraints += A[:, j]
        end
    end
  
    z = dot(C, x0)
    
    return x0, z
end
  
# The GRASP metaheuristic
function GRASP_construction(C, A, alpha=0.2, num_iterations=25)
    m, n = size(A)
    T = sum(A, dims=1)[:]
    prio = C ./ T
    prio[T .== 0] .= 0
    z_best = 0 
    x_best = zeros(Int, n)
    z_best_values = Float64[] 
    z_Ls_values = Float64[]  
    z_init_values = Float64[]

    for iter in 1:num_iterations
        # Initialize the solution
        x_init = zeros(Int, n)
        current_constraints = zeros(Float64, m)

        while true
            # Compute RCL based on the top `alpha` percent elements
            candidates = []
            for j in 1:n
                if x_init[j] == 0 && all(current_constraints .+ A[:, j] .<= 1)
                    push!(candidates, (j, prio[j]))
                end
            end

            # If there are no feasible candidates, stop
            if isempty(candidates)
                break
            end

            # Sort candidates by priority and form the RCL
            sorted_candidates = sort(candidates, by=x->x[2], rev=true)
            rcl_size = max(1, Int(round(alpha * length(sorted_candidates))))
            rcl = sorted_candidates[1:rcl_size]
            
            # Randomly select from the RCL
            selected = rand(rcl)[1]
            # Update the solution and constraints
            x_init[selected] = 1
            current_constraints += A[:, selected]
            #z0 = dot(C, x_init)
        end
        z_init = dot(C, x_init)
        push!(z_init_values, z_init)

        # Improving with local search
        x_Ls, z_Ls = local_search_kp(x_init, C, A)
        push!(z_Ls_values, z_Ls)

        if z_best < z_Ls
            x_best = x_Ls
            z_best = z_Ls
        end
        push!(z_best_values, z_best)

    end

    return x_best, z_best, z_init_values, z_Ls_values, z_best_values
end

# The Reactive GRASP metaheuristic
function reactive_GRASP_construction(C, A, alphas=[0.1, 0.2, 0.3, 0.4, 0.5], num_iterations=100)
    m, n = size(A)
    T = sum(A, dims=1)[:]
    prio = C ./ T
    prio[T .== 0] .= 0
    z_best = -Inf
    x_best = zeros(Int, n)
    z_best_values = Float64[] 
    z_Ls_values = Float64[]  
    z_init_values = Float64[]

    # Initialize solution quality tracking for each alpha
    alpha_probabilities = fill(1.0 / length(alphas), length(alphas))
    alpha_weights = zeros(Float64, length(alphas))

    for iter in 1:num_iterations
        # Select alpha based on probability distribution
        alpha_index = rand(Categorical(alpha_probabilities))
        alpha = alphas[alpha_index]
        
        # Initialize the solution for this iteration
        x_init = zeros(Int, n)
        current_constraints = zeros(Float64, m)
        
        while true
            # Compute the RCL based on the top `alpha` percent elements
            candidates = []
            for j in 1:n
                if x_init[j] == 0 && all(current_constraints .+ A[:, j] .<= 1)
                    push!(candidates, (j, prio[j]))
                end
            end

            # If there are no feasible candidates, stop
            if isempty(candidates)
                break
            end

            # Sort candidates by priority and form the RCL
            sorted_candidates = sort(candidates, by=x->x[2], rev=true)
            rcl_size = max(1, Int(round(alpha * length(sorted_candidates))))
            rcl = sorted_candidates[1:rcl_size]

            # Randomly select from the RCL
            selected = rand(rcl)[1]
            # Update the solution and constraints
            x_init[selected] = 1
            current_constraints += A[:, selected]
        end

        # Evaluate the current solution
        z_init = dot(C, x_init)
        # Update the alpha weights based on solution quality
        alpha_weights[alpha_index] += z_init
        push!(z_init_values, z_init)

        # Improving with local search
        x_Ls, z_Ls = local_search_kp(x_init, C, A)
        push!(z_Ls_values, z_Ls)

        if z_best < z_Ls
            x_best = x_Ls
            z_best = z_Ls
        end
        push!(z_best_values, z_best)

    end

    # Update alpha probabilities based on the accumulated weights
    total_weight = sum(alpha_weights)
    alpha_probabilities .= alpha_weights / total_weight
    
    return x_best, z_best, alpha_probabilities,  z_init_values, z_Ls_values, z_best_values
end

# The Reactive GRASP metaheuristic
function reactive_GRASP_construction2(C, A, alphas=[0.1, 0.2, 0.3, 0.4, 0.5], num_iterations=1000)
    m, n = size(A)
    T = sum(A, dims=1)[:]
    prio = C ./ T
    prio[T .== 0] .= 0

    # Initialize solution quality tracking for each alpha
    alpha_probabilities = fill(1.0 / length(alphas), length(alphas))
    alpha_weights = zeros(Float64, length(alphas))
    best_solution = nothing
    best_z = -Inf

    z_values = Float64[]
    iterations = Int[] 
    iteration_count = 0

    for iter in 1:num_iterations
        # Select alpha based on probability distribution
        alpha_index = rand(Categorical(alpha_probabilities))
        alpha = alphas[alpha_index]
        
        # Initialize the solution for this iteration
        x0 = zeros(Int, n)
        current_constraints = zeros(Float64, m)
        
        while true
            # Compute the RCL based on the top `alpha` percent elements
            candidates = []
            for j in 1:n
                if x0[j] == 0 && all(current_constraints .+ A[:, j] .<= 1)
                    push!(candidates, (j, prio[j]))
                end
            end

            # If there are no feasible candidates, stop
            if isempty(candidates)
                break
            end

            # Sort candidates by priority and form the RCL
            sorted_candidates = sort(candidates, by=x->x[2], rev=true)
            rcl_size = max(1, Int(round(alpha * length(sorted_candidates))))
            rcl = sorted_candidates[1:rcl_size]

            # Randomly select from the RCL
            selected = rand(rcl)[1]
            # Update the solution and constraints
            x0[selected] = 1
            current_constraints += A[:, selected]
        end
        
        # Evaluate the current solution
        z = dot(C, x0)
        
        # Update the alpha weights based on solution quality
        alpha_weights[alpha_index] += z
        if z > best_z   
            z_values = push!(z_values, z)
            iterations = push!(iterations, iter)
            best_z = z
            best_solution = copy(x0)
        end
    end

    # Update alpha probabilities based on the accumulated weights
    total_weight = sum(alpha_weights)
    alpha_probabilities .= alpha_weights / total_weight
    
    return best_solution, best_z, alpha_probabilities, z_values , iterations
end



# =========================================================================== #

# Improvement Metaheuristics (Improve the result of the construction Metaheuristics)

# Simple Descent

function local_search_0_1(x0, C, A, conflicts)
    m, n = size(A)
    best_x = copy(x0)
    best_z = dot(C, x0)
    current_constraints = sum(A[:, j] .* best_x[j] for j in 1:n)
    found_improvement = false

    while true
        improved = false
        non_selected_indices = findall(x -> x == 0, best_x)

        for j in non_selected_indices
            feasible = true
            temp_constraints = copy(current_constraints)
            
            # Check if adding item `j` would violate constraints
            for c in conflicts[j]
                temp_constraints[c] += 1
                if temp_constraints[c] > 1
                    feasible = false
                    break
                end
            end

            # Update the solution if feasible and improves objective
            if feasible
                candidate_z = best_z + C[j]
                if candidate_z > best_z
                    best_z = candidate_z
                    best_x[j] = 1
                    current_constraints = temp_constraints
                    improved = true
                    found_improvement = true
                    break
                end
            end
        end

        # Stop if no improvements were made in this iteration
        if !improved
            break
        end
    end

    return best_x, best_z, found_improvement
end

function local_search_1_1(x0, C, A, conflicts)
    m, n = size(A)
    best_x = copy(x0)
    best_z = dot(C, x0)
    current_constraints = sum(A[:, j] .* best_x[j] for j in 1:n)
    found_improvement = false

    while true
        improved = false
        selected_indices = findall(x -> x == 1, best_x)
        non_selected_indices = findall(x -> x == 0, best_x)

        for i in selected_indices
            candidate_constraints = copy(current_constraints)
            
            # Update constraints by removing element i
            for c in conflicts[i]
                candidate_constraints[c] -= 1
            end
            
            for j in non_selected_indices
                feasible = true
                temp_constraints = copy(candidate_constraints)
                
                # Check if adding element j is feasible
                for c in conflicts[j]
                    temp_constraints[c] += 1
                    if temp_constraints[c] > 1
                        feasible = false
                        break
                    end
                end

                # If feasible and improves objective, update solution
                if feasible
                    candidate_z = best_z + C[j] - C[i]
                    if candidate_z > best_z
                        best_z = candidate_z
                        best_x[i] = 0
                        best_x[j] = 1
                        current_constraints = temp_constraints
                        improved = true
                        found_improvement = true
                        break
                    end
                end
            end
            
            # Stop searching if improvement was found
            if improved
                break
            end
        end

        # Terminate if no improvement was made
        if !improved
            break
        end
    end

    return best_x, best_z, found_improvement
end

function local_search_2_1(x0, C, A, conflicts)
    m, n = size(A)
    best_x = copy(x0)
    best_z = dot(C, x0)
    current_constraints = sum(A[:, j] .* best_x[j] for j in 1:n)
    found_improvement = false

    while true
        improved = false
        selected_indices = findall(x -> x == 1, best_x)
        non_selected_indices = findall(x -> x == 0, best_x)

        for i in 1:(length(selected_indices) - 1)
            for j in (i + 1):length(selected_indices)
                candidate_constraints = copy(current_constraints)

                # Update constraints by removing items i and j
                for c in conflicts[selected_indices[i]]
                    candidate_constraints[c] -= 1
                end
                for c in conflicts[selected_indices[j]]
                    candidate_constraints[c] -= 1
                end

                for k in non_selected_indices
                    temp_constraints = copy(candidate_constraints)
                    feasible = true

                    # Check feasibility of adding item k
                    for c in conflicts[k]
                        temp_constraints[c] += 1
                        if temp_constraints[c] > 1
                            feasible = false
                            break
                        end
                    end

                    # Update solution if feasible and improves objective
                    if feasible
                        candidate_z = best_z + C[k] - (C[selected_indices[i]] + C[selected_indices[j]])
                        if candidate_z > best_z
                            best_z = candidate_z
                            best_x[selected_indices[i]] = 0
                            best_x[selected_indices[j]] = 0
                            best_x[k] = 1
                            current_constraints = temp_constraints
                            found_improvement = true
                            improved = true
                            break
                        end
                    end
                end

                # Exit outer loop if improvement was found
                if improved
                    break
                end
            end
            if improved
                break
            end
        end

        # Stop if no improvement was made
        if !improved
            break
        end
    end

    return best_x, best_z, found_improvement
end

# =========================================================================== #

function local_search_kp(x0, C, A)
    current_x = copy(x0)
    current_z = dot(C, current_x)
    conflicts = Dict(i => findall(A[:, i] .== 1) for i in 1:size(A, 2))

    best_x2, best_z2, found_improvement2 = local_search_2_1(current_x, C, A, conflicts)
    println("2-1 > ",best_z2)
    best_x1, best_z1, found_improvement1 = local_search_1_1(best_x2, C, A, conflicts)
    println("1-1 > ",best_z1)
    best_x, best_z, found_improvement = local_search_0_1(best_x1, C, A, conflicts)
    println("0-1 > ",best_z)
    println("--------------------------------------------------")

    return best_x, best_z
end

# =========================================================================== #

function tabu_search(x0, C, A; max_iters=100, tabu_tenure=7)
    m, n = size(A)
    current_x = copy(x0)
    best_x = copy(x0)
    current_z = dot(C, current_x)
    best_z = current_z

    # Initialize constraints based on x0
    current_constraints = sum(A[:, j] .* current_x[j] for j in 1:n)
    
    tabu_list = []

    conflicts = Dict(i => findall(A[:, i] .== 1) for i in 1:n)

    z_values = [current_z]
    iterations = Int[] 
    iteration_count = 0

    for iter in 1:max_iters
        println("> Iteration : ", iter)
        neighbors = []
        iteration_count += 1

        for j in 1:n
            # Create a copy of constraints for feasibility checking
            temp_constraints = copy(current_constraints)
            feasible = true
            
            # Check feasibility if we flip x[j]
            for c in conflicts[j]
                temp_constraints[c] += ifelse(current_x[j] == 1, -1, 1)
                if temp_constraints[c] > 1
                    feasible = false
                    break
                end
            end

            if feasible
                # Calculate new objective incrementally
                neighbor_z = current_z + (1 - 2 * current_x[j]) * C[j]
                push!(neighbors, (neighbor_z, j))
            end
        end

        # Select best non-tabu neighbor
        best_neighbor_z, best_move = -Inf, -1
        for (neighbor_z, move) in neighbors
            if !in(move, tabu_list) && neighbor_z > best_neighbor_z
                best_neighbor_z, best_move = neighbor_z, move
            end
        end
        
        # Update solution if a valid move is found
        if best_move != -1
            # Flip the selected variable
            current_x[best_move] = 1 - current_x[best_move]
            current_z = best_neighbor_z
            current_constraints = sum(A[:, j] .* current_x[j] for j in 1:n)  # or update incrementally if feasible
            
            # Update global best if necessary
            println("- Neighbor solution : ", current_z)
            if current_z > best_z
                best_x, best_z = copy(current_x), current_z
            end

            # Update tabu list
            push!(tabu_list, best_move)
            if length(tabu_list) > tabu_tenure
                popfirst!(tabu_list)
            end
        else
            println("- No non-tabu feasible neighbors found.")
            break
        end

        println("- Best solution value : ", best_z)
        println("---------------------------------------")
        # Record the best objective value at this iteration for plotting
        z_values = push!(z_values, best_z)
        iterations = push!(iterations, iteration_count)
    end
    
    return best_x, best_z, z_values, iterations
end

# =========================================================================== #
