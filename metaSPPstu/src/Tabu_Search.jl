
# --------------------------------------------------------------------------- #

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
            # Creating a copy of constraints for feasibility checking
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
                # Calculating new objective incrementally
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
        
        # Updating solution if a valid move is found
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
