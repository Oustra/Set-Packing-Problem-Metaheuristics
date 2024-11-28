
# --------------------------------------------------------------------------- #

# Improvement Metaheuristics (Improve the result of the construction Metaheuristics)

# Simple Descent

function local_search_0_1(x0, C, A, conflicts, non_selected_indices, temp_constraints)
    m, n = size(A)
    best_x = copy(x0)
    best_z = dot(C, x0)
    current_constraints = sum(A[:, j] .* best_x[j] for j in 1:n)
    found_improvement = false

    while true
        improved = false
        n_non_selected = 0

        # Non-selected indices
        for i in 1:n
            if best_x[i] == 0
                n_non_selected += 1
                non_selected_indices[n_non_selected] = i
            end
        end

        for idx in 1:n_non_selected
            j = non_selected_indices[idx]
            feasible = true
            
            copyto!(temp_constraints, current_constraints)
            
            for c in conflicts[j]
                temp_constraints[c] += 1
                if temp_constraints[c] > 1
                    feasible = false
                    break
                end
            end

            # Updating the solution
            if feasible
                candidate_z = best_z + C[j]
                if candidate_z > best_z
                    best_z = candidate_z
                    best_x[j] = 1
                    copyto!(current_constraints, temp_constraints)
                    improved = true
                    found_improvement = true
                    break 
                end
            end
        end

        if !improved 
            break 
        end
    end

    return best_x, best_z, found_improvement
end

# --------------------------------------------------------------------------- #

function local_search_1_1(x0, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
    m, n = size(A)
    best_x = copy(x0)
    best_z = dot(C, x0)
    current_constraints = sum(A[:, j] .* best_x[j] for j in 1:n)
    found_improvement = false

    while true
        improved = false
        n_selected = 0
        n_non_selected = 0

        # Non-selected indices
        for i in 1:n
            if best_x[i] == 1
                n_selected += 1
                selected_indices[n_selected] = i
            else
                n_non_selected += 1
                non_selected_indices[n_non_selected] = i
            end
        end

        for idx in 1:n_selected
            i = selected_indices[idx]
            copyto!(candidate_constraints, current_constraints)

            # Updating constraints
            for c in conflicts[i]
                candidate_constraints[c] -= 1
            end
            
            for jdx in 1:n_non_selected
                j = non_selected_indices[jdx]
                feasible = true

                copyto!(temp_constraints, candidate_constraints)

                for c in conflicts[j]
                    temp_constraints[c] += 1
                    if temp_constraints[c] > 1
                        feasible = false
                        break
                    end
                end

                # Update solution
                if feasible
                    candidate_z = best_z + C[j] - C[i]
                    if candidate_z > best_z
                        best_z = candidate_z
                        best_x[i] = 0
                        best_x[j] = 1
                        copyto!(current_constraints, temp_constraints)
                        improved = true
                        found_improvement = true

                        break 
                    end
                end
            end
            
            if improved 
                break 
            end 
        end

        if !improved 
            break 
        end 
    end

    return best_x, best_z, found_improvement
end

# --------------------------------------------------------------------------- #

function local_search_2_1(x0, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
    m, n = size(A)
    best_x = copy(x0)
    best_z = dot(C, x0)
    current_constraints = sum(A[:, j] .* best_x[j] for j in 1:n)
    found_improvement = false

    while true
        improved = false
        n_selected = 0
        n_non_selected = 0

        for i in 1:n
            if best_x[i] == 1
                selected_indices[n_selected += 1] = i
            else
                non_selected_indices[n_non_selected += 1] = i
            end
        end

        for i in 1:(n_selected - 1)
            i_idx = selected_indices[i]
            for j in (i + 1):n_selected
                j_idx = selected_indices[j]
                
                copyto!(candidate_constraints, current_constraints)
                for c in conflicts[i_idx]
                    candidate_constraints[c] -= 1
                end
                for c in conflicts[j_idx]
                    candidate_constraints[c] -= 1
                end

                # Check feasibility with non-selected items
                for k in 1:n_non_selected
                    k_idx = non_selected_indices[k]
                    feasible = true

                    copyto!(temp_constraints, candidate_constraints)
                    for c in conflicts[k_idx]
                        if (temp_constraints[c] += 1) > 1
                            feasible = false
                            break
                        end
                    end

                    # Update solution
                    if feasible
                        candidate_z = best_z + C[k_idx] - (C[i_idx] + C[j_idx])
                        if candidate_z > best_z
                            best_z = candidate_z
                            best_x[i_idx] = 0
                            best_x[j_idx] = 0
                            best_x[k_idx] = 1
                            copyto!(current_constraints, temp_constraints)
                            found_improvement = true
                            improved = true

                            break
                        end
                    end
                end

                if improved break end
            end
            if improved break end
        end

        if !improved break end
    end

    return best_x, best_z, found_improvement
end
       
# =========================================================================== #

function local_search_kp(x0, C, A)
    current_x = copy(x0)
    current_z = dot(C, current_x)
    
    # Pre-computing conflicts
    conflicts = [findall(view(A, :, i) .== 1) for i in 1:size(A, 2)]

    # Pre-allocating arrays
    m, n = size(A)
    selected_indices = Vector{Int}(undef, n)
    non_selected_indices = Vector{Int}(undef, n)
    candidate_constraints = similar(view(A, :, 1))
    temp_constraints = similar(candidate_constraints)

    best_x2, best_z2, found_improvement2 = local_search_2_1(current_x, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
    println("2-1 > ", best_z2)
    
    best_x1, best_z1, found_improvement1 = local_search_1_1(best_x2, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
    println("1-1 > ", best_z1)
        
    best_x, best_z, found_improvement = local_search_0_1(best_x1, C, A, conflicts, non_selected_indices, temp_constraints)
    println("0-1 > ", best_z)

    return best_x, best_z
end
