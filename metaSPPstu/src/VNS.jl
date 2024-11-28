# =========================================================================== #
include("Local_Search.jl")
# =========================================================================== #

# Improvement Metaheuristics (Improve the result of the construction Metaheuristics)
       
function vns(x0, C, A, max_iterations=100)
    current_x = copy(x0)
    current_z = dot(C, current_x)

    # Pre-compute conflicts
    conflicts = [findall(view(A, :, i) .== 1) for i in 1:size(A, 2)]

    # Pre-allocate arrays
    m, n = size(A)
    selected_indices = Vector{Int}(undef, n)
    non_selected_indices = Vector{Int}(undef, n)
    candidate_constraints = similar(view(A, :, 1))
    temp_constraints = similar(candidate_constraints)

    # Stopping condition
    iteration = 0
    while iteration < max_iterations
        # Missing chacking phase...
        iteration += 1
        found_improvement = false

        for neighborhood in [:local_search_2_1, :local_search_1_1, :local_search_0_1]
            if neighborhood == :local_search_2_1
                best_x, best_z, improved = local_search_2_1(current_x, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
            elseif neighborhood == :local_search_1_1
                best_x, best_z, improved = local_search_1_1(current_x, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
            elseif neighborhood == :local_search_0_1
                best_x, best_z, improved = local_search_0_1(current_x, C, A, conflicts, non_selected_indices, temp_constraints)
            end

            if improved && best_z > current_z
                current_x = copy(best_x)
                current_z = best_z
                found_improvement = true
                println("Improved in $neighborhood > $current_z")
                break
            end
        end

        if !found_improvement
            println("No improvement found in iteration $iteration. Stopping.")
            break
        end
    end

    return current_x, current_z
end

