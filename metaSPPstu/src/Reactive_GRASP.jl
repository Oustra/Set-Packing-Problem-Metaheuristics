# =========================================================================== #
using Distributions
# =========================================================================== #

# Fonctions pour résoudre le problème avec des heuristiques de construction

function reactive_GRASP_construction(C, A, alphas, num_iterations=100, N_alpha=10)
    m, n = size(A)
    T = sum(A, dims=1)[:]
    prio = C ./ T
    prio[T .== 0] .= 0
    zbetter = 0 
    xbetter = zeros(Int, n)

    # Variable for Local Search
    conflicts = [findall(view(A, :, i) .== 1) for i in 1:size(A, 2)]
    selected_indices = Vector{Int}(undef, n)
    non_selected_indices = Vector{Int}(undef, n)
    candidate_constraints = similar(view(A, :, 1))
    temp_constraints = similar(candidate_constraints)
    
    # Arrays to store values for simulation
    zconstruction = zeros(Float64, num_iterations)
    zamelioration = zeros(Float64, num_iterations)
    zbest = zeros(Float64, num_iterations)

    # Initialize solution quality tracking for each alpha
    pk = fill(1.0 / length(alphas), length(alphas))  # p_k
    qk = zeros(Float64, length(alphas))  # q_k intermediate weights
    alpha_solutions = zeros(Float64, length(alphas))  # Sum of solutions for each alpha
    alpha_counts = zeros(Int, length(alphas))  # Counts for averaging

    for iter in 1:num_iterations
        # Select alpha based on probability distribution
        alpha_index = rand(Categorical(pk))
        alpha = alphas[alpha_index]
        
        # Initialize the solution for this iteration
        x_init = zeros(Int, n)
        current_constraints = zeros(Float64, m)
        
        while true
            # Computing the RCL based on the top `alpha` percent elements
            candidates = []
            for j in 1:n
                if x_init[j] == 0 && all(current_constraints .+ A[:, j] .<= 1)
                    push!(candidates, (j, prio[j]))
                end
            end

            if isempty(candidates)
                break
            end

            # Sorting candidates by priority and form the RCL
            sorted_candidates = sort(candidates, by=x->x[2], rev=true)
            rcl_size = max(1, Int(round(alpha * length(sorted_candidates))))
            rcl = sorted_candidates[1:rcl_size]

            # Randomly select from the RCL
            selected = rand(rcl)[1]
            # Updating the solution and constraints
            x_init[selected] = 1
            current_constraints += A[:, selected]
        end

        z_init = dot(C, x_init)
        zconstruction[iter] = z_init

        # Updating the cumulative solution quality for this alpha
        alpha_solutions[alpha_index] += z_init
        alpha_counts[alpha_index] += 1
        
        # Improving with local search
        #x_Ls, z_Ls = local_search_kp(x_init, C, A) # Local search 2-1 1-1 0-1 on sequence
        x_Ls, z_Ls = local_search_1_1(x_init, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
        zamelioration[iter] = z_Ls
       
        # Update the best solution
        if zbetter < z_Ls
            xbetter = x_Ls
            zbetter = z_Ls
        end
        println("-> zbest = ",zbetter)
        zbest[iter] = zbetter

        # Best and Worst solution Found
        zBest = zbetter
        zWorst = minimum(zconstruction)

        # Update q_k and p_k every N_alpha iterations
        if iter % N_alpha == 0
            
            for k in 1:length(alphas)
                if alpha_counts[k] > 0
                    zAvg_k = alpha_solutions[k] / alpha_counts[k]
                    qk[k] = (zAvg_k - zWorst) / (zBest - zWorst)
                else
                    qk[k] = 0
                end
            end
            
            # Normalizing the probabilities
            total_weight = sum(qk)
            if total_weight > 0
                pk .= qk / total_weight
            else
                pk .= 1.0 / length(alphas)
            end
        end
    end
    
    return zconstruction, zamelioration, zbest, pk
end