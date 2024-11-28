
# --------------------------------------------------------------------------- #

# Functions to solve the problem with construction heuristics
  
# The GRASP metaheuristic

function GRASP_construction(C, A, alpha=0.2, num_iterations=25)
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
    
    # Initialize the solution
    x_init = zeros(Int, n)
    current_constraints = zeros(Float64, m)

    for iter in 1:num_iterations
        fill!(x_init, 0)
        fill!(current_constraints, 0.0)

        while true
            # Computing RCL based on the top `alpha` percent elements
            candidates = [(j, prio[j]) for j in 1:n if x_init[j] == 0 && all(current_constraints .+ A[:, j] .<= 1)]

            if isempty(candidates)
                break
            end

            # Sorting candidates by priority and form the RCL using a more efficient approach
            sorted_candidates = sort(candidates, by=x->x[2], rev=true)
            rcl_size = max(1, Int(round(alpha * length(sorted_candidates))))
            rcl = sorted_candidates[1:rcl_size]
            
            # Randomly select from the RCL using a more efficient method
            selected_index = rand(1:rcl_size)
            selected = rcl[selected_index][1]
            
            # Updating the solution and constraints
            x_init[selected] = 1
            current_constraints += A[:, selected]
        end
        
        z_init = dot(C, x_init)
        zconstruction[iter] = z_init
        
        # Improving with local search
        #x_Ls, z_Ls = local_search_kp(x_init, C, A)  # Local search 2-1 1-1 0-1 on sequence
        x_Ls, z_Ls = local_search_1_1(x_init, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
        zamelioration[iter] = z_Ls
        
        # Update the best solution
        if zbetter < z_Ls
            xbetter[:] = x_Ls[:]
            zbetter = z_Ls
        end
        
        println("-> zbest = ", zbetter)
        zbest[iter] = zbetter
    end

    # Return arrays compatible with simulation
    return zconstruction, zamelioration, zbest
end

