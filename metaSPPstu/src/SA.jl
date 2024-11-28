# =========================================================================== #
include("Local_Search.jl")
include("Greedy.jl")
# =========================================================================== #

# Simulated Annealing

function simulated_annealing(C, A, Tmax, Tmin, nbIterationSAF, p1, p2, C4)
    m, n = size(A)
    best_x = zeros(Int, n)
    best_z = 0
    current_x, _ = greedy_construction(C, A)
    current_z = dot(C, current_x)
    Temperature = Tmax

    conflicts = Dict(i => findall(A[:, i] .== 1) for i in 1:n)

    zinit = zeros(Int64, nbIterationSAF)
    zSals = zeros(Int64, nbIterationSAF)
    zbest = zeros(Int64, nbIterationSAF)

    # Initialize best solution with the initial solution
    if current_z > best_z
        best_z = current_z
        best_x = copy(current_x)
    end

    for iter in 1:nbIterationSAF
        zinit[iter] = current_z
        
        # Create a neighbor solution
        neighbor_x = copy(current_x)
        if rand() < p1
            neighbor_x = sa_local_search(neighbor_x, C, A, Temperature, conflicts, p2, C4)
        else
            neighbor_x, _ = local_search_kp(neighbor_x, C, A)
        end

        # Calculate the change in objective value
        neighbor_z = dot(C, neighbor_x)
        delta = neighbor_z - current_z

        # Acceptance criterion with normalized delta
        normalized_temp = C4 * Temperature * maximum(C)
        if delta > 0 || rand() < exp(delta / normalized_temp)
            current_x = copy(neighbor_x)
            current_z = neighbor_z

            # Update best solution if necessary
            if current_z > best_z
                best_z = current_z
                best_x = copy(current_x)
            end
        end

        zSals[iter] = current_z
        zbest[iter] = best_z

        # Cooling schedule
        Temperature *= C4

        # Dynamic restart mechanism
        if Temperature < Tmin
            Temperature = Tmax * 0.8  # Slightly reduced max temperature for diversification
            current_x, _ = greedy_construction(C, A)
            current_z = dot(C, current_x)
        end
    end

    return zinit, zSals, zbest
end

function sa_local_search(x, C, A, Temperature, conflicts, p2, C4)
    xtemp = copy(x)

    # Enhanced neighborhood exploration
    if rand() < p2
        xtemp = sa_exchange_1(xtemp, C, A, conflicts)
    else
        xtemp = sa_exchange_2(xtemp, C, A, conflicts)
    end

    # Improved acceptance criterion with normalized temperature
    delta = dot(C, xtemp) - dot(C, x)
    if delta > 0
        return xtemp
    else
        p = exp(delta / (C4 * Temperature * maximum(C)))  # Normalized temperature
        return rand() < p ? xtemp : x
    end
end

function sa_exchange_1(x, C, A, conflicts)
    n = length(x)
    # Filter only valid candidates with no conflicts
    feasible_candidates = filter(i -> i <= n && x[i] == 0 && 
        all(j -> j <= n && x[j] == 0, conflicts[i]), 1:n)

    if !isempty(feasible_candidates)
        selected = rand(feasible_candidates)
        x[selected] = 1
    end
    return x
end

function sa_exchange_2(x, C, A, conflicts)
    m, n = size(A)
    non_selected = findall(x .== 0)

    if length(non_selected) >= 2
        candidates = []
        candidate_values = []

        # Find pairs of non-conflicting elements and their combined values
        for i in 1:length(non_selected)
            for j in (i+1):length(non_selected)
                if all(A[:, non_selected[i]] .+ A[:, non_selected[j]] .<= 1)
                    push!(candidates, (non_selected[i], non_selected[j]))
                    push!(candidate_values, C[non_selected[i]] + C[non_selected[j]])
                end
            end
        end

        if !isempty(candidates)
            # Select pair with probability proportional to combined value
            probs = candidate_values / sum(candidate_values)
            selected_idx = findfirst(cumsum(probs) .>= rand())
            i, j = candidates[selected_idx]

            # Remove conflicting elements
            for k in findall(x .== 1)
                if any(A[:, k] .+ A[:, i] .> 1) || any(A[:, k] .+ A[:, j] .> 1)
                    x[k] = 0
                end
            end

            x[i] = 1
            x[j] = 1
        end
    end

    return x
end


  