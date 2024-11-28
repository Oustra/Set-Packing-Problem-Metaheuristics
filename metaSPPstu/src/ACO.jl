# =========================================================================== #
include("Greedy.jl")
include("Local_Search.jl")
using Distributions
# =========================================================================== #

function ant_colony_optimization(C, A, maxIter, maxAnt, rhoE, rhoD, phiInit)
    n_items = size(A, 2)                             # Number of items
    φ = fill(phiInit, n_items)                       # Initial pheromone levels
    best_x = zeros(Int, n_items)                     # Best solution found
    best_z = typemin(Int64)                                    # Best objective value
    zbest = zeros(Int64, maxIter * maxAnt)         # Global best values for each ant
    ziteration = zeros(Int64, maxIter * maxAnt)    # Iteration best values for each ant
    all_solutions = zeros(Int64, maxIter * maxAnt) # All solutions for plotting
    count = 0
      
    # Precomputed variables for Local-Search 
    conflicts = [findall(view(A, :, i) .== 1) for i in 1:size(A, 2)]
    selected_indices = Vector{Int}(undef, n_items)
    non_selected_indices = Vector{Int}(undef, n_items)
    candidate_constraints = similar(view(A, :, 1))
    temp_constraints = similar(candidate_constraints)

    for iter in 1:maxIter
        @show iter
        iter_best_x = zeros(Int, n_items)
        iter_best_z = -Inf
        iter_solutions = Float64[]

        for ant in 1:maxAnt
            # Decide between greedy or probabilistic construction
            if rand() < exploitation_probability(iter, maxIter)
                #solution, objective = greedy_construction(C, A)         # For static Construction (SpeedUp)
                solution, objective = greedy_construction_phi(φ, C, A)   # For dynamic Construction based on φ (Slow)
            else
                solution, objective = elaborate_solution_probabilistic(φ, C, A)
            end

            # Apply local search for refinement
            local_solution, local_objective = local_search_1_1(solution, C, A, conflicts, selected_indices, non_selected_indices, candidate_constraints, temp_constraints)
            
            count += 1
            all_solutions[count] = local_objective

            # Update iteration's best solution
            if local_objective > iter_best_z
                iter_best_x = local_solution
                iter_best_z = local_objective
            end
                  
            if iter_best_z > best_z
                #best_x = iter_best_x
                best_z = iter_best_z
            end
            
            if best_z > 0 zbest[count] = best_z end
            ziteration[count] = iter_best_z

        end

        # Update global best solution
        if iter_best_z > best_z
            #best_x = iter_best_x
            best_z = iter_best_z
        end

        # Update pheromones
        φ = manage_pheromones(φ, iter_best_x, rhoE, rhoD, iter, maxIter)

    end

    return zbest, ziteration, all_solutions
end

# --------------------------------------------------------------------------- #

function elaborate_solution_probabilistic(φ, C, A)
    n_items = length(φ)
    solution = zeros(Int, n_items)
    remaining_items = collect(1:n_items)
    max_iterations = 1000 
    iteration_count = 0

    while !isempty(remaining_items) && iteration_count < max_iterations
        iteration_count += 1
        # Compute selection probabilities using pheromone levels
        probabilities = φ[remaining_items] / sum(φ[remaining_items])
        #@show remaining_items
        # Use Distributions package for weighted sampling
        dist = Categorical(probabilities)  # Define the distribution
        selected_item = remaining_items[rand(dist)]  # Sample from the distribution
        solution[selected_item] = 1

        # Remove incompatible items based on constraints
        remaining_items = filter(i -> all(A[:, selected_item] + A[:, i] .≤ 1), remaining_items)
    end

    # Compute objective value
    objective = dot(C, solution)
    return solution, objective
end

# --------------------------------------------------------------------------- #

function manage_pheromones(φ, best_solution, rhoE, rhoD, iter, maxIter)
    # Evaporate pheromones
    φ = φ .* rhoE

    # Reinforce pheromones on the items in the best solution
    for i in findall(x -> x == 1, best_solution)
        φ[i] += rhoD
    end

    # Optional: Add perturbation to avoid stagnation
    if iter / maxIter > 0.8
        φ .+= rand(length(φ)) .* 0.1
    end

    return φ
end

# --------------------------------------------------------------------------- #

function exploitation_probability(iter, maxIter)
    # Dynamic probability: starts with exploration and shifts to exploitation
    return log10(iter + 1) / log10(maxIter + 1)
end

# --------------------------------------------------------------------------- #

