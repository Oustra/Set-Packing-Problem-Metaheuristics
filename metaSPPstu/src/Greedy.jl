
# --------------------------------------------------------------------------- #

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
  

# --------------------------------------------------------------------------- #

# Greedy Metaheuristic for ACO 
function greedy_construction_phi(φ, C, A)
    m, n = size(A)
    T = sum(A, dims=1)[:]

    prio = (φ .* C) ./ T
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