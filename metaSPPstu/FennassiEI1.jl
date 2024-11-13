# =========================================================================== #
# Compliant julia 1.x

# Using the following packages
using JuMP, GLPK
using LinearAlgebra
using CategoricalArrays
using Plots

include("src/loadSPP.jl")
include("src/setSPP.jl")
include("src/getfname.jl")
include("src/Heuristic.jl")


# =========================================================================== #

# Function to resolve a SPP instance identified by its filename
function resoudreSPP(fname)
    println("\n... Loading instance ", fname," ...")
    C, A = loadSPP(fname)
    #@show C
    #@show A

    # Choose the Metaheuristics to use------------------------------------------
        # Choose one Construction Metaheuristic 
            GLPK = false
            greedy = false
            GRASP = false
            reactive_GRASP = true
        # Choose one Improvement Metaheuritic
            local_search = false
            tabu_Search = false
    
    # Solving with GLPK----------------------------------
    if GLPK
        cpu_time = @elapsed begin
            println("=============================================================")
            println("... Solving with GLPK ...\n")
            solverSelected = GLPK.Optimizer
            spp = setSPP(C, A)
            set_optimizer(spp, solverSelected)
            optimize!(spp)
        end
        # Displaying the results
        println("GLPK > ", objective_value(spp))
        println("CPU time ", cpu_time)
    end

    # Solving with greedy construction heuristic----------------------------------
    if greedy 
        cpu_time= @elapsed begin
            println("\n=============================================================")
            println("... Solving with greedy construction heuristic ...")
            x0, z = greedy_construction(C, A)    
        end
        # Displaying the results
        println("Greedy > ", z)
        println("CPU time ", cpu_time)
    end

    # Solving with GRASP construction heuristic----------------------------------
    if GRASP 
        cpu_time = @elapsed begin
            println("\n=============================================================")
            println("... Solving with GRASP construction  heuristic ...")
            num_iterations = 25
            alpha = 0.2
            x0, z_best, z_init_values, z_Ls_values, z_best_values = GRASP_construction(C, A, alpha, num_iterations)
                 
        end
        # Displaying the results
        println("GRASP > ", z_best)
        println("CPU time ", cpu_time)
        # Plot the variation of z
        GRASP_plot = plot(z_best_values, label="meilleures solutions", color=:green, lw=3, legend=:bottomright)
        scatter!(1:num_iterations, z_Ls_values, label="toutes solutions améliorées", marker=:utriangle, color=:green)
        scatter!(1:num_iterations, z_init_values, label="toutes solutions construites", marker=:circle, color=:red)
        # Add lines connecting z_init to z_Ls for each iteration
        for iter in 1:num_iterations
            plot!([iter, iter], [z_init_values[iter], z_Ls_values[iter]], color=:black, lw=1, linestyle=:solid, label="")
        end
        
        # Labels and title
        xlabel!("Itérations | zbest = $z_best")
        ylabel!("valeurs de z(x)")
        title!("GRASP-SPP | zInit zLS zBest | $fname")
        savefig("res/GRASP_Plot.png")
        display(GRASP_plot)
    end

    # Solving with reactive GRASP construction heuristic----------------------------------
    if reactive_GRASP
        cpu_time = @elapsed begin
            println("\n=============================================================")
            println("... Solving with reactive_GRASP construction heuristic ...")
            alphas=[0.5, 0.6, 0.7, 0.8]
            num_iterations = 1000
            x0, z_best, alpha_probabilities, z_init_values, z_Ls_values, z_best_values = reactive_GRASP_construction(C, A, alphas, num_iterations)
        end

        # Displaying the results
        println("reactive GRASP > ", z_best)
        println("CPU time ", cpu_time)

        # Plot the alpha probanilities
        alpha_labels = ["α = $(a)" for a in alphas]
        plot2 = pie(alpha_labels, alpha_probabilities, title="Alpha Probabilities")   
        savefig("res/Reactive_GRASP_Proba.png")

        # Plot the variation of z
        reactive_GRASP_plot = plot(z_best_values, label="meilleures solutions", color=:green, lw=3, legend=:bottomright)
        scatter!(1:num_iterations, z_Ls_values, label="toutes solutions améliorées", marker=:utriangle, color=:green)
        scatter!(1:num_iterations, z_init_values, label="toutes solutions construites", marker=:circle, color=:red)
        # Add lines connecting z_init to z_Ls for each iteration
        for iter in 1:num_iterations
            plot!([iter, iter], [z_init_values[iter], z_Ls_values[iter]], color=:black, lw=1, linestyle=:solid, label="")
        end
        
        # Labels and title
        xlabel!("Itérations | zbest = $z_best")
        ylabel!("valeurs de z(x)")
        title!("reactive GRASP-SPP | zInit zLS zBest | $fname")
        savefig("res/Reactive_GRASP_Plot.png")
        display(reactive_GRASP_plot)
    end

    # Improving the result with a local-search heuristic k-p exchange-------------
    if local_search && (greedy || GRASP || reactive_GRASP)
        cpu_time = @elapsed begin
            println("\n=============================================================")
            println("... Improving with Local-Search heuristic ...")
            x0, z = local_search_kp(x0, C, A)
        end
        # Displaying the results
        println("Local Search > ", z)
        println("CPU time ", cpu_time)
    end

    # Improving the result with Tabu-search heuristic ----------------------------
    if tabu_Search && (greedy || GRASP || reactive_GRASP)
        cpu_time = @elapsed begin
            println("\n=============================================================")
            println("... Improving with Tabu-Search heuristic ...")
            x0, z, z_values, iterations = tabu_search(x0, C, A)
        end
        # Displaying the results
        println("Tabu Search > ", z)
        println("CPU time ", cpu_time)
        # Plot the variation of z
        tabo_plot = plot(1:length(z_values), z_values, label="Objective Value (z)", xlabel="Iteration", ylabel="z", title="Variation of z over iterations \n(instance = $fname)")
        savefig("res/Tabu_Search_Plot.png")
        display(tabo_plot)
    end

end

# =========================================================================== #

# Function to perform experimentation on multiple instances
function experimentationSPP()
    println("\n=============================================================")
    println("-> Collecting instance names...")
    target = "dat"
    fnames = getfname(target)
    fnames = first(fnames, 10)

    for fname in fnames
        full_path = "dat/" * fname
        if isfile(full_path)
            println("\nSolving instance:",full_path)
            resoudreSPP(full_path)
        else 
            println("Error: File not found at path ", full_path)
        end
    end
end



