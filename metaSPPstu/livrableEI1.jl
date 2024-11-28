# =========================================================================== #
include("src/loadSPP.jl")
include("src/getfname.jl")
include("src/Greedy.jl")
include("src/Local_Search.jl")
include("src/VNS.jl")
using LinearAlgebra
# =========================================================================== #

# Function to resolve a SPP instance identified by its filename
function resoudreSPP(fname="pb_100rnd0100.dat")
    println("\n... Loading instance ", fname," ...")
    C, A = loadSPP("dat/"*fname)
    
    # Solving with greedy construction heuristic----------------------------------
    cpu_time= @elapsed begin
        println("\n=============================================================")
        println("... Solving with greedy construction heuristic ...")
        x0, z = greedy_construction(C, A)    
    end
    # Displaying the results
    println("Greedy > ", z)
    println("CPU time ", cpu_time)
    
    # Improving the result with a local-search heuristic k-p exchange-------------
    cpu_time = @elapsed begin
        println("\n=============================================================")
        println("... Improving with Local-Search heuristic ...")
        x0, z = local_search_kp(x0, C, A)
    end
    # Displaying the results
    println("Local Search > ", z)
    println("CPU time ", cpu_time)
    
end

# =========================================================================== #

# Function to perform experimentation on multiple instances
function experimentationSPP()
    println("\n=============================================================")
    println("-> Collecting instance names...")
    target = "dat"
    fnames = ["didactic.txt","pb_100rnd0100.dat","pb_200rnd0100.dat","pb_500rnd1500.dat", "pb_1000rnd0700.dat", "pb_2000rnd0700.dat"]

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



