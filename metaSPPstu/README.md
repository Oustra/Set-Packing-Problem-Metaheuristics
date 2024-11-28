---------------------------------------------------------------------------------------------------
# metaSPPstu: Solving the Set Packing Problem
Implementation in Julia (compatible with Julia v1.x) of tools related to the Set Packing Problem (SPP) for educational purposes.

This implementation serves as the foundation for an exercise in the "Metaheuristics" course.

---------------------------------------------------------------------------------------------------
# Features of the Code

Support elements for implementing homework assignments in Julia v1.x as part of the "Metaheuristics" course for the first year of a Master’s in Computer Science. Updated for the 2024-2025 academic year.

- loadSPP.jl: Reads an SPP instance in OR-library format.
- setSPP.jl: Constructs a JuMP model for SPP.
- getfname.jl: Collects the names of non-hidden files in a given directory.
- experiment.jl: Protocol for conducting numerical experiments with graphical outputs.

---------------------------------------------------------------------------------------------------
# Data Used

The Data directory contains a selection of numerical SPP instances in OR-library format:

- didactic.dat
- pb_100rnd0100.dat
- ...
- pb_2000rnd0800.dat

---------------------------------------------------------------------------------------------------
# Available Metaheuristics 

The Source directory contains a selection of Metaheuritics :

- Greedy Construction 
- Local Search k-p Exchange (2-1) (1-1) (0-1)
- Reactive GRASP
- GRASP 
- Ant Colony Optimization
- Tabu Search
- VNS
- SA

---------------------------------------------------------------------------------------------------
# livrableEI1

Example usage (`'livrableEI1.jl`) with file paths matching a standard configuration for macOS/Windows:

- Load the file:
    include("livrableEI1.jl")

- Solve the problem:
    resoudreSPP("pb_100rnd0100.dat")
    resoudreSPP() # Default instance: "pb_100rnd0100.dat"

- Run an experiment:
    experimentationSPP()
    # On instances ["didactic.txt", "pb_100rnd0100.dat", "pb_200rnd0100.dat", "pb_500rnd1500.dat", "pb_1000rnd0700.dat", "pb_2000rnd0700.dat"]
    
---------------------------------------------------------------------------------------------------
# GRASP Metaheuristic

Example usage (experiment.jl) for running an experimental protocol on a simulated GRASP-SPP:

- Load the file:
    include("experiment.jl")

- Run a simulation:
    simulation() # Default instance: "pb_100rnd0100.dat"
    # On the instance "pb_100rnd0100.dat"

- Plot results are stored in the directory res/GRASP/.

---------------------------------------------------------------------------------------------------
# Reactive GRASP Metaheuristic

Example usage (experiment2.jl) for running an experimental protocol on a simulated reactive GRASP-SPP:

- Load the file:
    include("experiment2.jl")

- Run a simulation:
    simulation() # Default instance: "pb_100rnd0100.dat"
    # On the instance "pb_100rnd0100.dat"

- Plot results are stored in the directory res/ReactiveGRASP/.

---------------------------------------------------------------------------------------------------
# Ant Colony Metaheuristic

Example usage (experiment3.jl) for running an experimental protocol on a simulated ACO-SPP:

- Load the file:
    include("experiment3.jl")

- Run a simulation:
    simulation() # Default instance: "pb_100rnd0100.dat"
    # On the instance "pb_100rnd0100.dat"

- Plot results are stored in the directory res/ACO/.

---------------------------------------------------------------------------------------------------
# All Results

- The results for Zmin, Zavg, Zmax, and Tavg for all instances are stored in the file All-Results.txt.