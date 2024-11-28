# metaSPPstu: Solving the Set Packing Problem

The **Set Packing Problem (SPP)** is a classic combinatorial optimization problem where the goal is to select the maximum number of disjoint subsets from a given collection, ensuring that no two subsets overlap. It has wide-ranging applications, including resource allocation, scheduling, and network design.

This repository provides an implementation of tools and metaheuristics to solve the SPP using Julia.

---

## Features

- **SPP Data Handling**:
  - `loadSPP.jl`: Reads SPP instances in OR-library format.
  - `setSPP.jl`: Constructs a JuMP model for SPP.
  - `getfname.jl`: Collects names of non-hidden files in a given directory.
  - `experiment.jl`: Implements protocols for numerical experiments with graphical outputs.

---

## Data

The `dat/` directory contains various numerical instances of the SPP in OR-library format, such as:

- `didactic.dat`
- `pb_100rnd0100.dat`
- `...`
- `pb_2000rnd0800.dat`

---

## Implemented Metaheuristics

The `src/` directory includes the following metaheuristic algorithms for solving SPP:

- **Greedy Construction**
- **Local Search k-p Exchange**:
  - Variants: (2-1), (1-1), (0-1)
- **GRASP (Greedy Randomized Adaptive Search Procedure)**
- **Reactive GRASP**
- **Ant Colony Optimization (ACO)**
- **Tabu Search**
- **Variable Neighborhood Search (VNS)**
- **Simulated Annealing (SA)**

---

## How to Use

### `livrableEI1.jl`

This script provides example usage for solving SPP instances using Greedy construction and Local search for improvement. 

#### Load the file:
    include("livrableEI1.jl")

- Solve the problem:
    resoudreSPP("pb_100rnd0100.dat")
    resoudreSPP() # Default instance: "pb_100rnd0100.dat"

- Run an experiment:
    experimentationSPP()
    # On instances ["didactic.txt", "pb_100rnd0100.dat", "pb_200rnd0100.dat", "pb_500rnd1500.dat", "pb_1000rnd0700.dat", "pb_2000rnd0700.dat"]
    
### GRASP Metaheuristic

Example usage (experiment.jl) for running an experimental protocol on a simulated GRASP-SPP:

#### Load the file:
    include("experiment.jl")

#### Run a simulation:
    simulation() # Default instance: "pb_100rnd0100.dat"
    # On the instance "pb_100rnd0100.dat"

- Plot results are stored in the directory res/GRASP/.


### Reactive GRASP Metaheuristic

Example usage (experiment2.jl) for running an experimental protocol on a simulated reactive GRASP-SPP:

####  Load the file:
   - include("experiment2.jl")

#### Run a simulation:
   - simulation() # Default instance: "pb_100rnd0100.dat"
   - On the instance "pb_100rnd0100.dat"

#### Plot results are stored in the directory res/ReactiveGRASP/.


### Ant Colony Metaheuristic

Example usage (experiment3.jl) for running an experimental protocol on a simulated ACO-SPP:

#### Load the file:
  - include("experiment3.jl")

#### Run a simulation:
  - simulation() # Default instance: "pb_100rnd0100.dat"
  - On the instance "pb_100rnd0100.dat"

#### Plot results are stored in the directory res/ACO/.

---------------------------------------------------------------------------------------------------
### All Results

### The results for Zmin, Zavg, Zmax, and Tavg for all instances are stored in the file All-Results.txt.
