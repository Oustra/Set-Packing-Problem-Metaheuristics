# =========================================================================== #
using LinearAlgebra
using PyPlot
include("src/SA.jl")
include("src/loadSPP.jl")
include("src/ACO.jl")
# =========================================================================== #

function ACO_SPP(fname, nbIteration, maxAnt, rhoE, rhoD, phiInit)
    C, A = loadSPP(fname)
    return ant_colony_optimization(C, A, nbIteration, maxAnt, rhoE, rhoD, phiInit)
end

function plotRunACO(iname, zbest, ziteration, all_solutions, maxAnt)
    zbetter = maximum(zbest)
    figure("Examen d'un run", figsize=(6, 6))
    title("ACO-SPP | " * iname)
    xlabel("#Iterations x #Ants | zbest = $zbetter")
    ylabel("Valeurs de z(x)")
    ylim(0, zbetter + 10)

    nPoint = length(all_solutions)
    x = collect(1:nPoint)
    xticks([1,convert(Int64,ceil(nPoint/4)),convert(Int64,ceil(nPoint/2)), convert(Int64,ceil(nPoint/4*3)),nPoint])
    scatter(x, all_solutions, color="black", s=0.7, label="AllSolutions")
    plot(x, zbest, linewidth=3.0, color="red", label="BestGlobal")
    plot(x, ziteration, ls="--", marker="o", ms=2, color="black", lw=0.5, label="BestIteration")
    

    legend(loc="lower right", fontsize="small")
    savefig("res/ACO/ACO_Plot.jpg")
    println("Plot saved as res/ACO_Plot.jpg")
end

function plotCPUt(allfinstance, tmoy)
    figure("bilan CPUt tous runs",figsize=(6,6))
    title("ACO-SPP | tMoy")
    ylabel("CPUt moyen (s)")

    xticks(collect(1:length(allfinstance)), allfinstance, rotation=60, ha="right")
    margins(0.15)
    subplots_adjust(bottom=0.15,left=0.21)
    plot(collect(1:length(allfinstance)),tmoy,linestyle="--", lw=0.5, marker="o", ms=4, color="blue", label="tMoy")
    legend(loc=4, fontsize ="small")
    savefig("res/ACO/ACO_CPUt_Plot.jpg")
end

function plotZvalues(allfinstance, allfinstanceZ)
    figure("bilan zbest tous runs",figsize=(6,6))
    title("ACO-SPP | zbest")
    ylabel("z values")

    xticks(collect(1:length(allfinstance)), allfinstance, rotation=60, ha="right")
    margins(0.15)
    subplots_adjust(bottom=0.15,left=0.21)
    plot(collect(1:length(allfinstance)),allfinstanceZ,linestyle="--", lw=0.5, marker="o", ms=4, color="red", label="zbest")
    legend(loc=4, fontsize ="small")
    savefig("res/ACO/ACO_zbest_Plot.jpg")
end

function plotAnalyseACO(iname, x, zmoy, zmin, zmax)
    zbetter = maximum(zmax)
    figure("bilan tous runs",figsize=(6,6))
    title("ACO-SPP | \$z_{min}\$  \$z_{moy}\$  \$z_{max}\$ | " * iname)
    xlabel("#Iteration x #Ants | zbest = $zbetter")
    ylabel("valeurs de z(x)")
    ylim(0, zmax[end]+2)

    nPoint = length(x)
    intervalle = [reshape(zmoy,(1,nPoint)) - reshape(zmin,(1,nPoint)) ; reshape(zmax,(1, nPoint))-reshape(zmoy,(1,nPoint))]
    xticks(x)
    errorbar(x,zmoy,intervalle,lw=1, color="black", label="zMin zMax")
    plot(x,zmoy,linestyle="-", marker="o", ms=4, color="green", label="zMoy")
    legend(loc=4, fontsize ="small")
    savefig("res/ACO/ACO_Analyse_Plot.jpg")
end

function plotEvolutionMetrics(iname, x, zmoy, zmin, zmax, nbRunGrasp)
    zbetter = maximum(zmax)
    figure("ACO-SPP Analysis", figsize=(6, 6))
    title("ACO-SPP | zₘᵢₙ  zₘₒᵧ  zₘₐₓ | " * iname)
    xlabel("#Itérations x #Ants | nRun = $nbRunGrasp | zbest = $zbetter")
    ylabel("Valeurs de z(x)")

    plot(x, zmin, "r:", marker="o", ms=4, label="zMin")
    plot(x, zmoy, "b:", marker="o", ms=4, label="zMoy")
    plot(x, zmax, "g:", marker="o", ms=4, label="zMax")
    legend(loc="lower right", fontsize="small")
    ylim(minimum(zmin) - 2, maximum(zmax) + 2)
    savefig("res/ACO/ACO_Evolution.jpg")
end

# =========================================================================== #
# Simulation of an experimentation

function simulation()
    #allfinstance     =  ["dat/didactic.txt", "dat/pb_100rnd0100.dat", "dat/pb_100rnd0500.dat", "dat/pb_200rnd0100.dat",  "dat/pb_200rnd0700.dat", "dat/pb_200rnd1500.dat", "dat/pb_500rnd0600.dat", "dat/pb_500rnd0800.dat", "dat/pb_500rnd1500.dat", "dat/pb_1000rnd0700.dat", "dat/pb_2000rnd0700.dat"]
    allfinstance      =  ["dat/pb_100rnd0100.dat"]
    nbInstances       =  length(allfinstance)
    nbRunACO          =  3     # Number of runs ACO
    nbIteration       =  100   # Number of iteration ACO
    nbDivisionRun     =  10    # Number of division run ACO
    maxAnt            =  10    # Number of Ants
    rhoE              =  0.8
    rhoD              =  0.2
    phiInit           =  1.0

    all_solutions = zeros(Int64, nbIteration * maxAnt)
    zbest = zeros(Int64, nbIteration * maxAnt)  
    ziteration = zeros(Int64, nbIteration * maxAnt)              
    zbetter = 0

    x     = zeros(Int64, nbDivisionRun)
    zmax  = Matrix{Int64}(undef,nbInstances , nbDivisionRun); zmax[:] .= typemin(Int64)  # -Inf entier
    zmoy  = zeros(Float64, nbInstances, nbDivisionRun) # zero
    zmin  = Matrix{Int64}(undef,nbInstances , nbDivisionRun) ; zmin[:] .= typemax(Int64)  # +Inf entier
    tmoy  = zeros(Float64, nbInstances) 
    allfinstanceZ = zeros(Int64, nbInstances) #zero
    cpt = 0
    
    for division=1:nbDivisionRun
        x[division] = convert(Int64, ceil(nbIteration * maxAnt / nbDivisionRun * division))
    end

    for instance = 1:nbInstances
        print("-> ",allfinstance[instance]," : \n")
        for runACO in 1:nbRunACO
            @show runACO
            start = time()
            zbest, ziteration, all_solutions = ACO_SPP(allfinstance[instance], nbIteration, maxAnt, rhoE, rhoD, phiInit)
            tutilise = time()-start
            cpt+=1;
            
            for division=1:nbDivisionRun
                zmax[instance,division] = max(zbest[x[division]], zmax[instance,division])
                zmin[instance,division] = min(zbest[x[division]], zmin[instance,division])
                zmoy[instance,division] =  zbest[x[division]] + zmoy[instance,division]
            end

            tmoy[instance] = tmoy[instance] + tutilise
        end
        for division=1:nbDivisionRun
            zmoy[instance,division] =  zmoy[instance,division] /  nbRunACO
        end
        tmoy[instance] = tmoy[instance] / nbRunACO
        allfinstanceZ[instance] = maximum(zbest)
    end

    instancenb = nbInstances
    plotRunACO(allfinstance[instancenb], zbest, ziteration, all_solutions, maxAnt)
    plotAnalyseACO(allfinstance[instancenb], x, zmoy[instancenb,:], zmin[instancenb,:], zmax[instancenb,:] )
    plotEvolutionMetrics(allfinstance[instancenb], x, zmoy[instancenb,:], zmin[instancenb,:], zmax[instancenb,:], nbRunACO)
    plotZvalues(allfinstance, allfinstanceZ)
    plotCPUt(allfinstance, tmoy)
    
    println("Experimentation ACO-SPP avec :")
    println("  nbInstances       = ", nbInstances)
    println("  nbRunACO          = ", nbRunACO)
    println("  nbIterationACO    = ", nbIteration)
    println("  nbDivisionRun     = ", nbDivisionRun)
    println("  nbAnts            = ", maxAnt)
    println("  Graphiques de synthese...")

end