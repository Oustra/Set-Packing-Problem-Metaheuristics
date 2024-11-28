# =========================================================================== #
using LinearAlgebra
using PyPlot
include("src/GRASP.jl")
include("src/Local_Search.jl")
include("src/loadSPP.jl")
# =========================================================================== #

function graspSPP(fname, alpha, nbIterationGrasp)
    C, A = loadSPP(fname)
    return GRASP_construction(C, A, alpha, nbIterationGrasp)
end

function plotRunGrasp(iname,zinit, zls, zbest)
    zbetter = maximum(zbest)
    figure("Examen d'un run",figsize=(6,6))
    title("GRASP-SPP | \$z_{Init}\$  \$z_{LS}\$  \$z_{Best}\$ | " * iname)
    xlabel("Itérations | zbest = $zbetter")
    ylabel("valeurs de z(x)")
    ylim(0, zbetter+2)

    nPoint = length(zinit)
    x=collect(1:nPoint)
    xticks([1,convert(Int64,ceil(nPoint/4)),convert(Int64,ceil(nPoint/2)), convert(Int64,ceil(nPoint/4*3)),nPoint])
    plot(x,zbest, linewidth=2.0, color="green", label="meilleures solutions")
    plot(x,zls,ls="",marker="^",ms=2,color="green",label="toutes solutions améliorées")
    plot(x,zinit,ls="",marker=".",ms=2,color="red",label="toutes solutions construites")
    vlines(x, zinit, zls, linewidth=0.5)
    legend(loc=4, fontsize ="small")
    savefig("res/GRASP/GRASP_Plot.jpg")
end

function plotAnalyseGrasp(iname, x, zmoy, zmin, zmax)
    zbetter = maximum(zmax)
    figure("bilan tous runs",figsize=(6,6))
    title("GRASP-SPP | \$z_{min}\$  \$z_{moy}\$  \$z_{max}\$ | " * iname)
    xlabel("Itérations (pour nbRunGrasp) | zbest = $zbetter")
    ylabel("valeurs de z(x)")
    ylim(0, zmax[end]+2)

    nPoint = length(x)
    intervalle = [reshape(zmoy,(1,nPoint)) - reshape(zmin,(1,nPoint)) ; reshape(zmax,(1, nPoint))-reshape(zmoy,(1,nPoint))]
    xticks(x)
    errorbar(x,zmoy,intervalle,lw=1, color="black", label="zMin zMax")
    plot(x,zmoy,linestyle="-", marker="o", ms=4, color="green", label="zMoy")
    legend(loc=4, fontsize ="small")
    savefig("res/GRASP/GRASP_Analyse_Plot.jpg")
end

function plotEvolutionMetrics(iname, x, zmoy, zmin, zmax, nbRunGrasp)
    zbetter = maximum(zmax)
    figure("GRASP-SPP Analysis", figsize=(6, 6))
    title("GRASP-SPP | zₘᵢₙ  zₘₒᵧ  zₘₐₓ | " * iname)
    xlabel("Itérations | nRun = $nbRunGrasp | zbest = $zbetter")
    ylabel("Valeurs de z(x)")

    plot(x, zmin, "r:", marker="o", ms=4, label="zMin")
    plot(x, zmoy, "b:", marker="o", ms=4, label="zMoy")
    plot(x, zmax, "g:", marker="o", ms=4, label="zMax")
    legend(loc="lower right", fontsize="small")
    ylim(minimum(zmin) - 2, maximum(zmax) + 2)
    savefig("res/GRASP/GRASP_Evolution.jpg")
end

function plotCPUt(allfinstance, tmoy)
    figure("bilan CPUt tous runs",figsize=(6,6))
    title("GRASP-SPP | tMoy")
    ylabel("CPUt moyen (s)")

    xticks(collect(1:length(allfinstance)), allfinstance, rotation=60, ha="right")
    margins(0.15)
    subplots_adjust(bottom=0.15,left=0.21)
    plot(collect(1:length(allfinstance)),tmoy,linestyle="--", lw=0.5, marker="o", ms=4, color="blue", label="tMoy")
    legend(loc=4, fontsize ="small")
    savefig("res/GRASP/GRASP_CPUt_Plot.jpg")
end

function plotZvalues(allfinstance, allfinstanceZ)
    figure("bilan zbest tous runs",figsize=(6,6))
    title("GRASP-SPP | zbest")
    ylabel("z values")

    xticks(collect(1:length(allfinstance)), allfinstance, rotation=60, ha="right")
    margins(0.15)
    subplots_adjust(bottom=0.15,left=0.21)
    plot(collect(1:length(allfinstance)),allfinstanceZ,linestyle="--", lw=0.5, marker="o", ms=4, color="red", label="zbest")
    legend(loc=4, fontsize ="small")
    savefig("res/GRASP/GRASP_zbest_Plot.jpg")
end

# =========================================================================== #
# Simulation of an experimentation

function simulation()
    #allfinstance     =  ["dat/didactic.txt", "dat/pb_100rnd0100.dat", "dat/pb_100rnd0500.dat", "dat/pb_200rnd0100.dat",  "dat/pb_200rnd0700.dat", "dat/pb_200rnd1500.dat", "dat/pb_500rnd0600.dat", "dat/pb_500rnd0800.dat", "dat/pb_500rnd1500.dat", "dat/pb_1000rnd0700.dat", "dat/pb_2000rnd0700.dat"]
    allfinstance      =  ["dat/pb_100rnd0100.dat"]
    nbInstances       =  length(allfinstance)
    nbRunGrasp        =  3    # Number of runs GRASP
    nbIterationGrasp  =  200  # Number of iterations GRASP
    nbDivisionRun     =  10   # Number of division run GRASP
    alpha             =  0.5

    zinit = zeros(Int64, nbIterationGrasp) # zero
    zls   = zeros(Int64, nbIterationGrasp) # zero
    zbest = zeros(Int64, nbIterationGrasp) # zero
 
    x     = zeros(Int64, nbDivisionRun)
    zmax  = Matrix{Int64}(undef,nbInstances , nbDivisionRun); zmax[:] .= typemin(Int64)  # -Inf entier
    zmoy  = zeros(Float64, nbInstances, nbDivisionRun) # zero
    zmin  = Matrix{Int64}(undef,nbInstances , nbDivisionRun) ; zmin[:] .= typemax(Int64)  # +Inf entier
    tmoy  = zeros(Float64, nbInstances)  # zero
    allfinstanceZ = zeros(Float64, nbInstances) #zero

    for division=1:nbDivisionRun
        x[division] = convert(Int64, ceil(nbIterationGrasp / nbDivisionRun * division))
    end

    println("Experimentation GRASP-SPP avec :")
    println("  nbInstances       = ", nbInstances)
    println("  nbRunGrasp        = ", nbRunGrasp)
    println("  nbIterationGrasp  = ", nbIterationGrasp)
    println("  nbDivisionRun     = ", nbDivisionRun)
    println("  alphaValue        = ", alpha)
    println(" ")
    cpt = 0

    for instance = 1:nbInstances
        
        print("-> ",allfinstance[instance]," : \n")
        for runGrasp = 1:nbRunGrasp

            start = time()
            zinit, zls, zbest = graspSPP(allfinstance[instance], alpha, nbIterationGrasp)
            tutilise = time()-start
            cpt+=1
            println("GRASP Run: ",cpt%10)

            for division=1:nbDivisionRun
                zmax[instance,division] = max(zbest[x[division]], zmax[instance,division])
                zmin[instance,division] = min(zbest[x[division]], zmin[instance,division])
                zmoy[instance,division] =  zbest[x[division]] + zmoy[instance,division]
            end 
            tmoy[instance] = tmoy[instance] + tutilise

        end
        for division=1:nbDivisionRun
             zmoy[instance,division] =  zmoy[instance,division] /  nbRunGrasp
        end
        tmoy[instance] = tmoy[instance] / nbRunGrasp
        allfinstanceZ[instance] = maximum(zbest)
        println(" ")

    end

    instancenb = nbInstances
    plotRunGrasp(allfinstance[instancenb], zinit, zls, zbest)
    plotAnalyseGrasp(allfinstance[instancenb], x, zmoy[instancenb,:], zmin[instancenb,:], zmax[instancenb,:] )
    plotEvolutionMetrics(allfinstance[instancenb], x, zmoy[instancenb,:], zmin[instancenb,:], zmax[instancenb,:], nbRunGrasp)
    plotZvalues(allfinstance, allfinstanceZ)
    plotCPUt(allfinstance, tmoy)

    println("Experimentation GRASP-SPP avec :")
    println("  nbInstances       = ", nbInstances)
    println("  nbRunGrasp        = ", nbRunGrasp)
    println("  nbIterationGrasp  = ", nbIterationGrasp)
    println("  nbDivisionRun     = ", nbDivisionRun)
    println("  alphaValue        = ", alpha)
    println("  Graphiques de synthese...")
end
