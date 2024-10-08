# title: Create phylogenetic trees to analyze constant population
# author: alexander stein

using SomaticEvolution, Random, AbstractTrees
using Statistics: mean, std, quantile!
using Distributions: Poisson, pdf
using SpecialFunctions: loggamma
using CairoMakie
using Colors

using ColorBrewer
ColPaired = palette("Paired", 12)

# Define some RGB colors for the plots
pastelorange = colorant"rgb(255,150,79)"
banana = colorant"rgb(254,225,53)"
mustard = colorant"rgb(234,170,0)"
darkgreen = colorant"rgb(19,122,99)"
ferngreen = colorant"rgb(79,121,66)"
orangebrown = colorant"rgb(110,71,21)"
darkcyan = colorant"rgb(0,100,100)"
pastelpink = colorant"rgb(222,165,164)"

#########################
### Define parameters ###
#########################
b_1 = 1.0
d_1 = 0.0

mu = 2.0
omega = 2*mu    # effective mutation rate

expNd = 1e4
td = log((b_1)/(b_1-d_1)*expNd + d_1/b_1)/(b_1-d_1)



realizations = 500

###############################
### Analysis of simulations ###
###############################

include("/Users/stein02/Desktop/Phylogenetics/phylo/src/simu_analysis.jl")

################################
### Mathematical predictions ###
################################

include("/Users/stein02/Desktop/Phylogenetics/phylo/src/math_analysis.jl")


################
### Plotting ###
################

function plot_SFS(SFS_xaxis, SFS_simu_single, sSFS_simu, SFS_pred)
    # Simulations
    SFS_mean = mean.(sSFS_simu)
    SFS_sigma = std.(sSFS_simu)

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], xscale=log10, yscale=log10,
        xlabel = "Number of cells, 𝑘", xlabelsize = 25, xticklabelsize = 25,
        ylabel = "Number of mutations, 𝑆ₖ", ylabelsize = 25, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    #xlims!(ax, 0.8, 250)
    ylims!(ax, 0.8, maximum(SFS_mean)*1.2)

    lines!(ax, SFS_xaxis, SFS_mean, color=darkcyan, linewidth=6.0, label = "simulations")
    lowerband = SFS_mean-SFS_sigma
    lowerband = [maximum([x,0.5]) for x in lowerband]
    upperband = SFS_mean+SFS_sigma
    upperband = [maximum([x,0.5]) for x in upperband]
    band!(ax, SFS_xaxis, lowerband, upperband, color = (darkcyan, 0.2))

    scatter!(ax, SFS_xaxis, SFS_simu_single, color = ferngreen, markersize = 15, label = "single realization")

    lines!(ax, SFS_xaxis, SFS_pred, color=pastelorange, linewidth=6.0, linestyle = :dot, label = "prediction")

    #Label(f[1, 1, TopLeft()], "C", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/ExpGrowth/SFS_exp_fixTime.png", f)
end

function plot_SFS_variance(SFS_xaxis, sSFS_simu, SFS_VAR_pred)
    # Simulation variance
    SFS_sigma = std.(sSFS_simu)
    SFS_var = SFS_sigma .* SFS_sigma
    SFS_var = [maximum([x,0.5]) for x in SFS_var]

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], xscale=log10, yscale=log10,
        xlabel = "# of cells, 𝑘", xlabelsize = 25, xticklabelsize = 25,
        ylabel = "var. in # of mutations, Var[𝑆ₖ]", ylabelsize = 25, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"])
    )

    scatter!(ax, SFS_xaxis, SFS_var, color = ferngreen, markersize = 15, label = "simulations")
    lines!(ax, SFS_xaxis, SFS_VAR_pred, color=pastelorange, linewidth=6.0, linestyle = :dot, label = "prediction")
    #lines!(ax, SFS_xaxis, yaxis_pred2, color=:blue, linewidth=6.0, linestyle = :dot, label = L"\text{ prediction}")

    ylims!(ax, 0.2, maximum(SFS_VAR_pred)*2.0)

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/ExpGrowth/SFS_exp_var_fixTime.png", f)
end

function plot_Sk_dist(Sk_simu, k)

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1],
        xlabel = "# of mutations in 20 cells, 𝑆₂₀", xlabelsize = 25, xticklabelsize = 25,
        ylabel = "distribution", ylabelsize = 25, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    hist!(ax, Sk_simu, bins = 15, normalization = :pdf, color = ferngreen, label = "simulations")
    #SFS_k_pred = compound_poisson(SFS_k_xaxis, lam1, 2*mu)
    lines!(ax, Sk_xaxis, Sk_pred2, color=pastelorange, linewidth=6.0, linestyle = :dot, label = "prediction")

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/ExpGrowth/Sk_dist_fixTime.png", f)
end



function plot_tMB(tMB_xaxis, tMB_simu, tMB_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Total mutational burden, 𝘉", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25,
    ytickformat = "{:.6f}"
    #ytickformat = "{:8.1e}"
    #ytickformat = values -> [rich("$(value*1e5)", " x 10",superscript("-5")) for value in values]
    #values -> [L"\sqrt{%$(value^2)}" for value in values])
    )

    hist!(ax, tMB_simu, bins = 20, normalization = :pdf, color=(ColPaired[2],0.8), label = "simulation")
    lines!(ax, tMB_xaxis, tMB_pred, color=:black, linewidth=7.0, linestyle = (:dot,1.5), label = "prediction")

    ylims!(0.0,2e-5)

    #axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figures/ExpGrowth/tMB_fixTime.png", f)
end

function plot_scMBdist(MB_xaxis, sMB_yaxis_series, MB_yaxis_single, MB_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "# mutations, 𝘫", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "# cells with 𝘫 mutations, 𝑀ⱼ", ylabelsize = 30, yticklabelsize = 25
    )

    # Simulations
    MB_mean = mean.(sMB_yaxis_series)
    MB_sigma = std.(sMB_yaxis_series)

    #println("Mean of sMB = ", sum(MB_xaxis.*MB_mean)/sum(MB_mean))
    #println("total number of divisions = ", sum(MB_yaxis_single))

    lines!(ax, MB_xaxis, MB_mean, color=darkcyan, linewidth=6.0, label = "simulations")
    lowerband = MB_mean-MB_sigma
    lowerband = [maximum([x,0.5]) for x in lowerband]
    upperband = MB_mean+MB_sigma
    upperband = [maximum([x,0.5]) for x in upperband]
    band!(ax, MB_xaxis, lowerband, upperband, color = (darkcyan, 0.2))

    scatter!(ax, MB_xaxis, MB_yaxis_single, color = ferngreen, label = nothing, markersize = 15)

    lines!(ax, MB_xaxis, MB_pred, color=pastelorange, linewidth=6.0, linestyle = :dot, label = "prediction")

    axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figures/ExpGrowth/scMBdist_fixTime.png", f)
end


############
### Main ###
############ 

# From the phylogenetic trees we obtain 
#   SFS_collection
#   tMB_collection
#   scMB_collection

### Create variables to store measurements
SFS_collection = []
tMB_collection = Int32[]
scMB_collection = []

Nd = Int64[]

### Run simulations 
for i in 1:realizations
    # Grow tree to Nd
    input = BranchingInput(tmax=td, Nmax=1e12, birthrate=b_1, deathrate=d_1, μ=mu, ploidy=1)
    rng = Random.seed!(i)
    simulation = runsimulation(SimpleTreeCell, input, rng)
    
    # Find root of the tree
    #randomnode = simulation.output.cells[1]
    randomnode = simulation.output.singlemodule[1]
    root = getroot(randomnode)

    # Extract SFS and tMB data
    sitefreqs = get_sitefreqs(root)
    #sitefreqs = get_branchwidths(root)
    push!(SFS_collection, sitefreqs)
    tMB = length(sitefreqs)
    push!(tMB_collection, tMB)
    
    # Extract scMB data
    mb_dist = get_scMBs(root)
    push!(scMB_collection, mb_dist)

    # other information of this population
    push!(Nd,size(simulation.output)[1])
    #println(size(simulation.output)[1])

end

###############################
### Site frequency spectrum ###
###############################

# We have a common xaxis
Nmax = maximum(Nd)
Naverage = mean(Nd)
SFS_xaxis = 1:Nmax

# Compute the prediction
SFS_pred = [predict_SFS_ExpGrowth(x,b_1,d_1,Nmax,omega) for x in SFS_xaxis]

# Create histogram from simulation data
SFS_simu = []
for sitefreqs in SFS_collection
    SFS = create_histo(sitefreqs,SFS_xaxis)
    push!(SFS_simu, SFS)
end
# Group together according to frequencies k
sortedSFS_simu = [Int32[0 for i in 1:realizations] for k in 1:Nmax]
for (i,SFS) in enumerate(SFS_simu)
    for (k,Sk) in enumerate(SFS)
        sortedSFS_simu[k][i] = Sk
    end
end

plot_SFS(SFS_xaxis, SFS_simu[1], sortedSFS_simu, SFS_pred)

# Prediction and plot for the variance
SFS_VAR_pred = SFS_pred*(1+mu)#*(b_1+d_1)/(b_1-d_1)
plot_SFS_variance(SFS_xaxis, sortedSFS_simu, SFS_VAR_pred)

# Prediction and plot for Sk distribution
k = 20
Sk_xaxis, Sk_pred2 = predict_SkDist(k,sortedSFS_simu[k],Nmax,mu)
plot_Sk_dist(sortedSFS_simu[k], k)

###############################
### Total mutational burden ###
###############################

# Create the prediction
tMB_pred = predict_tMB_ExpGrowth(b_1,d_1,expNd,omega)
# Evaulate simulation
tMBmean = mean(tMB_collection)
tMB_std = std(tMB_collection)
# Output
println("predicted tMB = ", tMB_pred)
println("measured tMB = ", tMBmean, " ± ", tMB_std)

# Plot the histogram
# First create xaxis according to simulations
tMB_xaxis = minimum(tMB_collection):maximum(tMB_collection)
# Compute the corresponding yaxis
D = Poisson(omega*(expNd-1))
tMB_pred_histo = [pdf(D,B) for B in tMB_xaxis]
#tMB_pred_histo = predict_tMB_dist(tMB_xaxis, b_1, d_1, tMB_pred/mu,mu)


plot_tMB(tMB_xaxis, tMB_collection, tMB_pred_histo)


##################################################
### Single cell mutational burden distribution ###
##################################################
#=
# We have a common xaxis
highestburden = maximum(maximum.(scMB_collection))
lowestburden = minimum(minimum.(scMB_collection))
scMB_xaxis = lowestburden:highestburden

# Create histogram from simulation data
scMB_simu = []
for scMB in scMB_collection
    scMB_dist = create_histo(scMB, scMB_xaxis)
    push!(scMB_simu, scMB_dist)
end

# Group together according to frequencies k
sortedscMB_simu = [Int32[0 for i in 1:realizations] for k in scMB_xaxis]
# Iterate over realizations
for (r,scMB) in enumerate(scMB_simu)
    # Iterate over scMB_xaxis
    for (j,Mj) in enumerate(scMB)
        sortedscMB_simu[j][r] = Mj
    end
end


# Create the prediction

genAverage = predict_averageGeneration(b_1,d_1,Nd)

scMB_pred = predict_scMBdist(scMB_xaxis, genAverage, mu, Nd)
#scMB_compPoiss = Nd*dist_compound_poisson(scMB_xaxis, genAverage, mu)

# Make the plot
plot_scMBdist(scMB_xaxis, sortedscMB_simu, scMB_simu[1], scMB_pred)
=#

###############################
### Experimental code space ###
###############################




