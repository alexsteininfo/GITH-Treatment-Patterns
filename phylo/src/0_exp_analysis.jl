# title: Create phylogenetic trees to analyze constant population
# author: alexander stein

using SomaticEvolution, Random, AbstractTrees
using Statistics: mean, std, quantile!
using Distributions: Poisson, pdf
using LandauDistribution
using SpecialFunctions: loggamma
using CairoMakie
using Colors

# Define some RGB colors for the plots
pastelorange = colorant"rgb(255,150,79)"
banana = colorant"rgb(254,225,53)"
mustard = colorant"rgb(234,170,0)"
darkgreen = colorant"rgb(19,122,99)"
ferngreen = colorant"rgb(79,121,66)"
orangebrown = colorant"rgb(110,71,21)"
darkcyan = colorant"rgb(0,100,100)"
pastelpink = colorant"rgb(222,165,164)"

using ColorBrewer
ColPaired = palette("Paired", 12)

#########################
### Define parameters ###
#########################
b_1 = 1.0
d_1 = 0.0

mu = 2.0
omega = 2*mu    # effective mutation rate

Nd = Int(1e5)

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
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Number of mutations, 𝑆ₖ", ylabelsize = 30, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    #xlims!(ax, 0.8, 250)
    ylims!(ax, 0.8, maximum(SFS_mean)*1.2)

    scatter!(ax, SFS_xaxis, SFS_simu_single, color = ColPaired[9], markersize = 18, label = "single simulation")

    lines!(ax, SFS_xaxis, SFS_mean, color=ColPaired[2], linewidth=7.0, label = "avg. simulations")
    lowerband = SFS_mean-SFS_sigma
    lowerband = [maximum([x,0.5]) for x in lowerband]
    upperband = SFS_mean+SFS_sigma
    upperband = [maximum([x,0.5]) for x in upperband]
    band!(ax, SFS_xaxis, lowerband, upperband, color = (ColPaired[2], 0.4))

    lines!(ax, SFS_xaxis, SFS_pred, color=:black, linewidth=7.0, linestyle = (:dot, 1.5), label = "prediction")

    #Label(f[1, 1, TopLeft()], "C", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/ExpGrowth/SFS_exp.png", f)
end

function plot_SFS_variance(SFS_xaxis, sSFS_simu, SFS_VAR_pred)
    # Simulation variance
    SFS_sigma = std.(sSFS_simu)
    SFS_var = SFS_sigma .* SFS_sigma
    SFS_var = [maximum([x,0.5]) for x in SFS_var]

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], xscale=log10, yscale=log10,
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Variance, Var[𝑆ₖ]", ylabelsize = 30, yticklabelsize = 25#,
        #yticks = ([10^0,10^2,10^4,10^6])
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e2,1e4,1e6],[L"10^0",L"10^2",L"10^4",L"10^6"])
    )

    scatter!(ax, SFS_xaxis, SFS_var, color = (ColPaired[2], 0.6), markersize = 15, label = "simulations")
    lines!(ax, SFS_xaxis, SFS_VAR_pred, color = :black, linewidth=7.0, linestyle = (:dot, 1.5), label = "prediction")
    #lines!(ax, SFS_xaxis, yaxis_pred2, color=:blue, linewidth=6.0, linestyle = :dot, label = L"\text{ prediction}")

    ylims!(ax, 0.2, maximum(SFS_VAR_pred)*2.0)

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/ExpGrowth/SFS_exp_var.png", f)
end

function plot_Sk_dist(Sk_simu, k)

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1],
        xlabel = "Number of mut. in 20 cells, 𝑆₂₀", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Density", ylabelsize = 30, yticklabelsize = 25,
        xticks = ([850,950,1050])
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    hist!(ax, Sk_simu, bins = 15, normalization = :pdf, color = (ColPaired[2],0.8), label = "simulations")
    #SFS_k_pred = compound_poisson(SFS_k_xaxis, lam1, 2*mu)
    lines!(ax, Sk_xaxis, Sk_pred2, color= :black, linewidth=7.0, linestyle = (:dot, 1.5), label = "prediction")

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/ExpGrowth/Sk_dist.png", f)
end



function plot_tMB(tMB_xaxis, tMB_simu, tMB_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Total mutational burden, 𝘉", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25,
    xticks = ( [683000, 688000, 693000] )
    #yticks = ([0.0000,0.0002,0.0004,0.0006])
    )

    #xlims!(397800, 402500)

    hist!(ax, tMB_simu, bins = 20, normalization = :pdf, color=(ColPaired[2],0.8), label = "simulations")
    lines!(ax, tMB_xaxis, tMB_pred, color=:black, linewidth=7.0, linestyle = (:dot,1.5), label = "comp. Pois.")

    #axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figures/ExpGrowth/tMB.png", f)
end

function plot_scMBdist(MB_xaxis, sMB_yaxis_series, MB_yaxis_single, MB_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Number of mutations, 𝘫", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Number of cells, 𝑀ⱼ", ylabelsize = 30, yticklabelsize = 25
    )

    # Simulations
    MB_mean = mean.(sMB_yaxis_series)
    MB_sigma = std.(sMB_yaxis_series)

    #println("Mean of sMB = ", sum(MB_xaxis.*MB_mean)/sum(MB_mean))
    #println("total number of divisions = ", sum(MB_yaxis_single))

    scatter!(ax, MB_xaxis, MB_yaxis_single, color = ColPaired[9], markersize = 15, label = "single simulation")

    lines!(ax, MB_xaxis, MB_mean, color=ColPaired[2], linewidth=6.0, label = "avg. simulations")
    lowerband = MB_mean-MB_sigma
    lowerband = [maximum([x,0.5]) for x in lowerband]
    upperband = MB_mean+MB_sigma
    upperband = [maximum([x,0.5]) for x in upperband]
    band!(ax, MB_xaxis, lowerband, upperband, color = (ColPaired[2], 0.2))

    lines!(ax, MB_xaxis, MB_pred, color=:black, linewidth=7.0, linestyle = (:dot,1.5), label = "prediction")

    axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figures/ExpGrowth/scMBdist.png", f)
end


function plot_subpopsizes(kappas_simu, kesslers_simu, xaxis, kappas_pred, kesslers_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xticks = ([0,1000,2000,3000,4000]),
    xlabel = "Number of cells", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25
    )

    hist!(kappas_simu, color=(:green,0.4), normalization=:pdf, bins=600, label="P(κ)")
    hist!(kesslers_simu, color=(:blue,0.4), normalization=:pdf, bins=600, label="P(𝘙)")

    lines!(xaxis, kappas_pred, color=:green)
    lines!(xaxis, kesslers_pred, color=:blue)

    xlims!(0,5500)

    axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figures/ExpGrowth/subpopsizes.png", f)
end

function plot_exponents(exponents)
    f = Figure()
    ax = Axis(f[1,1],
    #xticks = ([0,1000,2000,3000,4000]),
    xlabel = "Exponent, α", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25
    )

    hist!(exponents, color = (ColPaired[2], 0.6), normalization=:pdf, bins=20)
    #axislegend(ax, labelsize=25, framecolor = :white)

    save("figures/ExpGrowth/exponenthist.png", f)
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

### Run simulations 
for i in 1:realizations
    # Grow tree to Nd
    input = BranchingInput(Nmax=Nd, birthrate=b_1, deathrate=d_1, μ=mu, ploidy=1)
    #rng = Random.seed!(i)
    rng = Random.seed!()
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

end

###############################
### Site frequency spectrum ###
###############################

# We have a common xaxis
SFS_xaxis = 1:Nd

# Compute the prediction
SFS_pred = [predict_SFS_ExpGrowth(x,b_1,d_1,Nd,omega) for x in SFS_xaxis]

# Create histogram from simulation data
SFS_simu = []
for sitefreqs in SFS_collection
    SFS = create_histo(sitefreqs,SFS_xaxis)
    push!(SFS_simu, SFS)
end
# Group together according to frequencies k
sortedSFS_simu = [Int32[0 for i in 1:realizations] for k in 1:Nd]
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
Sk_xaxis, Sk_pred2 = predict_SkDist(k,sortedSFS_simu[k],Nd,mu)
plot_Sk_dist(sortedSFS_simu[k], k)

###############################
### Total mutational burden ###
###############################

# Create the prediction
tMB_pred = predict_tMB_ExpGrowth(b_1,d_1,Nd,omega)
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
#tMB_pred_histo = predict_tMB_dist(tMB_xaxis,b_1,d_1,tMB_pred/mu,mu)
D = Poisson(tMB_pred)
tMB_pred_histo = [pdf(D,B) for B in tMB_xaxis]

plot_tMB(tMB_xaxis, tMB_collection, tMB_pred_histo)


##################################################
### Single cell mutational burden distribution ###
##################################################

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


#########################
### Fitting exponents ###
#########################

kmin = 100
exponents = [get_exponent(SFS, kmin) for SFS in SFS_simu]

println("Mean exponent = ", mean(exponents))
println("Std exponent = ", std(exponents))

plot_exponents(exponents)


###############################
### Experimental code space ###
###############################

kappas_sim = maximum.(SFS_collection)
kesslers_sim = sum.(SFS_collection)

xaxis_popsizes = 1:Nd
sigma_inv = b_1/(b_1-d_1)
sigma = (b_1-d_1)/b_1
numres = 2*mu*Nd*sigma_inv
numres0 = 2*mu*Nd
kappas_pred = [numres/kappa^2*(1-1/kappa)^(numres-1) for kappa in xaxis_popsizes]
Lan = Landau(1.0)
kesslers_pred = [1/numres*pdf(Lan, m/numres - log(numres0) ) for m in xaxis_popsizes]

plot_subpopsizes(kappas_sim, kesslers_sim, xaxis_popsizes, kappas_pred, kesslers_pred)
