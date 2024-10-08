# title: Create phylogenetic trees to analyze resistance dynamics
# author: alexander stein

#include("/Users/stein02/Desktop/github/TreeStatistics.jl/src/SomaticEvolution.jl")
#using .SomaticEvolution
# using Revise
# Comment: Use SomaticEvolution.function to use functionalities of SomaitcEvolution to avoid conflicts

using SomaticEvolution, Random, AbstractTrees
using Statistics: mean, std, quantile!
using Distributions: Poisson, Binomial, DiscreteNonParametric, pdf
using SpecialFunctions: loggamma    # for the Consul Poisson distribution
using CairoMakie
using Colors
using Integrals
using LandauDistribution


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

# Parameters before treatment
b_s1 = 1.0
d_s1 = 0.0
b_r1 = 1.0
d_r1 = 0.0

# rate of neutral mutations
mu = 2.0
omega = 2*mu
# rate of resistance mutations
nu = 2e-3
# resistance mutation distribution
dist1 = DiscreteNonParametric([0.0],[1.0])

# Maximum number of subclones
maxclones = 1000

# other parameters
Nd = Int(1e5)

timevector = Float64[]

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

function plot_subpopsizes(kappas_simu, kesslers_simu, xaxis, kappas_pred, kesslers_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xticks = ([0,1000,2000,3000]),
    xlabel = "Number of cells", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25
    )

    hist!(kappas_simu, color=(:green,0.4), normalization=:pdf, bins=600, label="P(κ)")
    hist!(kesslers_simu, color=(:blue,0.4), normalization=:pdf, bins=600, label="P(𝘙)")

    lines!(xaxis, kappas_pred, color=:green)
    lines!(xaxis, kesslers_pred, color=:blue)

    xlims!(0,3500)

    axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figures/compare/subpopsizes222.png", f)
end

function plot_sizedifference_simulated(samples)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Difference in size", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25#, #yscale=log10#,
    #xticks = ([-2e5, 0.0, 2e5, 4e5], [rich("-2×10", superscript("5")), "0", rich("2×10", superscript("5")), rich("4×10", superscript("5"))]),
    #yticks = ([0, 5e-6, 1e-5, 1.5e-5, 2.0e-5], ["0", rich("5.0×10", superscript("-6")), rich("1.0×10", superscript("-5")), rich("1.5×10", superscript("-5")),rich("2.0×10", superscript("-5"))])
    )

    bins = [i*0.02*1e4 - 1e4 for i in 1:100]
    reds = [(ColPaired[6],0.7) for i in 1:49]
    blues = [(ColPaired[2],0.7) for i in 1:50]
    col_redblue = [reds;blues]
    hist!(ax, samples, bins=bins, normalization = :pdf, color = col_redblue)

    xlims!(-3e3,5e3)


    save("figures/compare/SizeDifference_simulated.png", f)    
end


############
### Main ###
############

# From the phylogenetic trees we obtain 
#   SFS_collection_pre, SFS_collection_new
#   tMB_collection_pre, tMB_collection_new
#   sMB_collection_pre, sMB_collection_new

SFS_collection = []
tMB_collection = Int32[]
scMB_collection = []

# Analysis of clonesizes
kappa_tilde = Int32[]   # data taken from clonesizes
t_tilde = Float64[]

CSD_collection = Any[]
ancestry_collection = Any[]

# Grow population to mutant arrival size
#input = BranchingInput(Nmax=Nmut, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
#simulation = runsimulation(SimpleTreeCell, input)

for i in 1:realizations
    # Set a random variable 
    rng = Random.seed!(i)

    # Grow population to detection size
    input = BranchingInput(Nmax=Nd, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
    selection = SelectionDistribution(dist1, nu, maxclones)
    global simulation = runsimulation(SimpleTreeCell, input, selection, rng)

    #println(treebreadth(root))
    #println(getsubclonesizes(simulation))
    
    # Extract and save distributions
    #randomnode = simulation.output.cells[1]
    randomnode = simulation.output.singlemodule[1]
    root = getroot(randomnode)
    global t_detec = age(simulation.output)
    global SenPop = getsubclonesizes(simulation)[1]

    i == 1 ? push!(timevector, t_detec) : nothing
    #println("Subclone sizes:")
    modules = simulation.output.singlemodule
    #println(treebreadth(root))
    #println(getsubclonesizes(simulation))

    # SFS
    #sitefreqs_pre, sitefreqs_new = get_sitefreqs(root, t_detec)
    sitefreqs_pre, sitefreqs_new = get_sitefreqs_selection(root)
    # after phase 1, all mutations occured before t_detec
    # thus, sitefreq_pre will be full and sitefreq_new will be empty
    push!(SFS_collection, sitefreqs_pre)

    #sitefreqs = getallelefreq(simulation)

    # tMB
    tMB = length(sitefreqs_pre)
    push!(tMB_collection, tMB)
    
    # scMB
    mb_dist = get_scMBs(root)
    push!(scMB_collection, mb_dist)

    ### clone size information
    
    clonesizes = getsubclonesizes(simulation)
    push!(CSD_collection, clonesizes)
    ancestry = Int32[]
    for subclone in simulation.output.subclones
        push!(ancestry, subclone.parentid)
    end
    push!(ancestry_collection, ancestry)
    biggestclone_size = maximum(clonesizes)
    #push!(kappa_tilde, biggestclone_size)
    index = findfirst(item -> item == biggestclone_size, clonesizes)
    biggestclone_time = simulation.output.subclones[index].mutationtime
    push!(t_tilde, biggestclone_time)
    

end

### Consider only resistant clones, not the sensitive one.
CSD_collection2 = [collection[2:end] for collection in CSD_collection]

### Kappa versus Kessler
kappas_sim = maximum.(CSD_collection2)
kesslers_sim = sum.(CSD_collection2)

xaxis_popsizes = 1:Nd
sigma_inv = b_s1/(b_s1-d_s1)
sigma = (b_s1-d_s1)/b_s1
numres = nu*Nd*sigma_inv
numres0 = nu*Nd
kappas_pred = [numres/kappa^2*(1-1/kappa)^(numres-1) for kappa in xaxis_popsizes]
Lan = Landau(1.0)
kesslers_pred = [1/numres*pdf(Lan, m/numres - log(numres0) ) for m in xaxis_popsizes]

plot_subpopsizes(kappas_sim, kesslers_sim, xaxis_popsizes, kappas_pred, kesslers_pred)

### First versus second clone
firsts = [collection[2] for collection in CSD_collection]
seconds = [collection[3] for collection in CSD_collection]

difference = firsts.-seconds
plot_sizedifference_simulated(difference)




