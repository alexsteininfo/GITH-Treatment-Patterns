# title: Create phylogenetic trees to analyze resistance dynamics
# author: alexander stein

#include("/Users/stein02/Desktop/github/TreeStatistics.jl/src/SomaticEvolution.jl")
#using .SomaticEvolution
#using Revise
# Comment: Use SomaticEvolution.function to use functionalities of SomaitcEvolution to avoid conflicts

using SomaticEvolution
using Random, AbstractTrees
using Distributions: DiscreteNonParametric, pdf
using Statistics: mean, std, quantile!
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

# Parameters before treatment
b_s1 = 1.0
d_s1 = 0.0
b_r1 = 1.0
d_r1 = 0.0

# Parameters during treatment
b_s2 = 1.0
d_s2 = 3.0
b_r2 = 1.0
d_r2 = 0.0

# rate of neutral mutations
mu = 2.0
omega = 2*mu
# rate of resistance mutations
nu = 1e-3
# resistance mutation distribution
dist1 = DiscreteNonParametric([0.0],[1.0])
dist2 = DiscreteNonParametric([0.0],[1.0])
# Maximum number of subclones
maxclones = 1000

# other parameters
Nd = Int(1e4)
Nmax_vector = Nd*[1.2]
Nf = Int(Nmax_vector[end])
measurements = length(Nmax_vector) + 1

timevector = Float64[]

realizations = 5

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

function plot_SFS_loglog(SFS_xaxis, SFS_single1, SFS_single2)
    
    f = Figure()
    ax = Axis(f[1,1], xscale=log10, yscale=log10,
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Number of mutations, 𝑆ₖ", ylabelsize = 30, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    SFS_single1_ad = [max(Sk,0.5) for Sk in SFS_single1]
    SFS_single2_ad = [max(Sk,0.5) for Sk in SFS_single2]

    # Plot a single trajectory of the end state
    scatter!(ax, SFS_xaxis, SFS_single1_ad, color = ColPaired[9], markersize = 15, label = "homog. treatment response")
    scatter!(ax, SFS_xaxis, SFS_single2_ad, color = ColPaired[6], markersize = 12, marker = :x, label = "heterog. treatment response")

    ylims!(0.6,max(SFS_single1[1],SFS_single2[1])*1.2)

    axislegend(ax, labelsize=25, framecolor = :white, position = :rt)
    
    save("figures/Resis/SFS_resisVSnorm.png", f)
end

function plot_SFS_linear(SFS_xaxis, SFS_single1, SFS_simu_single2)
    
    f = Figure()
    ax = Axis(f[1,1], #xscale=log10, yscale=log10,
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Number of mutations, 𝑆ₖ", ylabelsize = 30, yticklabelsize = 25,
        xticks = ([0,10000,20000]),
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    # Plot a single trajectory of the end state
    scatter!(ax, SFS_xaxis, SFS_single1, color = ColPaired[9], markersize = 18, label = "SFS 1")
    scatter!(ax, SFS_xaxis, SFS_single2, color = ColPaired[1], markersize = 18, label = "SFS 1")

    axislegend(ax, labelsize=25, framecolor = :white, position = :rt)
    
    save("figures/Resis/SFS_resis_inf_linear.png", f)
end



############
### Main ###
############


SFS_collection_res = Any[[] for j in 1:measurements]

for i in 1:realizations
    # Set a random variable 
    rng = Random.seed!(i)

    # Grow population to detection size
    input = BranchingInput(Nmax=Nd, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
    selection = SelectionDistribution(dist1, nu, maxclones)
    global simulation = runsimulation(SimpleTreeCell, input, selection, rng)
    
    # Extract and save distributions
    #randomnode = simulation.output.cells[1]
    randomnode = simulation.output.singlemodule[1]
    root = getroot(randomnode)
    #global t_detec = age(simulation.output)

    modules = simulation.output.singlemodule

    ### SFS
    # after phase 1, all mutations occured before t_detec
    # thus, sitefreq_pre will be full and sitefreq_new will be empty
    #sitefreqs_pre, sitefreqs_new = get_sitefreqs_selection(root)
    #push!(SFS_collection_pre[1], sitefreqs_pre)
    #push!(SFS_collection_new[1], sitefreqs_new)

    sitefreqs = get_sitefreqs(root)
    push!(SFS_collection_res[1], sitefreqs)


    # Simlate the decrease

    for (j, Nmax) in enumerate(Nmax_vector)
        # set a random seed
        rng = Random.seed!(Int32(i+j*Nf))

        # simulate 1 more time unit
        input = BranchingInput(Nmax=Nmax, birthrate=b_s2, deathrate=d_s2, μ=mu, ploidy=1)
        selection = SelectionDistribution(dist2, nu, maxclones)
        global simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input, selection, rng)

        # Find root of the tree
        randomnode = simulation.output.singlemodule[1]
        root = getroot(randomnode)

        modules = simulation.output.singlemodule
        
        # SFS
        #sitefreqs_pre, sitefreqs_new = get_sitefreqs_selection(root)
        #push!(SFS_collection_pre[j+1], sitefreqs_pre)
        #push!(SFS_collection_new[j+1], sitefreqs_new)

        sitefreqs = get_sitefreqs(root)
        push!(SFS_collection_res[j+1], sitefreqs)

    end
end


SFS_collection_sen = []

for i in 1:realizations
    # Set a random variable 
    rng = Random.seed!(realizations+i)

    # Grow population to detection size
    input = BranchingInput(Nmax=Nf, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
    global simulation = runsimulation(SimpleTreeCell, input, rng)
    
    # Extract and save distributions
    #randomnode = simulation.output.cells[1]
    randomnode = simulation.output.singlemodule[1]
    root = getroot(randomnode)
    #global t_detec = age(simulation.output)

    modules = simulation.output.singlemodule

    ### SFS
    sitefreqs = get_sitefreqs(root)
    push!(SFS_collection_sen, sitefreqs)
end


###############################
### Site frequency spectrum ###
###############################

### We have a common xaxis
SFS_xaxis = 1:Nf

### Create histogram from simulation data
# Initilization
SFS_simu_res = Any[Any[] for it in 1:measurements]
# Iterate over measurements
for (m,SFS_subcollection) in enumerate(SFS_collection_res)
    # Iterae over realizations
    for sitefreqs in SFS_subcollection
        SFS = create_histo(sitefreqs,SFS_xaxis)
        push!(SFS_simu_res[m], SFS)
    end
end


### Within one measurement, group together according to frequencies k
# Initilization
sortedSFS_simu_res = [[Int32[0 for i in 1:realizations] for k in 1:Nmax_vector[end]] for j in 1:measurements]
# Iterate over measurements
for (j,SFS_subsimu) in enumerate(SFS_simu_res)
    for (i,SFS) in enumerate(SFS_subsimu)
        for (k,Sk) in enumerate(SFS)
            sortedSFS_simu_res[j][k][i] = Sk
        end
    end
end

# Create histogram from simulation data
SFS_simu_sen = []
for sitefreqs in SFS_collection_sen
    SFS = create_histo(sitefreqs, SFS_xaxis)
    push!(SFS_simu_sen, SFS)
end
# Group together according to frequencies k
sortedSFS_simu_sen = [Int32[0 for i in 1:realizations] for k in 1:Nf]
for (i,SFS) in enumerate(SFS_simu_sen)
    for (k,Sk) in enumerate(SFS)
        sortedSFS_simu_sen[k][i] = Sk
    end
end


### Make the plots
plot_SFS_loglog(SFS_xaxis, SFS_simu_sen[1], SFS_simu_res[2][1])
#plot_SFS_linear(SFS_xaxis, sortedSFS_simu_tot, SFS_simu_single, SFS_pred_tot, CSD_collection)
#plot_SFS_residuals(SFS_xaxis, sortedSFS_simu_tot, SFS_simu_single, SFS_pred_tot, CSD_collection)

#plot_SFS_pres(SFS_xaxis, sortedSFS_simu_pre, SFS_pred_pre, sortedSFS_simu_new, SFS_pred_new, sortedSFS_simu_tot, SFS_pred_tot)


