# title: Create phylogenetic trees to analyze resistance dynamics
# author: alexander stein

using SomaticEvolution, Random, AbstractTrees
using Statistics: mean, std, quantile!
using Distributions: Poisson, Binomial, pdf
using SpecialFunctions: loggamma    # for the Consul Poisson distribution
using CairoMakie
using Colors
using Integrals

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

# Parameters before treatment
b_s1 = 1.0
d_s1 = 0.0
b_r1 = 1.0
d_r1 = 0.0

# Parameters during treatment
b_s2 = 1.0
d_s2 = 1.2

# b_r2 and d_r2 are implemented via selection
#b_r2 = 1.0
#d_r2 = 0.0

# rate of neutral mutations
mu = 4.0
omega = 2*mu

# Define input variables for simulations
mutant_time = [8.0, Inf]
mutant_selection = [0.0, 0.0]

mutant_time2 = [Inf, Inf]
mutant_selection2 = [0.0, 0.0]

numclones = length(mutant_time)+1

# other parametersxn
Nd = 10000
Nmax_vector = Nd*[1.1, 1.2, 1.3, 1.4, 1.5]
measurements = length(Nmax_vector) + 1

timestep = Inf
timevector = Float32[]

realizations = 20

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

function plot_SFS_loglog(SFS_xaxis, SFS_simu, SFS_simu_single, SFS_pred)
    
    # Simulations
    SFS_means = [mean.(Sk) for Sk in SFS_simu] # iterate over measurements and take mean over realizations
    SFS_means = [[maximum([1e-3,f]) for f in SFS] for SFS in SFS_means ]

    #SFS_histo_x, SFS_histo_y = ScatterToHisto(SFS_xaxis, SFS_means[end], 500)
    #SFS_histo_y = [maximum([1e-3,f]) for f in SFS_histo_y ]

    SFS_pred = [[maximum([0.1,f]) for f in SFS] for SFS in SFS_pred ]

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], xscale=log10, yscale=log10,
        xlabel = "# of cells, 𝑘", xlabelsize = 25, xticklabelsize = 25,
        ylabel = "exp. # of mutations, E[𝑆ₖ]", ylabelsize = 25, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    #xlims!(ax, 0.8, 250)
    #ylims!(ax,1/realizations/100, maximum(SFS_means[end])*2.0)
    ylims!(ax,0.5, maximum(SFS_means[end])*2.0)

    # Plot the expected SFS of the end state
    lines!(ax, SFS_xaxis, SFS_means[end], color=pastelorange, linewidth=4.0)
    #lines!(ax, SFS_histo_x, SFS_histo_y, color=pastelorange, linewidth=4.0)

    # Plot a single trajectory of the end state
    scatter!(ax, SFS_xaxis, SFS_simu_single[end], color=darkcyan)
    
    save("figures/Resis/SFS_resis_loglog.png", f)
end

function plot_SFS_linear(SFS_xaxis, SFS_simu, SFS_simu_single, SFS_pred)
    
    # Simulations
    SFS_means = [mean.(Sk) for Sk in SFS_simu] # iterate over measurements and take mean over realizations
    SFS_means = [[maximum([1e-3,f]) for f in SFS] for SFS in SFS_means ]

    #SFS_histo_x, SFS_histo_y = ScatterToHisto(SFS_xaxis, SFS_means[end], 500)
    #SFS_histo_y = [maximum([1e-3,f]) for f in SFS_histo_y ]


    SFS_pred = [[maximum([0.1,f]) for f in SFS] for SFS in SFS_pred ]

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], #xscale=log10, yscale=log10,
        xlabel = "# of cells, 𝑘", xlabelsize = 25, xticklabelsize = 25,
        ylabel = "exp. # of mutations, E[𝑆ₖ]", ylabelsize = 25, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    #xlims!(ax, 0.8, 250)
    ylims!(ax,0, 100)
    #ylims!(ax,0.5, maximum(SFS_means[end])*2.0)

    #lines!(ax, SFS_xaxis, SFS_means[1], color=darkcyan, linewidth=6.0)
    #scatter!(ax, SFS_xaxis, SFS_means[1], color=darkcyan)
    for m in 1:measurements
        #lines!(ax, SFS_xaxis, SFS_means[m], color=pastelorange, linewidth=4.0)
        if(m==0)
            scatter!(ax, SFS_xaxis, SFS_simu_single[m], color=:black)
        else
            scatter!(ax, SFS_xaxis, SFS_simu_single[m], color=(:black,m*0.15))
        end
        
        #lines!(ax, SFS_histo_x, SFS_histo_y, color=pastelorange, linewidth=4.0)
        #scatter!(ax, SFS_xaxis, SFS_means[m], color=pastelorange)
        #lines!(ax, SFS_xaxis, SFS_pred[m], color=pastelorange, linewidth=6.0, linestyle = :dot)
    end

    # Plot a single trajectory of the end state
    #scatter!(ax, SFS_xaxis, SFS_simu_single, color=darkcyan)
    
    save("figures/Resis/SFS_resis_linear.png", f)
end

function plot_scMBdist(MB_xaxis, sMB_yaxis_series, MB_yaxis_single, MB_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "# mutations, 𝘫", xlabelsize = 25, xticklabelsize = 25,
    ylabel = "# cells with 𝘫 mutations, 𝑀ⱼ", ylabelsize = 25, yticklabelsize = 25
    )

    # Simulations
    MB_mean = mean.(sMB_yaxis_series)
    MB_sigma = std.(sMB_yaxis_series)

    #println("Mean of sMB = ", sum(MB_xaxis.*MB_mean)/sum(MB_mean))
    #println("total number of divisions = ", sum(MB_yaxis_single))

    lines!(ax, MB_xaxis, MB_mean, color=darkcyan, linewidth=6.0)
    lowerband = MB_mean-MB_sigma
    lowerband = [maximum([x,0.5]) for x in lowerband]
    upperband = MB_mean+MB_sigma
    upperband = [maximum([x,0.5]) for x in upperband]
    band!(ax, MB_xaxis, lowerband, upperband, color = (darkcyan, 0.2))

    scatter!(ax, MB_xaxis, MB_yaxis_single, color = ferngreen, label = nothing, markersize = 15)

    lines!(ax, MB_xaxis, MB_pred, color=pastelorange, linewidth=6.0, linestyle = :dot)

    save("figures/Resis/scMBdist.png", f)
end


############
### Main ###
############

# From the phylogenetic trees we obtain 
#   SFS_collection_pre, SFS_collection_new
#   tMB_collection_pre, tMB_collection_new
#   sMB_collection_pre, sMB_collection_new

scMB_collection = Any[[] for j in 1:measurements]

tMB_collection_pre = Any[Float64[] for j in 1:measurements]
tMB_collection_new = Any[Float64[] for j in 1:measurements]
tMB_collection_tot = Any[Float64[] for j in 1:measurements]


SFS_collection_pre = Any[[] for j in 1:measurements]
SFS_collection_new = Any[[] for j in 1:measurements]
SFS_collection_tot = Any[[] for j in 1:measurements]


for i in 1:realizations
    # Grow to maximal population size
    #input = BranchingInput(Nmax=Nd, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
    #selection = SelectionDistribution(Exponential(0.2), 0.2, 50)
    #selection = SelectionPredefined(1.0,1.0)
    #rng = Random.seed!(i)
    #simulation = runsimulation(SimpleTreeCell, input, selection, rng)

    rng = Random.seed!(i)
    input = BranchingInput(Nmax=Nd, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
    selection = SelectionPredefined(mutant_selection, mutant_time)
    simulation = runsimulation(SimpleTreeCell, input, selection, rng)
    
    # Extract and save distributions
    #randomnode = simulation.output.cells[1]
    randomnode = simulation.output.singlemodule[1]
    root = getroot(randomnode)
    global t_detec = age(simulation.output)

    i == 1 ? push!(timevector, t_detec) : nothing
    #println("Subclone sizes:")
    modules = simulation.output.singlemodule
    #println(getsubclonesizes(simulation))

    # SFS
    sitefreqs_pre, sitefreqs_new = get_sitefreqs(root, t_detec)
    # after phase 1, all mutations occured before t_detec
    # thus, sitefreq_pre will be full and sitefreq_new will be empty
    push!(SFS_collection_pre[1], sitefreqs_pre)
    push!(SFS_collection_new[1], sitefreqs_new)

    # tMB
    tMB_pre = length(sitefreqs_pre)
    tMB_new = length(sitefreqs_new)
    push!(tMB_collection_pre[1], tMB_pre)
    push!(tMB_collection_new[1], tMB_new)
    
    # scMB
    mb_dist = get_scMBs(root)
    push!(scMB_collection[1], mb_dist)

    # Simlate the decrease

    for (j, Nmax) in enumerate(Nmax_vector)
        # simulate 1 more time unit
        #input = BranchingInput(Nmax=Nmax, birthrate=b_s2, deathrate=d_s2, μ=mu, ploidy=1)
        #simulation = runsimulation(SimpleTreeCell, input)
        #simulation = runsimulation(SimpleTreeCell, WellMixed, simulation.output, input)
        #simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input)

        input = BranchingInput(Nmax=Nmax, tmax=t_detec+j*timestep, birthrate=b_s2, deathrate=d_s2, μ=mu, ploidy=1)
        selection = SelectionPredefined(mutant_selection2, mutant_time2)
        simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input, selection)

        # Find root of the tree
        #randomnode = simulation.output.cells[1]
        randomnode = simulation.output.singlemodule[1]
        root = getroot(randomnode)

        i == 1 ? push!(timevector, age(simulation.output)) : nothing
        modules = simulation.output.singlemodule
        #println(getsubclonesizes(simulation))
        
        # SFS
        sitefreqs_pre, sitefreqs_new = get_sitefreqs(root, t_detec)
        push!(SFS_collection_pre[j+1], sitefreqs_pre)
        push!(SFS_collection_new[j+1], sitefreqs_new)

        # tMB
        tMB_pre = length(sitefreqs_pre)
        tMB_new = length(sitefreqs_new)
        push!(tMB_collection_pre[j+1], tMB_pre)
        push!(tMB_collection_new[j+1], tMB_new)

        # scMB
        mb_dist = get_scMBs(root)
        push!(scMB_collection[j+1], mb_dist)
    end
end


###############################
### Site frequency spectrum ###
###############################

### We have a common xaxis
SFS_xaxis = 1:Int(Nmax_vector[end])

### Create histogram from simulation data
# Initilization
SFS_simu_pre = Any[Any[] for it in 1:measurements]
# Iterate over measurements
for (m,SFS_subcollection) in enumerate(SFS_collection_pre)
    # Iterae over realizations
    for sitefreqs in SFS_subcollection
        SFS = create_histo(sitefreqs,SFS_xaxis)
        push!(SFS_simu_pre[m], SFS)
    end
end

# Initilization
SFS_simu_new = Any[Any[] for it in 1:measurements]
# Iterate over measurements
for (m,SFS_subcollection) in enumerate(SFS_collection_new)
    # Iterae over realizations
    for sitefreqs in SFS_subcollection
        SFS = create_histo(sitefreqs,SFS_xaxis)
        push!(SFS_simu_new[m], SFS)
    end
end

# Create a vector of SFS_yaxis for single simulations
SFS_simu_single = []
for m in 1:measurements
    SFS = SFS_simu_pre[m][1] .+ SFS_simu_new[m][1]
    push!(SFS_simu_single, SFS)
end

### Within one measurement, group together according to frequencies k
# Initilization
sortedSFS_simu_pre = [[Int32[0 for i in 1:realizations] for k in 1:Nmax_vector[end]] for j in 1:measurements]
# Iterate over measurements
for (j,SFS_subsimu) in enumerate(SFS_simu_pre)
    for (i,SFS) in enumerate(SFS_subsimu)
        for (k,Sk) in enumerate(SFS)
            sortedSFS_simu_pre[j][k][i] = Sk
        end
    end
end

# Initilization
sortedSFS_simu_new = [[Int32[0 for i in 1:realizations] for k in 1:Nmax_vector[end]] for j in 1:measurements]
# Iterate over measurements
for (j,SFS_subsimu) in enumerate(SFS_simu_new)
    for (i,SFS) in enumerate(SFS_subsimu)
        for (k,Sk) in enumerate(SFS)
            sortedSFS_simu_new[j][k][i] = Sk
        end
    end
end

# the total site frequency spectrum
sortedSFS_simu_tot = sortedSFS_simu_pre .+ sortedSFS_simu_new

### Prepare predictions

SFS_pred_pre = Any[]
SFS_pred_new = Any[]
SFS_pred_tot = Any[]

# Compute the prediction after exponential growth
    # after phase 1, all mutations occured before t_detec
    # thus, sitefreq_pre will be full and sitefreq_new will be empty
push!(SFS_pred_pre, [predict_SFS_ExpGrowth(x,b_s1,d_s1,Nd,omega) for x in SFS_xaxis])
push!(SFS_pred_new, [0 for x in SFS_xaxis])

# Compute the predictions after decreasing population
for Nmax in Nmax_vector
    # Determine deterministic timing
    tN = log(Nmax/Nd)/(b_s2-d_s2)
    # Preexisting mutations
    #SFS_yaxis_pre = predict_SFS_IncrPre(SFS_pred_pre[1], SFS_xaxis, b_s2, d_s2, tN)
    SFS_yaxis_pre = [1 for x in SFS_xaxis]
    push!(SFS_pred_pre, SFS_yaxis_pre)
    #SFS_yaxis_new = [predict_SFS_NumInt(x, Nd, Nmax, b_s2, d_s2, omega) for x in SFS_xaxis]
    SFS_yaxis_new = [1 for x in SFS_xaxis]
    push!(SFS_pred_new, SFS_yaxis_new)
end

SFS_pred_tot = SFS_pred_pre .+ SFS_pred_new

### Make the plots
plot_SFS_loglog(SFS_xaxis, sortedSFS_simu_tot, SFS_simu_single, SFS_pred_tot)
plot_SFS_linear(SFS_xaxis, sortedSFS_simu_tot, SFS_simu_single, SFS_pred_tot)

#plot_SFS_pres(SFS_xaxis, sortedSFS_simu_pre, SFS_pred_pre, sortedSFS_simu_new, SFS_pred_new, sortedSFS_simu_tot, SFS_pred_tot)


##################################################
### single cell mutational burden distribution ###
##################################################

### Let us analyse the distribution at the end first
scMB_collection_end = scMB_collection[end]

# We have a common xaxis
highestburden = maximum(maximum.(scMB_collection_end))
lowestburden = minimum(minimum.(scMB_collection_end))
scMB_xaxis = lowestburden:highestburden

# Create histogram from simulation data
scMB_simu = []
for scMB in scMB_collection_end
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
genAverage = predict_averageGeneration(b_s1,d_s1,Nd)
tf = log(Nmax_vector[end]/Nd)/(b_s2-d_s2)

#scMB_compPoiss = Nmax_vector[end]*dist_compound_poisson(scMB_xaxis, genAverage+2*b_s2*tf, mu)
scMB_pred = predict_scMBdist(scMB_xaxis, genAverage+2*b_s2*tf, mu, Nmax_vector[end])


# Make the plot
#plot_scMBdist(scMB_xaxis, sortedscMB_simu, scMB_simu[1], scMB_pred)

