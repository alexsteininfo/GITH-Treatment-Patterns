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
d_s2 = 1.5
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
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Number of mutations, 𝑆ₖ", ylabelsize = 30, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    #xlims!(ax, 0.8, 250)
    #ylims!(ax,1/realizations/100, maximum(SFS_means[end])*2.0)
    ylims!(ax,0.5, maximum(SFS_means[end])*2.0)

    # Plot a single trajectory of the end state
    scatter!(ax, SFS_xaxis, SFS_simu_single[end], color = ColPaired[9], markersize = 18, label = "single simulation")

    # Plot the expected SFS of the end state
    lines!(ax, SFS_xaxis, SFS_means[end], color=ColPaired[2], linewidth=7.0, label = "avg. simulations")
    #lines!(ax, SFS_histo_x, SFS_histo_y, color=pastelorange, linewidth=4.0)

    # Plot the predicted SFS
    #SFS_xaxis_new = Float64[i for i in SFS_xaxis]
    #lines!(ax, SFS_xaxis, SFS_pred_tot, color=:black, linewidth=6.0, linestyle=:dot)
    lines!(ax, SFS_xaxis, SFS_pred_tot, color=:black, linewidth=7.0, linestyle = (:dot, 1.5), label="prediction")

    axislegend(ax, labelsize=25, framecolor = :white, position = :rt)
    
    save("figures/Resis/SFS_resis_inf_loglog.png", f)
end

function plot_SFS_linear(SFS_xaxis, SFS_simu, SFS_simu_single, SFS_pred, CSD_collection)
    
    # Simulations
    SFS_means = [mean.(Sk) for Sk in SFS_simu] # iterate over measurements and take mean over realizations
    SFS_means = [[maximum([1e-3,f]) for f in SFS] for SFS in SFS_means ]

    #SFS_histo_x, SFS_histo_y = ScatterToHisto(SFS_xaxis, SFS_means[end], 500)
    #SFS_histo_y = [maximum([1e-3,f]) for f in SFS_histo_y ]


    SFS_pred = [[maximum([0.1,f]) for f in SFS] for SFS in SFS_pred ]

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], #xscale=log10, yscale=log10,
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Number of mutations, 𝑆ₖ", ylabelsize = 30, yticklabelsize = 25,
        xticks = ([0,10000,20000]),
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    #xlims!(ax, 0.8, 250)
    ylims!(ax,0, 50)
    #ylims!(ax,0.5, maximum(SFS_means[end])*2.0)

    scatter!(ax, SFS_xaxis, SFS_simu_single[end], color=darkcyan, label = "SFS")
    vlines!(ax, CSD_collection[1], label = "CSD")

    axislegend(ax, labelsize=25, framecolor = :white, position = :rt)
    
    save("figures/Resis/SFS_resis_inf_linear.png", f)
end

function plot_SFS_residuals(SFS_xaxis, SFS_simu, SFS_simu_single, SFS_pred, CSD_collection)
    
    # Simulations
    SFS_means = [mean.(Sk) for Sk in SFS_simu] # iterate over measurements and take mean over realizations
    SFS_means = [[maximum([1e-3,f]) for f in SFS] for SFS in SFS_means ]

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], #xscale=log10, yscale=log10,
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "S_k-E[S_k], 𝑆ₖ", ylabelsize = 30, yticklabelsize = 25,
        xticks = ([0,10000,20000]),
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    #xlims!(ax, 0.8, 250)
    ylims!(ax,-1, 50)
    #ylims!(ax,0.5, maximum(SFS_means[end])*2.0)

    residuals = SFS_simu_single[end] - SFS_pred

    scatter!(ax, SFS_xaxis, residuals, color=darkcyan, label = "SFS")
    #vlines!(ax, CSD_collection[1], label = "CSD")

    axislegend(ax, labelsize=25, framecolor = :white, position = :rt)
    
    save("figures/Resis/SFS_resis_resid.png", f)
end

function plot_scMBdist(MB_xaxis, sMB_yaxis_series, MB_yaxis_single, MB_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Number mutations, 𝘫", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Number cells, 𝑀ⱼ", ylabelsize = 30, yticklabelsize = 25
    )

    # Simulations
    MB_mean = mean.(sMB_yaxis_series)
    MB_sigma = std.(sMB_yaxis_series)

    #println("Mean of sMB = ", sum(MB_xaxis.*MB_mean)/sum(MB_mean))
    #println("total number of divisions = ", sum(MB_yaxis_single))

    scatter!(ax, MB_xaxis, MB_yaxis_single, color = ColPaired[9], markersize = 15, label = "single simulation")

    #=
    lines!(ax, MB_xaxis, MB_mean, color=ColPaired[2], linewidth=6.0, label = "avg. simulations")
    lowerband = MB_mean-MB_sigma
    lowerband = [maximum([x,0.5]) for x in lowerband]
    upperband = MB_mean+MB_sigma
    upperband = [maximum([x,0.5]) for x in upperband]
    band!(ax, MB_xaxis, lowerband, upperband, color = (darkcyan, 0.2))

    lines!(ax, MB_xaxis, MB_pred, color=:black, linewidth=7.0, linestyle = (:dot,1.5), label = "prediction")
    =#

    save("figures/Resis/scMBdist_resis_inf.png", f)
end

function plot_peaksovertime(h_peaks, h_peaks2, t_tilde, t_predicted)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Peak height, 𝘩", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Arrival time, 𝘵", ylabelsize = 30, yticklabelsize = 25
    )

    scatter!(ax, h_peaks, t_tilde, color=:black, linewidth=6.0, linestyle = :dot, label="using 𝘩(𝘬)")
    scatter!(ax, h_peaks2, t_tilde, color=:red, linewidth=6.0, linestyle = :dot, label="using 𝘩(𝜅)")
    lines!(ax, h_peaks2, t_predicted, color=pastelorange, linewidth=4.0, label="prediction")

    #xlims!(ax, 0.9*minimum(h_peaks), 1.1*maximum(h_peaks))
    #ylims!(ax, 2.0, 15.0)
    axislegend(ax, labelsize=25, framecolor = :white, position = :rb)

    save("figures/Resis/peakanalysis.png", f)
end

function plot_kappaVSk(k_tilde, kappa_tilde)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Real clone size, 𝜅", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Pseudo clone size, 𝘬", ylabelsize = 30, yticklabelsize = 25
    )

    scatter!(ax, kappa_tilde, k_tilde, color=:black, linewidth=6.0, linestyle = :dot)
    lines!(ax, kappa_tilde, kappa_tilde, color=pastelorange, linewidth=4.0, label="prediction")

    #xlims!(ax, 0.9*minimum(h_peaks), 1.1*maximum(h_peaks))
    #ylims!(ax, 2.0, 15.0)
    axislegend(ax, labelsize=25, framecolor = :white, position = :rb)

    save("figures/Resis/kVSkappa.png", f)
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

SFS_collection_sen = Any[[] for j in 1:measurements]
SFS_collection_res = Any[[] for j in 1:measurements]

SFS_collection_tot = Any[[] for j in 1:measurements]

# Analysis of single trajectories
#kappa_tilde = Int32[]   # data taken from clonesizes
t_tilde = Float64[]

k_tilde = Int32[]       # data taken from the SFS
h_peaks = Float64[]

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
    push!(SFS_collection_pre[1], sitefreqs_pre)
    push!(SFS_collection_new[1], sitefreqs_new)

    #sitefreqs = getallelefreq(simulation)

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
        # set a random seed
        rng = Random.seed!(Int32(i+j*Nf))

        # simulate 1 more time unit
        input = BranchingInput(Nmax=Nmax, birthrate=b_s2, deathrate=d_s2, μ=mu, ploidy=1)
        selection = SelectionDistribution(dist2, nu, maxclones)
        global simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input, selection, rng)

        # Find root of the tree
        randomnode = simulation.output.singlemodule[1]
        root = getroot(randomnode)

        i == 1 ? push!(timevector, age(simulation.output)) : nothing
        modules = simulation.output.singlemodule
        #println(treebreadth(root))
        #println(getsubclonesizes(simulation))
        
        # SFS
        #sitefreqs_pre, sitefreqs_new = get_sitefreqs(root, t_detec)
        sitefreqs_pre, sitefreqs_new = get_sitefreqs_selection(root)
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

        # clone size information
        if(j == measurements-1)
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

    end
end


function create_variantsizes(typesizes, ancestry)
    clone = length(typesizes)
    variantsizes = copy(typesizes)
    for anc in reverse(ancestry)
        variantsizes[anc] += variantsizes[clone]
        clone -= 1
    end
    return variantsizes
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

SFS_simu_tot = SFS_simu_pre .+ SFS_simu_new

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

# expected time 
ResPop = Nd-SenPop
tf = log(Nf/ResPop)/(b_r2-d_r2)

# Compute predicted CSD before and after treatment
CSD_initial = [predict_SFS_ExpGrowth(x,b_s1,d_s1,Nd,nu) for x in 1:Nd]
CSD_new = [predict_CSD_new(k, Nd, Nf, b_s2, d_s2, b_r2, d_r2, nu, tf) for k in 1:Nf]
CSD_pre = predict_CSD_pre(CSD_initial, Nd, Nf, b_r2, d_r2, tf)


#CSD_tot = CSD_new .+ CSD_pre
CSD_tot = CSD_pre
#CSD_tot = CSD_new 
#CSD_tot = [predict_SFS_ExpGrowth(x,b_s1,d_s1,Nf,nu) for x in 1:Nf]

# Then compute the neutral tails of resistant subclones
SFS_neutral = predict_SFS_resis(CSD_tot, Nf, b_r1, d_r1, omega)
#SFS_peaks = predict_SFS_resis2(CSD_tot, t_detec, tf, b_s1, b_s2, b_r1, d_r1, Nf)
#SFS_pred_tot = SFS_neutral .+ SFS_peaks
#SFS_pred_tot = SFS_neutral

#SFS_pred_tot = [omega*Nf/(k^2+k) for k in SFS_xaxis]

SFS_pred_tot = [predict_SFS_ExpGrowth(x,b_r1,d_r1,Nf,omega) for x in 1:Nf]

### Make the plots
plot_SFS_loglog(SFS_xaxis, sortedSFS_simu_tot, SFS_simu_single, SFS_pred_tot)
plot_SFS_linear(SFS_xaxis, sortedSFS_simu_tot, SFS_simu_single, SFS_pred_tot, CSD_collection)
plot_SFS_residuals(SFS_xaxis, sortedSFS_simu_tot, SFS_simu_single, SFS_pred_tot, CSD_collection)

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
#println("Average mutation count for first resistant cell = ", mu*predict_averageGeneration(b_s1,d_s1,Nmut))

#scMB_compPoiss = Nmax_vector[end]*dist_compound_poisson(scMB_xaxis, genAverage+2*b_s2*tf, mu)
#scMB_pred = predict_scMBdist(scMB_xaxis, genAverage+2*b_s2*age(simulation), mu, Nmax_vector[end])
scMB_pred = predict_scMBdist(scMB_xaxis, genAverage+2*b_s2*(age(simulation)-t_detec), mu, Nmax_vector[end])


# Make the plot
plot_scMBdist(scMB_xaxis, sortedscMB_simu, scMB_simu[1], scMB_pred)



#####################################
### Analysing single trajectories ###
#####################################

# Subjects of analysis are
# kappa_tilde   # read out from simulation
# t_tilde

# k_tilde       # inferred from SFS
# h_peaks

function createallpossiblesums(vector, threshold)
    summands = []
    v1 = vector[1]
    for v2 in vector
        sum = v1+v2
        push!(summands, sum)
    end
    return summands
end

# compute a corrected clone size distribution and kappa_tilde
kappa_tilde = Int32[]
CSD_collection_corrected = Any[]
for (CSD,ances) in zip(CSD_collection, ancestry_collection)
    CSD_corrected = create_variantsizes(CSD, ances[2:end])
    push!(CSD_collection_corrected, CSD_corrected)
    push!(kappa_tilde, maximum(CSD_corrected[2:end]))
end

# now make a pure inference from the SFS
k_tilde = Int32[]       # data taken from the SFS
h_peaks = Int32[]

for SFS in SFS_simu_tot[end]
    threshold = 8.0
    index = findlast(item -> item > threshold, SFS)
    push!(k_tilde, index)
    push!(h_peaks,SFS[index])
end

# now measure peaks at kappa_tilde
h_peaks2 = Int32[]
for (SFS, kappa) in zip(SFS_simu_tot[end], kappa_tilde)
    push!(h_peaks2, SFS[kappa])
end

h_peaks3 = Int32[]
for (SFS, kappa) in zip(SFS_simu_tot[end], kappa_tilde)
    maximum = kappa+50
    cum_height = sum(SFS[kappa:end])
    push!(h_peaks3, cum_height)
end

t_predicted = [(h/mu + 2)/(2*b_s1) for h in h_peaks3]

plot_kappaVSk(kappa_tilde, k_tilde)
plot_peaksovertime(h_peaks, h_peaks3, t_tilde, t_predicted)


#clonesizes2 = CSD_collection[1]
#ancestry2 = ancestry_collection[1]
#realclonesizes2 = create_variantsizes(clonesizes2, ancestry2[2:end])


#for i in 1:100
#    #println(k_tilde[i], " and ", kappa_tilde[i])
#    println(k_tilde[i]-kappa_tilde[i])
#end

#for subclone in simulation.output.subclones
#    println(subclone.size)
#end


