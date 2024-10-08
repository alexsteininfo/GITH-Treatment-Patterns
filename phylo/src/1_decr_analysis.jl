# title: Create phylogenetic trees to analyze decreasing population
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

using ColorBrewer
ColPaired = palette("Paired", 12)

#########################
### Define parameters ###
#########################
b_1 = 1.0
d_1 = 0.0

b_2 = 1.0
d_2 = 3.0

mu = 2.0
omega = 2*mu

Nd = Int(1e5)
#Nmin_vector = Nd*[0.5,0.1]
Nmin_vector = Nd*[0.64,0.32,0.16,0.08,0.04,0.02,0.01]
Nmin_vector = floor.(Int, Nmin_vector)
measurements = length(Nmin_vector) + 1

realizations = 100

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

function plot_SFS_pres(SFS_xaxis, SFS_simu_pre, SFS_pred_pre, SFS_simu_new, SFS_pred_new, SFS_simu_tot, SFS_pred_tot)

    ### Prepare simulation data
    SFS_means_pre = [mean.(Sk) for Sk in SFS_simu_pre]
    SFS_means_new = [mean.(Sk) for Sk in SFS_simu_new]
    SFS_means_tot = [mean.(Sk) for Sk in SFS_simu_tot]

    ### The plot itself
    f = Figure()
    ax = Axis(f[1,1], xscale=log10, yscale=log10,
        xlabel = "Number of cells, 𝑘", xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Number of mutations, 𝑆ₖ", ylabelsize = 30, yticklabelsize = 25#,
        #xticks = ([1e0,1e1,1e2,1e3,1e4],[L"10^0",L"10^1",L"10^2",L"10^3",L"10^4"]),
        #yticks = ([1e0,1e1,1e2,1e3],[L"10^0",L"10^1",L"10^2",L"10^3"])
    )

    xlims!(ax, 0.8, 500)
    ylims!(ax, 0.5, maximum(SFS_means_tot[1])*2.0)

    lines!(ax, SFS_xaxis, SFS_means_pre[end], color=ColPaired[7], linewidth=8.0, label = "preexisting")
    lines!(ax, SFS_xaxis, SFS_pred_pre[end], color=:black, linewidth=7.0, linestyle = (:dot,1.5))

    lines!(ax, SFS_xaxis, SFS_means_new[end], color=ColPaired[3], linewidth=8.0, label = "newly emerging")
    lines!(ax, SFS_xaxis, SFS_pred_new[end], color=:black, linewidth=7.0, linestyle = (:dot,1.5))

    lines!(ax, SFS_xaxis, SFS_means_tot[end], color=ColPaired[1], linewidth=8.0, label = "total")
    lines!(ax, SFS_xaxis, SFS_pred_tot[end], color=:black, linewidth=7.0, linestyle = (:dot,1.5), label="precictions")

    #someotherprediction = [1e4*k^(-2.0) for k in SFS_xaxis]
    #lines!(ax, SFS_xaxis, someotherprediction, color=:red, linewidth=3.0, linestyle = (:dot,0.5), label="power law")
    #lines!(ax, SFS_xaxis, someotherprediction/10, color=:red, linewidth=3.0, linestyle = (:dot,0.5))
    #lines!(ax, SFS_xaxis, someotherprediction/100, color=:red, linewidth=3.0, linestyle = (:dot,0.5))

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/Decr/SFS_decr_pres.png", f)
end


function plot_SFS_all(SFS_xaxis, SFS_simu, SFS_pred)
    
    # Simulations
    SFS_means = [mean.(Sk) for Sk in SFS_simu] # iterate over measurements and take mean over realizations
    SFS_means = [[maximum([0.1,f]) for f in SFS] for SFS in SFS_means ]

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
    ylims!(ax, 0.5, maximum(SFS_means[1])*2.0)

    for m in 1:measurements
        lines!(ax, SFS_xaxis, SFS_means[m], color=darkcyan, linewidth=6.0)
        lines!(ax, SFS_xaxis, SFS_pred[m], color=pastelorange, linewidth=6.0, linestyle = :dot)
    end
    
    save("figures/Decr/SFS_decr_all.png", f)
end

function plot_SFS_pre(SFS_xaxis, SFS_simu, SFS_pred)
    
    # Simulations
    SFS_means = [mean.(Sk) for Sk in SFS_simu] # iterate over measurements and take mean over realizations
    SFS_means = [[maximum([0.1,f]) for f in SFS] for SFS in SFS_means ]

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
    ylims!(ax, 0.5, maximum(SFS_means[1])*2.0)

    lines!(ax, SFS_xaxis, SFS_means[1], color=darkcyan, linewidth=6.0, label = "simulations")
    lines!(ax, SFS_xaxis, SFS_pred[1], color=pastelorange, linewidth=6.0, linestyle = :dot, label = "prediction")

    for m in 2:measurements
        lines!(ax, SFS_xaxis, SFS_means[m], color=darkcyan, linewidth=6.0)
        lines!(ax, SFS_xaxis, SFS_pred[m], color=pastelorange, linewidth=6.0, linestyle = :dot)
    end

    axislegend(ax, labelsize=25, framecolor = :white)
    
    save("figures/Decr/SFS_decr_pre.png", f)
end

function plot_SFS_new(SFS_xaxis, SFS_simu, SFS_pred)
    
    # Simulations
    SFS_means = [mean.(Sk) for Sk in SFS_simu] # iterate over measurements and take mean over realizations
    SFS_means = [[maximum([0.1,f]) for f in SFS] for SFS in SFS_means ]

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
    ylims!(ax, 0.5, maximum(SFS_means[2])*2.0)

    for m in 1:measurements
        lines!(ax, SFS_xaxis, SFS_means[m], color=darkcyan, linewidth=6.0)
        lines!(ax, SFS_xaxis, SFS_pred[m], color=pastelorange, linewidth=6.0, linestyle = :dot)
    end
    
    save("figures/Decr/SFS_decr_new.png", f)
end

function plot_mutburdenovertime(tMB_simu_pre, tMB_simu_new, tMB_pred_pre, tMB_pred_new)

    # Prepare simulation data for the plot
    tMB_mean_pre = [mean(it) for it in tMB_simu_pre]
    tMB_mean_new = [mean(it) for it in tMB_simu_new]
    
    tMB_simu_all = tMB_simu_pre .+ tMB_simu_new
    tMB_mean_all = [mean(it) for it in tMB_simu_all]

    tMb_pred_all = tMB_pred_pre .+ tMB_pred_new

    # The plot itself
    f = Figure()
    ax = Axis(f[1,1],
        xlabel = rich("Treatment time, 𝑡", subscript("f", font = :italic),font = :regular), xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Total mutational burden, 𝘉", ylabelsize = 30, yticklabelsize = 25#,
        #xticks = ([0,25,50,75,100],[L"0",L"25",L"50",L"75",L"100"]),
        #yticks = ([0,200,400,600,800],[L"0",L"200",L"400",L"600",L"800"])
    )
    
    sizescale = [Nd; Nmin_vector]
    timescale = log.(sizescale/Nd)/(b_2-d_2)

    linearscale = [tMB_mean_all[1]*size/Nd for size in sizescale]

    # An exponentially growing population with the same size
    #lines!(timescale, linearscale, linewidth = 3.0, color=:black, linestyle=:dash)

    # Plot the predictions
    lines!(timescale, tMB_pred_pre, linewidth = 3.0, color=ColPaired[7], label = "preexisting")
    lines!(timescale, tMB_pred_new, linewidth = 3.0, color=ColPaired[3], label = "newly emerging")
    lines!(timescale, tMb_pred_all, linewidth = 3.0, color=ColPaired[1], label = "total")

    # Plot the simulation results
    scatter!(timescale, tMB_mean_pre, markersize = 18, color = ColPaired[8], marker = :utriangle)
    scatter!(timescale, tMB_mean_new, markersize = 18, color = ColPaired[4], marker = :dtriangle)
    scatter!(timescale, tMB_mean_all, markersize = 18, color = ColPaired[2])

    #Label(f[1, 1, TopLeft()], "L", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    #axislegend(ax, labelsize=32, framecolor = :white, position = :rc)
    axislegend(ax, labelsize=25, framecolor = :white)

    save("figures/Decr/tMB_time.png", f)
end

function plot_scMBdist(MB_xaxis, scMB_yaxis_series, MB_yaxis_single, MB_pred)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Number of mutations, 𝘫", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Number of cells, 𝑀ⱼ", ylabelsize = 30, yticklabelsize = 25
    )

    # Simulations
    MB_mean = mean.(scMB_yaxis_series)
    MB_sigma = std.(scMB_yaxis_series)

    #println("Mean of scMB = ", sum(MB_xaxis.*MB_mean)/sum(MB_mean))
    #println("total number of divisions = ", sum(MB_yaxis_single))

    scatter!(ax, MB_xaxis, MB_yaxis_single, color = ColPaired[9], markersize = 15, label = "single simulation")

    lines!(ax, MB_xaxis, MB_mean, color=ColPaired[2], linewidth=8.0, label = "avg. simulations")
    lowerband = MB_mean-MB_sigma
    lowerband = [maximum([x,0.0]) for x in lowerband]
    upperband = MB_mean+MB_sigma
    upperband = [maximum([x,0.0]) for x in upperband]
    band!(ax, MB_xaxis, lowerband, upperband, color = (ColPaired[2], 0.2))

    lines!(ax, MB_xaxis, MB_pred, color=:black, linewidth=7.0, linestyle = (:dot,1.5), label = "prediction")

    #axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figures/Decr/scMBdist.png", f)
end

function plot_exp_scMB(timescale, scMB_simu, scMB_pred)
    # The plot itself
    f = Figure()
    ax = Axis(f[1,1],
        xlabel = rich("Treatment time, 𝑡", subscript("f", font = :italic),font = :regular), xlabelsize = 30, xticklabelsize = 25,
        ylabel = "Average scMB,〈 𝑗 〉", ylabelsize = 30, yticklabelsize = 25,
        yticks = ([45, 50, 55])
        #xticks = ([0,25,50,75,100],[L"0",L"25",L"50",L"75",L"100"]),
        #yticks = ([0,200,400,600,800],[L"0",L"200",L"400",L"600",L"800"])
    )


    # Plot the simulation results
    scatter!(timescale, scMB_simu, markersize = 18, color = ColPaired[2], label = "simulations")
    # Plot the predictions
    lines!(timescale, scMB_pred, linewidth = 5.0, color=:black, linestyle = (:dot,1.5), label = "prediction")

    ylims!(ax, 42, 58)

    axislegend(ax, labelsize=25, framecolor = :white, position=:lt)

    save("figures/Decr/mutburden_time.png", f)
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

    save("figures/Decr/exponenthist.png", f)
end

############
### Main ###
############

# From the phylogenetic trees we obtain 
#   SFS_collection_pre, SFS_collection_new
#   tMB_collection_pre, tMB_collection_new
#   scMB_collection_pre, scMB_collection_new

scMB_collection = Any[[] for j in 1:measurements]

tMB_collection_pre = Any[Float64[] for j in 1:measurements]
tMB_collection_new = Any[Float64[] for j in 1:measurements]
tMB_collection_tot = Any[Float64[] for j in 1:measurements]


SFS_collection_pre = Any[[] for j in 1:measurements]
SFS_collection_new = Any[[] for j in 1:measurements]
SFS_collection_tot = Any[[] for j in 1:measurements]


for i in 1:realizations
    #println("Realization = ", i)
    # Grow to maximal population size
    input = BranchingInput(Nmax=Nd, birthrate=b_1, deathrate=d_1, μ=mu, ploidy=1)
    rng = Random.seed!(i)
    simulation = runsimulation(SimpleTreeCell, input, rng)
    
    # Extract and save distributions
    #randomnode = simulation.output.cells[1]
    randomnode = simulation.output.singlemodule[1]
    root = getroot(randomnode)
    global t_detec = age(simulation.output)

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

    for (j, Nmin) in enumerate(Nmin_vector)
        #println("Nmin = ", Nmin)
        # simulate 1 more time unit
        input = BranchingInput(Nmax=Nmin, birthrate=b_2, deathrate=d_2, μ=mu, ploidy=1)
        #simulation = runsimulation(SimpleTreeCell, input)
        #simulation = runsimulation(SimpleTreeCell, WellMixed, simulation.output, input)
        simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input)

        # Find root of the tree
        #randomnode = simulation.output.cells[1]
        randomnode = simulation.output.singlemodule[1]
        root = getroot(randomnode)
        
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
SFS_xaxis = 1:Nd

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


### Within one measurement, group together according to frequencies k
# Initilization
sortedSFS_simu_pre = [[Int32[0 for i in 1:realizations] for k in 1:Nd] for j in 1:measurements]
# Iterate over measurements
for (j,SFS_subsimu) in enumerate(SFS_simu_pre)
    for (i,SFS) in enumerate(SFS_subsimu)
        for (k,Sk) in enumerate(SFS)
            sortedSFS_simu_pre[j][k][i] = Sk
        end
    end
end

# Initilization
sortedSFS_simu_new = [[Int32[0 for i in 1:realizations] for k in 1:Nd] for j in 1:measurements]
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
push!(SFS_pred_pre, [predict_SFS_ExpGrowth(x,b_1,d_1,Nd,omega) for x in SFS_xaxis])
push!(SFS_pred_new, [0 for x in SFS_xaxis])

# Compute the predictions after decreasing population
for Nmin in Nmin_vector
    # Preexisting mutations
    SFS_yaxis_pre = predict_SFS_DecrPre(SFS_pred_pre[1], SFS_xaxis, b_2, d_2, Nd, Nmin)
    push!(SFS_pred_pre, SFS_yaxis_pre)
    # Newly arising mutations
    SFS_yaxis_new = [predict_SFS_NumInt(x, Nd, Nmin, b_2, d_2, omega) for x in SFS_xaxis]
    push!(SFS_pred_new, SFS_yaxis_new)
end

SFS_pred_tot = SFS_pred_pre .+ SFS_pred_new

### Make the plots

#plot_SFS_pre(SFS_xaxis, sortedSFS_simu_pre, SFS_pred_pre)
#plot_SFS_new(SFS_xaxis, sortedSFS_simu_new, SFS_pred_new)
#plot_SFS_all(SFS_xaxis, sortedSFS_simu_tot, SFS_pred_tot)


plot_SFS_pres(SFS_xaxis, sortedSFS_simu_pre, SFS_pred_pre, sortedSFS_simu_new, SFS_pred_new, sortedSFS_simu_tot, SFS_pred_tot)

###############################
### total mutational burden ###
###############################

# Make the predictions
tMB_pred_pre = Float32[predict_tMB_ExpGrowth(b_1,d_1,Nd,omega)]
tMB_pred_new = Float32[0.0]

for Nmin in Nmin_vector
    push!(tMB_pred_pre, predict_tMB_DecrPre(SFS_pred_pre[1],b_2,d_2,Nd,Nmin)) 
    push!(tMB_pred_new, predict_tMB_DecrNew(b_2,d_2,Nd,Nmin,omega))
end

# Make the plot
plot_mutburdenovertime(tMB_collection_pre, tMB_collection_new, tMB_pred_pre, tMB_pred_new)

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
genAverage = predict_averageGeneration(b_1,d_1,Nd)
tf = log(Nmin_vector[end]/Nd)/(b_2-d_2)

#scMB_compPoiss = Nmin_vector[end]*dist_compound_poisson(scMB_xaxis, genAverage2*b_2*tf, mu)
scMB_pred = predict_scMBdist(scMB_xaxis, genAverage+2*b_2*tf, mu, Nmin_vector[end])
# Make the plot
plot_scMBdist(scMB_xaxis, sortedscMB_simu, scMB_simu[1], scMB_pred)


### Now, compute the average scMB
scMB_average = []
# Iterate over measurement times
for (i,scMB) in enumerate(scMB_collection)
    if(i==1)
        averages = sum.(scMB)/Nd
    else
        averages = sum.(scMB)/Nmin_vector[i-1]
    end
    push!(scMB_average, averages) 
end
scMB_ExpAverage = sum.(scMB_average)/realizations

# Work out the prediction
sizescale = [Nd; Nmin_vector]
timescale = log.(sizescale/Nd)/(b_2-d_2)
scMB_ExpPred = [mu*(genAverage + 2*b_2*t_2) for t_2 in timescale]

plot_exp_scMB(timescale, scMB_ExpAverage, scMB_ExpPred)


#########################
### Fitting exponents ###
#########################

kmin = 100
#exponents = [get_exponent(SFS, kmin) for SFS in SFS_simu]
exponents = [get_exponent(SFS, kmin) for SFS in SFS_simu_pre[end]]

println("Mean exponent = ", mean(exponents))
println("Std exponent = ", std(exponents))

plot_exponents(exponents)

