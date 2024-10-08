# Plot size difference between first and second subclone

using CairoMakie
using Colors
using Distributions
using ColorBrewer
ColPaired = palette("Paired", 12)


function T1_random(b,d,nu)
    u = rand()
    r = b-d
    return 1/r*log(1-r/(b*nu)*log(1-u))
end

function T2_random(b,d,nu,t1)
    u = rand()
    r = b-d
    return 1/r*log(exp(r*t1)-r/(b*nu)*log(1-u))
end

function SizeDifference_random(b,d,nu,td)
    t1 = T1_random(b,d,nu)
    t2 = T2_random(b,d,nu,t1)

    t2 > td ? println("ATTENTION") : nothing
    
    W = Exponential((b-d)/b)
    W1 = rand(W)
    W2 = rand(W)
    return W1*exp((b-d)*(td-t1)) - W2*exp((b-d)*(td-t2))
end

b = 0.14
d = 0.13
nu = 1e-7
td = 2555.0

samplesize = 1000000

samples = [SizeDifference_random(b,d,nu,td) for i in 1:samplesize]


function plot_compare_linear(samples)
    sample_max = maximum(samples)
    sample_min = minimum(samples)
    samples_pos = samples[samples .> 0]
    samples_neg = samples[samples .< 0]
    
    #println(length(samples_pos))
    #println(length(samples_neg))
    
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Difference in size", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25, #yscale=log10#,
    xticks = ([-2e5, 0.0, 2e5, 4e5], [rich("-2×10", superscript("5")), "0", rich("2×10", superscript("5")), rich("4×10", superscript("5"))]),
    yticks = ([0, 5e-6, 1e-5, 1.5e-5, 2.0e-5], ["0", rich("5.0×10", superscript("-6")), rich("1.0×10", superscript("-5")), rich("1.5×10", superscript("-5")),rich("2.0×10", superscript("-5"))])
    )

    bins = [i*0.02*1e6 - 1e6 for i in 1:100]
    reds = [(ColPaired[6],0.7) for i in 1:49]
    blues = [(ColPaired[2],0.7) for i in 1:50]
    col_redblue = [reds;blues]
    hist!(ax, samples, bins=bins, normalization = :pdf, color = col_redblue)

    xlims!(-3e5,5e5)


    save("figures/compare/SizeDifference_linear.png", f)    
end

function plot_compare_log(samples)
    sample_max = maximum(samples)
    sample_min = minimum(samples)
    samples_pos = samples[samples .> 0]
    samples_neg = samples[samples .< 0]
    
    println(length(sample_max))
    println(length(sample_min))
    
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = "Difference in size", xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Density", ylabelsize = 30, yticklabelsize = 25, yscale=log10,
    xticks = ([-2e7, 0.0, 2e7, 4e7])
    )

    bins = [i*0.02*1e8 - 1e8 for i in 1:100]
    reds = [(ColPaired[6],0.7) for i in 1:49]
    blues = [(ColPaired[2],0.7) for i in 1:50]
    col_redblue = [reds;blues]
    hist!(ax, samples, bins=bins, normalization = :pdf, color = col_redblue)

    xlims!(-3e7,5e7)


    save("figures/compare/SizeDifference_log.png", f)    
end

plot_compare_linear(samples)
plot_compare_log(samples)


