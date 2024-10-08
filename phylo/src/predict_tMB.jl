# Comparing exact and approximate tMB

##########################
### tMB approximations ###
##########################

using CairoMakie
using Colors

include("/Users/stein02/Desktop/Phylogenetics/phylo/src/math_analysis.jl")

using ColorBrewer
ColPaired = palette("Paired", 12)

##################
### parameters ###
##################

### Exponential growth
b1 = 1.0
d1 = 0.5

mu = 2.0
omega = 2*mu

Nd = Int(1e5)

tMB_detect = predict_tMB_ExpGrowth(b1,d1,Nd,omega)

### Decreasing population
b2_decr = 1.0
d2_decr = 3.0

### Constant population
b2_const = 1.0
d2_const = 1.0

### Increasing population
b2_incr = 1.0
d2_incr = 0.5

#######################
### Approximate tMB ###
#######################

function predict_tMB_new(b2,d2,Nd,tf,m)
    L = exp((b2-d2)*tf)
    alpha = (d2*L-d2)/(b2*L-d2)
    
    #tf = ln( (1-alpha)*Nf/Nd ) / (b2-d2)
    Nf = 1/(1-alpha)*Nd*exp((b2-d2)*tf)
    
    tMB = 2*m*Nf*(b2/d2)*log( b2/(b2-d2) - 1/(1-alpha)/Nd*d2/(b2-d2) )
    
    return tMB
end

function predict_tMB_pre(b1,d1,b2,d2,Nf,tf,m)
    L = exp((b2-d2)*tf)
    alpha = (d2*L-d2)/(b2*L-d2)
    return tMB_detect*(alpha-1)/alpha*log((1-alpha))
    #return 2*m*Nd*b1/(b1-d1)*log((1-alpha))
end

function predict_tMB_pre(b1,d1,b2,Nf,tf,m)
    alpha = (b2*tf)/(1+b2*tf)
    return tMB_detect*(alpha-1)/alpha*log((1-alpha))
    #return 2*m*Nd*b1/(b1-d1)*log((1-alpha))
end


function predict_tMB_pre2(b1,d1,b2,d2,Nf,tf,m)
    L = exp((b2-d2)*tf)
    alpha = (d2*L-d2)/(b2*L-d2)

    part1 = (alpha-1)*Nd-alpha
    part2 = log(alpha*(1/Nd-1) + 1)
    part3 = - alpha*log(Nd)
    part4 = alpha*(Nd-1)

    #sum = (part1*part2+part3)/part4

    sum = ((-alpha + (-1 + alpha)*Nd)*log(1 + alpha*(-1 + 1/Nd)) - alpha*log(Nd))/(alpha*(-1 + Nd))
    
    return 2*m*Nd*sum
end

function predict_tMB_pre2(b1,d1,b2,Nf,tf,m)
    alpha = (b2*tf)/(1+b2*tf)
    
    part1 = (alpha-1)*Nd-alpha
    part2 = log(alpha*(1/Nd-1) + 1)
    part3 = - alpha*log(Nd)
    part4 = alpha*(Nd-1)

    #sum = (part1*part2+part3)/part4

    sum = ((-alpha + (-1 + alpha)*Nd)*log(1 + alpha*(-1 + 1/Nd)) - alpha*log(Nd))/(alpha*(-1 + Nd))
    
    return 2*m*Nd*sum
end

#((-alpha + (-1 + alpha) Nd) Log[1 + alpha (-1 + 1/Nd)] - alpha Log[Nd])/(alpha (-1 + Nd))

#############################
### Decreasing population ###
#############################

### tMB at detection
tMB_det = predict_tMB_ExpGrowth(b1, d1, Nd, omega)

### SFS at detection
SFS_xaxis = 1:Nd
SFS_det = [predict_SFS_ExpGrowth(x,b1,d1,Nd,omega) for x in SFS_xaxis]

### Decreasing population
Nf_decr0 = Int(0.5*Nd)
tf_decr0 = log(Nf_decr0/Nd)/(b2_decr-d2_decr)

println("tMB_pre_decr_exact = ", predict_tMB_DecrPre(SFS_det, b2_decr, d2_decr, Nd, Nf_decr0))
println("tMB_pre_decr_approx = ", predict_tMB_pre(b1, d1, b2_decr, d2_decr, Nf_decr0, tf_decr0, mu))
println("tMB_pre_decr_approx2 = ", predict_tMB_pre2(b1, d1, b2_decr, d2_decr, Nf_decr0, tf_decr0, mu))


size_decr = [Nd*(1-0.5*i/500) for i in 1:500]
time_decr = [log(Nf/Nd)/(b2_decr-d2_decr) for Nf in size_decr]

tMB_pre_decr_exact = [predict_tMB_DecrPre(SFS_det, b2_decr, d2_decr, Nd, Nf_decr) for (Nf_decr, tf_decr) in zip(size_decr, time_decr)]
tMB_pre_decr_approx = [predict_tMB_pre(b1, d1, b2_decr, d2_decr, Nf_decr, tf_decr, mu) for (Nf_decr, tf_decr) in zip(size_decr, time_decr)]
tMB_pre_decr_approx2 = [predict_tMB_pre2(b1, d1, b2_decr, d2_decr, Nf_decr, tf_decr, mu) for (Nf_decr, tf_decr) in zip(size_decr, time_decr)]


### Constant population
Nf_const0 = Nd
tf_const0 = 2.0

println("tMB_pre_const_exact = ", predict_tMB_pre(SFS_det, b2_const, tf_const0))
println("tMB_pre_const_approx = ", predict_tMB_pre(b1, d1, b2_const, Nf_const0, tf_const0, mu))
println("tMB_pre_const_approx2 = ", predict_tMB_pre2(b1, d1, b2_const, Nf_const0, tf_const0, mu))

time_const = [i*0.01 for i in 1:500]

tMB_pre_const_exact = [predict_tMB_pre(SFS_det, b2_const, tf_const) for tf_const in time_const]
tMB_pre_const_approx = [predict_tMB_pre(b1, d1, b2_const, Nf_const0, tf_const, mu) for tf_const in time_const]
tMB_pre_const_approx2 = [predict_tMB_pre2(b1, d1, b2_const, Nf_const0, tf_const, mu) for tf_const in time_const]

### Increasing population
Nf_incr0 = Int(2.0*Nd)
tf_incr0 = log(Nf_incr0/Nd)/(b2_incr-d2_incr)

println("tMB_pre_incr_exact = ", predict_tMB_IncrPre(SFS_det, b2_incr, d2_incr, Nd, Nf_incr0))
println("tMB_pre_incr_approx = ", predict_tMB_pre(b1, d1, b2_incr, d2_incr, Nf_incr0, tf_incr0, mu))
println("tMB_pre_incr_approx2 = ", predict_tMB_pre2(b1, d1, b2_incr, d2_incr, Nf_incr0, tf_incr0, mu))


size_incr = [Nd*(1+1.0*i/500) for i in 1:500]
time_incr = [log(Nf/Nd)/(b2_incr-d2_incr) for Nf in size_incr]

#predict_tMB_IncrPre(SFS_pred_pre[1], b_2, d_2,Nd,Nmax)

tMB_pre_incr_exact = [predict_tMB_IncrPre(SFS_det, b2_incr, d2_incr, Nd, Nf_incr) for (Nf_incr, tf_incr) in zip(size_incr, time_incr)]
tMB_pre_incr_approx = [predict_tMB_pre(b1, d1, b2_incr, d2_incr, Nf_incr, tf_incr, mu) for (Nf_incr, tf_incr) in zip(size_incr, time_incr)]
tMB_pre_incr_approx2 = [predict_tMB_pre2(b1, d1, b2_incr, d2_incr, Nf_incr, tf_incr, mu) for (Nf_incr, tf_incr) in zip(size_incr, time_incr)]


predict_tMB_IncrPre(SFS_det, b2_incr, d2_incr, Nd, size_incr[end])
predict_tMB_pre(b1, d1, b2_incr, d2_incr, size_incr[end], time_incr[end], mu)

##########################
### Plotting functions ###
##########################

function plot_compare(timeaxis, tMB_exact, tMB_approx, tMB_approx2, label)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = rich("Treatment time, 𝑡", subscript("f", font = :italic),font = :regular), xlabelsize = 30, xticklabelsize = 25,
    ylabel = rich("tMB, E[𝘉", superscript("(pre)"), "]"), ylabelsize = 30, yticklabelsize = 25
    )

    lines!(timeaxis, tMB_exact, color=ColPaired[8], label = "exact", linewidth=7.0, linestyle = (:dot,1.5))
    lines!(timeaxis, tMB_approx, color=ColPaired[7], label = "approx.", linewidth=7.0)
    #lines!(timeaxis, tMB_approx2, color=:black, label = "approx2.", linewidth=4.0)

    #axislegend(ax, labelsize=25, framecolor = :white)
    #Label(f[1, 1, TopLeft()], "E", fontsize = 36, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save(label, f)    
end

plot_compare(time_decr, tMB_pre_decr_exact, tMB_pre_decr_approx, tMB_pre_decr_approx2, "figures/compare/decrease.png")
plot_compare(time_const, tMB_pre_const_exact, tMB_pre_const_approx, tMB_pre_decr_approx2, "figures/compare/constant.png")
plot_compare(time_incr, tMB_pre_incr_exact, tMB_pre_incr_approx, tMB_pre_decr_approx2, "figures/compare/increase.png")
