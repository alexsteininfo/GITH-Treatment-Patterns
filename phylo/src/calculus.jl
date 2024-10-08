# title: some calculations for the manuscript
# author: alexander stein

include("/Users/stein02/Desktop/Phylogenetics/phylo/src/math_analysis.jl")

# Parameters before treatment
b_1 = 1.0
d_1 = 0.0

# Parameters during treatment
b_2 = 0.25
tf = 365

# rate of neutral mutations
mu = 2.0
omega = 2*mu

Nd = 1e10
kmin = 1e8

SFS_xaxis = 1:Int(Nd)

#S0 = [predict_SFS_ExpGrowth(x,b_1,d_1,Nd,omega) for x in SFS_xaxis]
S0 = [predict_SFS_ExpGrowth_PureBirth(x,b_1,Nd,omega) for x in SFS_xaxis]

predict_SFS_ExpGrowth_PureBirth(k,b,N,omega)

Sk = 0
for kp in SFS_xaxis
    Sk += SFS_initial[kp]*pmf_CritBranch(kmin,b_2,t,kp)
end


Sk_pre = predict_SFS_ConstPre(SFS_xaxis,S0,b_2,time)
Sk_new = predict_SFS_const_new(10^4, b_2, omega, tf)

# Compute the predictions after decreasing population
for time in timesteps
    # Preexisting mutations
    #SFS_yaxis_pre = SFS_pred_pre_consul(SFS_xaxis, SFS_pred_pre[1],b_2,time)
    #SFS_yaxis_pre = SFS_pred_pre_spa(SFS_xaxis,SFS_pred_pre[1],b_2,time)
    #SFS_yaxis_pre = predict_SFS_pre_nspa(SFS_xaxis,SFS_pred_pre[1],b_2,time)
    SFS_yaxis_pre = predict_SFS_ConstPre(SFS_xaxis,SFS_pred_pre[1],b_2,time)
    push!(SFS_pred_pre, SFS_yaxis_pre)
    SFS_yaxis_new = [predict_SFS_const_new(k, b_2, omega, time) for k in SFS_xaxis]
    push!(SFS_pred_new, SFS_yaxis_new)
end

SFS_pred_tot = SFS_pred_pre .+ SFS_pred_new

S0 = []



















