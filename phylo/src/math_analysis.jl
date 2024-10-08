###########################################
### Functions for mathematical analysis ###
###########################################

# author: alexander stein
# comment: mathematical predictions for
#           (1) site frequnecy spectra
#           (2) total mutaitonal burden
#           (3) single cell mutational burden
#           and more

include("/Users/stein02/Desktop/Phylogenetics/phylo/src/pmf_BirthDeath.jl")
using SpecialFunctions: loggamma
using Distributions: Poisson, pdf

##############################
### Generic math functions ###
##############################

# Compute pmf of the Consul Poisson distribution
# Input: evalution point k, parameters alpha and beta,
#           type of parameterization meanandvariance
# Output: pmf for value k
function pmf_consul_poisson(k, alpha, beta, meanandvariance=false)
    if(meanandvariance==false)
        # alpha and beta are parameters of the Consul Poisson distributions
        mu = alpha
        lambda = beta
        # make sure we are well defined space
        if(mu<=0 || lambda<0 || lambda >= 1)
            println("ERROR: Consul Poisson is not valid with given parameters!!!")
        end
        # compute return value
        pm = log(mu)-(k*lambda+mu)+(k-1)*log(k*lambda+mu)-loggamma(k+1)
        return exp(pm)
    else
        # alpha and beta are mean and variance of the Consul Poisson distributions
        # that we transform into appropriate parameters
        mu = alpha^(3/2)/sqrt(beta)
        lambda = (sqrt(beta)-sqrt(alpha))/sqrt(beta)
        # make sure we are well defined space
        if(mu<=0 || lambda<0 || lambda >= 1)
            println("ERROR: Consul Poisson is not valid! with given parameters!!")
        end
        # and then compute as before
        pm = log(mu)-(k*lambda+mu)+(k-1)*log(k*lambda+mu)-loggamma(k+1)
        return exp(pm)
    end
end

# Compute distribution of a compound Poisson process by sampling random numbers
# Input: xaxis of interest, Poisson rates for the compound process: lam1 and lam2
#           samplesize that defines the accurancy
# Ouput: yaxis, which the (approximate) pmf of the process
function dist_compound_poisson(xaxis, lam1, lam2, samplesize=100000)
    # Define the two Poisson processes
    P1 = Poisson(lam1)
    P2 = Poisson(lam2)
    # Initilization
    dist = [0 for x in xaxis]
    # Iterate over many times and average
    for i in 1:samplesize
        # Compute the random sum
        p1 = rand(P1)
        value = sum(rand(P2,p1))
        index = findall(x->x==value, xaxis)
        dist[index] .+= 1
        #if(value-xaxis[1] >= 0 && value-xaxis[1] < xaxis[end])
        #    abundance[value-xaxis[1]] += 1
        #end
    end
    return dist/samplesize
end


##########################
### Exponential growth ###
########################## 

# Compute the predicted site frequency spectrum in the fixed time limit
#           using the time to grow to expected size N
# Input:
# Output: single value for E[S_k]
function predict_SFS_ExpGrowth(k,b,d,N,omega)
    p0=d/b
    sum = 0
    for j in 0:1000         # Compute maximum 1000 terms
        add = p0^j/((k+j)*(k+j+1))
        sum += add
        if(add<=1.0e-10)    # stop at this precision
            break
        end
    end
    return omega*N*sum
end

# Compute the predicted site frequency spectrum in the fixed time limit
#           using the time to grow to expected size N
# Input:
# Output: single value for E[S_k]
function predict_SFS_ExpGrowth_PureBirth(k,b,N,omega)
    sum = 1/(k^2+k)
    return omega*N*sum
end

# Compute the predicted total mutational burden in the fixed time limit
#           using the time to grow to size N
# Input:
# Output: single value for E[B]
function predict_tMB_ExpGrowth(b,d,N,omega)
    if(d==0)    # pure birth process
        return omega*(N-1)
    else        # using Gunnarsson's result
        rho = d/b
        sigma = 1-rho
        MB = -log(sigma-rho/N)/rho
        return omega*N*MB
    end
end

# Compute distribution of tMB after exponential growth
# Input: xaxis for distribution, birth and death rate, expected number of divisions
#           mutation rate
# Output: yaxis to xaxis of tMB distribution in form of probability density
function predict_tMB_dist(tMB_xaxis,b,d,R,m)
    tMB_pred_histo = Float64[]
    if(d_1==0)
        D = Poisson(omega*(Nd-1))
        tMB_pred_histo = [pdf(D,B) for B in tMB_xaxis]
        #tMB_pred_histo = [consul_poisson2(B,tMB_pred,tMB_pred*(1+mu)) for B in tMB_xaxis]
    else
        tMB_pred_histo = dist_compound_poisson(tMB_xaxis, tMB_pred/m, m)
        #tMB_pred_histo = [consul_poisson2(B,tMB_pred,tMB_pred) for B in tMB_xaxis]
    end
    return tMB_pred_histo
end


# Compute the distribution of a single value in the site frequency spectrum
# Input: 
# Output:
function predict_SkDist(k,Sk_simu,N,mu)
    # Define xaxis based on simulations
    min = minimum(Sk_simu)
    max = maximum(Sk_simu)
    xaxis = Int32[i for i in min:max]
    
    # Compute the predicted E[Sk]/m
    expWk = predict_SFS_ExpGrowth(k,b_1,d_1,N,2.0)
    expSk = expWk*mu
    #yaxis = Float32[pmf_consul_poisson(k, expSk, expSk*(1+mu), true) for x in xaxis]
    yaxis = dist_compound_poisson(xaxis, expWk, mu, 1000000)
    return xaxis, yaxis
end


# Compute scMB pdf for pure-birth process according to Cheek and Johnston
# Input:
# Output:
function getscMB_ExactPMF(j,b,time)
    pm = 2^j/(exp(b*time)-1)
    sum = 0 
    for i in 0:j
        logadd = i*log((b*time))-loggamma(i+1)
        sum += exp(logadd)
    end
    return pm*(1-exp(-b*time)*sum)
end


# Compute scMB pdf for pure-birth process according to our formula
# Input:
# Output:
function predict_scMBdist(scMB_xaxis, averageGen, m, N)
    # Initilization
    scMB_yaxis = Float64[]
    # Define generation distribution
    D1 = Poisson(averageGen)
    maxgen = 4*averageGen
    # Compute the sum for Mj and add to scMB_yaxis
    for j in scMB_xaxis
        Mj = 0
        for l in 1:maxgen
            D2 = Poisson(m*l)
            Mj += pdf(D1,l)*pdf(D2,j)
        end
        push!(scMB_yaxis, Mj)
    end
    return N*scMB_yaxis
end

# Predict average generation in the population in the fixed size limit
# Input
# Output
function predict_averageGeneration(b,d,N)
    # prediction for the birth-death process
    genAverage = 0
    if(d==0)
        genAverage = 2*(0.577+log(N)-1)
    else
        tN = (0.577+log(N*(b-d)/b))/(b-d)  # first time hitting size N
        genAverage = 2*b*tN - 2
    end
    return genAverage
end



#############################
### Decreasing population ###
#############################


# Predict SFS or pre-existing mutations for decreasing population given an initial SFS
# Input: SFS xaxis
# Output: SFS yaxis
function predict_SFS_DecrPre(SFS_initial, SFS_xaxis, b, d, Ni, Nf)
    # Determine time without conditioning on survival
    tf = log(Nf/Ni)/(b-d)
    # Initilization
    SFS_yaxis = Float64[0 for Sk in SFS_xaxis]
    # Compute SFS for all k in SFS_xaxis
    for k in SFS_xaxis
        Sk = 0
        for kp in SFS_xaxis
            Sk += SFS_initial[kp]*pmf_BirthDeath(k,b,d,tf,kp)
        end
        SFS_yaxis[k] = Sk
    end

    return SFS_yaxis
end


# Compute SFS numerically given the pmf for birth-death process 
#           and neglecting conditioning on survival
# Input: interested site frequency k, initial pop size Ni
# Output: 
function predict_SFS_NumInt(k, Ni, Nf, b, d, omega)
    # Determine time without conditioning on survival
    tf = log(Nf/Ni)/(b-d)
    # Solve the integral
    #integrand(t,p) = exp((b-d)*t)*pmf_BirthDeath_IntialOne(k,b,d,tf-t) pmf_BirthDeath(k,b,d,t,kp)
    integrand(t,p) = exp((b-d)*t)*pmf_BirthDeath(k,b,d,tf-t,1)
    domain = (0, tf)
    prob = IntegralProblem(integrand, domain)
    sol = solve(prob, HCubatureJL(); reltol = 1e-3)
    return omega*Ni*sol.u
end

# Compute SFS numerically given the pmf for birth-death process 
#           and neglecting conditioning on survival
# Input: interested site frequency k, initial pop size Ni
# Output: 
function predict_tMB_DecrNew(b,d,Ni,Nf,omega)
    MB = -b/d*log((b-d)/(b-d*Ni/Nf))
    return omega*Nf*MB
end

# Compute SFS numerically given the pmf for birth-death process 
#           and neglecting conditioning on survival
# Input: interested site frequency k, initial pop size Ni
# Output: 
function predict_tMB_DecrPre(SFS_initial,b,d,Ni,Nf)
    # Determine time without conditioning on survival
    tf = log(Nf/Ni)/(b-d)
    # Compute the tMB
    MB = 0
    for (k,Sk) in enumerate(SFS_initial)
        MB += Sk*(1-pmf_BirthDeath(0,b,d,tf,k))
    end
    return MB
end


###########################
### Constant population ###
###########################

# Compte the predicted SFS
# Input:
# Output:
function predict_SFS_ConstPre(SFS_xaxis,SFS_initial,b,t)
    # Initilization
    SFS_yaxis = Float64[0 for i in SFS_xaxis]
    # Sample from the initial SFS following the normalized saddle point approximation
    # Compute SFS for all k in SFS_xaxis
    for k in SFS_xaxis
        Sk = 0
        for kp in SFS_xaxis
            Sk += SFS_initial[kp]*pmf_CritBranch(k,b,t,kp)
        end
        SFS_yaxis[k] = Sk
    end

    return SFS_yaxis
end

# Compte the SFS using the saddle point approximation method for the pmf
# Input:
# Output:
function SFS_pred_pre_spa(SFS_xaxis,SFS_initial,b,t)
    SFS_yaxis = Float64[0 for i in SFS_xaxis]
    for k in SFS_xaxis
        Sk = 0
        for (K, SK_0) in enumerate(SFS_initial)
            Sk += SK_0*pmf_Critical_SPA(k,b,t,K)
        end
        SFS_yaxis[k] = Sk
    end
    return SFS_yaxis
end

# Compte the SFS using the normalized saddle point approximation method for the pmf
# Input:
# Output:
function predict_SFS_pre_nspa(SFS_xaxis,SFS_initial,b,t)
    # Initilization
    SFS_yaxis = Float64[0 for i in SFS_xaxis]
    # Sample from the initial SFS following the normalized saddle point approximation
    for (K, SK_0) in enumerate(SFS_initial)
        # We compute pmf_dist up to Nd, which leads to approximate normalization for K near Nd
        pmf_dist = dist_Critical_NSPA(SFS_xaxis,b,t,K)
        # We ignore all k > Nd
        for k in SFS_xaxis
            SFS_yaxis[k] += SK_0*pmf_dist[k]
        end
    end

    return SFS_yaxis
end

function predict_SFS_pre_consul(SFS_xaxis, SFS_initial,b,t)
    SFS_yaxis = Float64[0 for Sk in SFS_xaxis]

    for k in SFS_xaxis
        Sk = 0
        for kp in SFS_xaxis
            if(kp>=k)
                #Sk += SFS_initial[kp]*pmf_critical_saddlepoint(k,b,t,kp)
                mean = kp
                sigma2 = 2*kp*b*t
                Sk += SFS_initial[kp]*consul_poisson2(k, mean, sigma2)
            end
        end
        SFS_yaxis[k] = Sk
    end
    return SFS_yaxis
end

predict_SFS_const_new(k, b, omega, time) = omega*Nd/k*(b*time/(1+b*time))^k

function predict_tMB_pre(SFS_initial, b, t)
    MB = 0
    for (k,Sk) in enumerate(SFS_initial)
        p0k = (b*t/(1 + b*t))^k
        MB += Sk*(1-p0k)
    end
    return MB
end

predict_tMB_new(omega,N,b,t) = omega*N*log(1+b*t)


#############################
### Increasing population ###
#############################

# Compute the SFS exact for pure birth process and approximate for birth-death process
# Input:
# Outout:
function predict_SFS_IncrPre(SFS_initial, SFS_xaxis, b, d, t)
    SFS_yaxis = Float64[0 for Sk in SFS_xaxis]

    for k in SFS_xaxis
        Sk = 0
        for kp in SFS_xaxis
            #Sk += SFS_initial[kp]*pmf_BirthDeath_Consul(k,b,d,t,kp)
            Sk += SFS_initial[kp]*pmf_BirthDeath(k,b,d,t,kp)
        end
        SFS_yaxis[k] = Sk
    end
    return SFS_yaxis
end

# Compute SFS numerically given the pmf for birth-death process 
#           and neglecting conditioning on survival
# Input: interested site frequency k, initial pop size Ni
# Output: 
function predict_tMB_IncrNew(b,d,Ni,Nf,omega)
    MB = -b/d*log((b-d)/(b-d*Ni/Nf))
    return omega*Nf*MB
end

# Compute SFS numerically given the pmf for birth-death process 
#           and neglecting conditioning on survival
# Input: interested site frequency k, initial pop size Ni
# Output: 
function predict_tMB_IncrPre(SFS_initial,b,d,Ni,Nf)
    # Determine time without conditioning on survival
    tf = log(Nf/Ni)/(b-d)
    # Compute the tMB
    MB = 0
    for (k,Sk) in enumerate(SFS_initial)
        MB += Sk*(1-pmf_BirthDeath(0,b,d,tf,k))
    end
    return MB
end


###############################
### Emergence of resistance ###
###############################

# Compute CSD numerically given the pmf for birth-death process 
# Input:
# Output: 
function predict_CSD_new(k, Ni, Nf, b_s2, d_s2, b_r2, d_r2, nu, tf)
    # Solve the integral
    integrand(t,p) = exp((b_s2-d_s2)*t)*pmf_BirthDeath(k,b_r2,d_r2,tf-t,1)
    domain = (0, tf)
    prob = IntegralProblem(integrand, domain)
    sol = solve(prob, HCubatureJL(); reltol = 1e-3)
    return nu*Ni*sol.u
end

# Compute CSD 
# Input: 
# Output: 
function predict_CSD_pre(CSD_initial, Nd, Nf, b_r2, d_r2, t)
    CSD_yaxis = Float64[0 for Sk in 1:Nf]
    for k in 1:Nf
        Sk = 0
        for kp in 1:Nd
            Sk += CSD_initial[kp]*pmf_BirthDeath(k,b_r2,d_r2,t,kp)
        end
        CSD_yaxis[k] = Sk
    end
    return CSD_yaxis
end

function predict_SFS_resis(CSD, Nf, b_r, d_r, omega)
    SFS_xaxis = 1:Nf
    SFS_yaxis = Float64[]
    for k in SFS_xaxis
        Sk = 0
        for kappa in SFS_xaxis
             Sk += CSD[kappa]*predict_SFS_ExpGrowth(k,b_r,d_r,kappa,omega)
        end
        push!(SFS_yaxis, Sk)
    end
    return SFS_yaxis
end

function predict_ClonalMutations(kappa, b_s1, b_s2, b_r, d_r, td, tf)
    # auxiliary variables
    rho = d_r/b_r
    sigma = 1-rho
    # compute texpected living time
    t_kappa = log(sigma*kappa+rho)/(b_r-d_r)

    # compute the expected clonal mutations
    mut = 0
    if(t_kappa > tf+td)
        tau_d0 = td+tf
        mut = 2*b_s1*tau_d0
    elseif(t_kappa > tf)
        tau_d1 = td+tf-t_kappa
        mut = 2*b_s1*tau_d1
    else
        tau_d2 = tf-t_kappa
        mut = 2*b_s1*td + 2*b_s2*tau_d2
    end
    mut<0 ? println("ERROR in predict_ClonalMutations!!!") : nothing
    return mut
end

function predict_SFS_resis2(CSD, td, tf, b_s1, b_s2, b_r, d_r, Nf)
    SFS_xaxis = 1:Nf
    SFS_yaxis = Float64[]
    for k in SFS_xaxis
        Sk = 0
        for kappa in SFS_xaxis
             Sk = CSD[kappa]*predict_ClonalMutations(kappa, b_s1, b_s2, b_r, d_r, td, tf)*omega
        end
        push!(SFS_yaxis, Sk)
    end
    return SFS_yaxis
end



