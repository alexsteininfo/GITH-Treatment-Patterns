#############################################################
### Probability mass functions of the birth-death process ###
#############################################################

# author: alexander stein
# comment: exact and approximate values for the probability mass functions
#           of the birth-death process p_n(t) given b, d and a


using SpecialFunctions: loggamma
using Distributions: Binomial, NegativeBinomial, pdf

##############################
### All functions combined ###
##############################

# A combined function for all functions
# Input: birth rate b, death rate d, time t, initial population size a,
#           and final population size n
# Output: probability to have size n after time t
function pmf_BirthDeath(n,b,d,t,a)
    ### Check for valid input
    if(b<0 || d<0 || t<0 || a<1 || n<0)
        println("ERROR: Invalid Input!!!")
    end
    
    ### Select the right formula based on parameters 
    if(b!=d)
        if(d==0)        # pure birth process
            return pmf_PureBirth(n,b,t,a)
        elseif(b==0)    # pure death process
            return pmf_PureDeath(n,d,t,a)
        elseif(a==1)    # general b and d but starting from 1 cell
            return pmf_BirthDeath_InitialOne(n,b,d,t)
        elseif(n<100 || a<100)    # most general solution, exact
            return pmf_BirthDeath_Tavare(n,b,d,t,a)
        else                    # most general solution, SPA
            return pmf_BirthDeath_SPA(n,b,d,t,a)
        end
    else
        println("Warning: b=d, we recommend using pmf_CritBranch!!!")
        return pmf_CritBranch(n,b,t,a)
    end
end

function pmf_CritBranch(n,b,t,a)
    ### Check for valid input
    if(b<0 || t<0 || a<1 || n<0)
        println("ERROR: Invalid Input!!!")
        return 0    
    ### Select the right formula based on parameters 
    elseif(b>0)         # non-trivial process
        if(a==1)    # starting from 1 cell
            return pmf_CritBranch_InitialOne(n,b,t)
        elseif(n<50 || a<50)    # most general solution, exact
            return pmf_CritBranch_Tavare(n,b,t,a)
        else                    # most general solution, SPA
            return pmf_CritBranch_SPA(n,b,t,a)
        end
    else            # trivial process
        println("Warning: b=0 and d=0 is a trivial process!!!")
        if(a==n)
            return 1
        else
            return 0
        end
    end
end


###############################
### Exact formulas for pmfs ###
###############################

# Compute general pmf for pure-birth process as given in eqn. (8.15) in 
#           "The Elements of Stochstic Processes" by Norman Bailey
# Input: birth rate b, time t and initial population size a
# Output: single value for the probability mass to have size n
function pmf_PureBirth(n,b,t,a)
    successprob = exp(-b*t)
    P = NegativeBinomial(a,successprob)
    return pdf(P,a+n)
end


# Compute general pmf for pure-death process as given in eqn. (8.30) in 
#           "The Elements of Stochstic Processes" by Norman Bailey
# Input: death rate d, time t and initial population size a
# Output: single value for the probability mass to have size n
function pmf_PureDeath(n,d,t,a)
    P = Binomial(a,exp(-d*t))
    return pdf(P,n)
end


# Compute pmf for birth-death process starting from 1 cell as given in eqn. (8.46)
#           in "The Elements of Stochstic Processes" by Norman Bailey
# Input: variables for the birth death process: birth rate b, death rate d,
#           time t
# Output: single value for the probability mass
function pmf_BirthDeath_InitialOne(n,b,d,t)
    # Define auxiliary variables
    L = exp((b-d)*t)
    alpha = (d*L-d)/(b*L-d)
    beta = (b*L-b)/(b*L-d)
    # return the value according to the formula
    if(n==0)
        return alpha
    else
        return (1-alpha)*(1-beta)*beta^(n-1)
    end
end


# Compute general pmf for birth-death process as given in eqn. (8.47) in 
#           "The Elements of Stochstic Processes" by Norman Bailey
# Input: variables for the birth death process: birth rate b, death rate d,
#           time t, and initial population size a
# Output: single value for the probability mass
function pmf_BirthDeath_Bailey(n,b,d,t,a)
    # Define auxiliary variables
    L = exp((b-d)*t)
    alpha = (d*L-d)/(b*L-d)
    beta = (b*L-b)/(b*L-d)
    gamma = abs(1-alpha-beta)
    # Intilize probability mass
    pm = 0
    # Sum of Bailey using loggamma function to avoid machine errors
    if(n==0)
        return alpha^a
    elseif(a+n<50)  # We can compute binomial coefficients exactly
        for j in 0:minimum([a,n])
            part1 = binomial(a,j)
            part2 = binomial(a+n-j-1,a-1)
            part3 = alpha^(a-j)
            part4 = beta^(n-j)
            part5 = (1-alpha-beta)^j
            pm += part1*part3*part4*part5*part2
        end
        return pm
    else        # We use the x = exp(log(x)) trick and use loggamma function
        if(1-alpha-beta<0)
            for j in 0:minimum([a,n])
                part1 = loggamma(a+1)-loggamma(j+1)-loggamma(a-j+1)
                part2 = loggamma(a+n-j)-loggamma(a)-loggamma(n-j+1)
                part3 = (a-j)*log(alpha)
                part4 = (n-j)*log(beta)
                part5 = j*log(gamma)
                pm += (-1)^j*exp(part1+part2+part3+part4+part5)
            end
        elseif(1-alpha-beta>0)
            for j in 0:minimum([a,n])
                part1 = loggamma(a+1)-loggamma(j+1)-loggamma(a-j+1)
                part2 = loggamma(a+n-j)-loggamma(a)-loggamma(n-j+1)
                part3 = (a-j)*log(alpha)
                part4 = (n-j)*log(beta)
                part5 = j*log(gamma)
                pm += exp(part1+part2+part3+part4+part5)
            end
        end
        return pm
    end
end


# Compute general pmf for birth-death process as given in eqn. (11) in 
#           "The linear birth-death process: an inferential retrospective" 
#           by Simon Tavare
# Input: variables for the birth death process: birth rate b, death rate d,
#           time t, and initial population size a
# Output: single value for the probability mass
function pmf_BirthDeath_Tavare(n,b,d,t,a)
    ### Define auxiliary variables
    L = exp((b-d)*t)
    alpha = (d*L-d)/(b*L-d)
    beta = (b*L-b)/(b*L-d)
    
    ### Compute pstar
    pstar = Float64[]
    # For j=0
    pstar0 = 0.0
    if(n==0)
        pstar0 = 1.0
    else
        pstar0 = 0.0
    end
    # For j>0
    for j in 1:minimum([n,a])
        D = NegativeBinomial(j, 1-beta)
        push!(pstar, pdf(D,n-j))
    end
    
    ### Compute p
    # First, for j=0
    p = alpha^a*pstar0
    # Then for j>0
    for j in 1:minimum([n,a])
        #part1 = log(binomial(a,j))
        part1 = loggamma(a+1)-loggamma(j+1)-loggamma(a-j+1)
        part2 = j*log(1-alpha) + (a-j)*log(alpha)
        p += exp(part1+part2+log(pstar[j]))
    end
    return p
end


# Compute pmf for critical birth-death process starting from 1 cell as given
#           in eqn. (8.53) in "The Elements of Stochstic Processes" by Norman Bailey
# Input: birth rate b (that equals the death rate) and time t
# Output: single value for the probability mass to have size n
function pmf_CritBranch_InitialOne(n,b,t)
    if(n==0)
        return b*t/(1+b*t)
    else
        part1 = (n-1)*log(b*t)
        part2 = (n+1)*log(1+b*t)
        return exp(part1-part2)
        #return (b*t)^(n-1)/(1+b*t)^(n+1)
    end
end


# Compute pmf for critical birth-death process for n < 4 computed
#           using Taylor expansion from eqn. (8.52) in 
#           "The Elements of Stochstic Processes" by Norman Bailey
# Input: birth rate b (that equals the death rate) and time t
# Output: single value for the probability mass to have size n
function pmf_CritBranch_MyAttempt(n,b,t,a)
    # Implement function for n=0,1,2,3,4
    if(n==0)
        return (b*t/(1 + b*t))^a
    elseif(n==1)
        part1 = (b*t/(1 + b*t))^a
        part2 = a/(b*t+1)
        part3 = 1/(b*t)
        return part1*part2*part3
    elseif(n==2)
        part1 = (b*t/(1 + b*t))^a
        part2 = a/(b*t+1)^2
        part3 = (a+2*(b*t)^2-1)/(2*(b*t)^2)
        return part1*part2*part3
    elseif(n==3)
        part1 = (b*t/(1 + b*t))^a
        part2 = a/(b*t+1)^3
        part3 = (a^2+6*(b*t)^2-3*a+6*(b*t)^4-6*(b*t)^2+2)/(6*(b*t)^3)
        println("Likely contains a typo in pmf for n=3!!!")
        return part1*part2*part3
    elseif(n==4)
        part1 = (b*t/(1 + b*t))^a
        part2 = a/(b*t+1)^4
        part3 = (a^3+12*a^2*(b*t)^2-6*a^2+36*a*(b*t)^4-36*a*(b*t)^2+11*a+24*(b*t)^6-36*(b*t)^4+24*(b*t)^2-6)/(24*(b*t)^4)
        return part1*part2*part3
    else 
        println("ERROR: pmf_CritBranch not valid for given parameters!!!")
        return 0
    end
end

# Compute pmf for critical birth-death process as given in 
#           "The Elements of Stochstic Processes" by Norman Bailey
# Input: birth rate b (that equals the death rate) and time t
# Output: single value for the probability mass to have size n
function pmf_CritBranch_Bailey(n,b,t,a)
    # Define auxiliary variables
    alpha = b*t/(1+b*t)
    beta = alpha

    # Intilize probability mass
    pm = 0
    # Sum of Bailey using loggamma function to avoid machine errors
    if(n==0)
        return alpha^a
    elseif(a+n<50)  # We can compute binomial coefficients exactly
        for j in 0:minimum([a,n])
            part1 = binomial(a,j)
            part2 = binomial(a+n-j-1,a-1)
            part3 = alpha^(a-j)
            part4 = beta^(n-j)
            part5 = (1-alpha-beta)^j
            pm += part1*part3*part4*part5*part2
        end
        return pm
    else        # We use the x = exp(log(x)) trick and use loggamma function
        if(1-alpha-beta<0)
            for j in 0:minimum([a,n])
                part1 = loggamma(a+1)-loggamma(j+1)-loggamma(a-j+1)
                part2 = loggamma(a+n-j)-loggamma(a)-loggamma(n-j+1)
                part3 = (a-j)*log(alpha)
                part4 = (n-j)*log(beta)
                part5 = j*log(gamma)
                pm += (-1)^j*exp(part1+part2+part3+part4+part5)
            end
        elseif(1-alpha-beta>0)
            for j in 0:minimum([a,n])
                part1 = loggamma(a+1)-loggamma(j+1)-loggamma(a-j+1)
                part2 = loggamma(a+n-j)-loggamma(a)-loggamma(n-j+1)
                part3 = (a-j)*log(alpha)
                part4 = (n-j)*log(beta)
                part5 = j*log(gamma)
                pm += exp(part1+part2+part3+part4+part5)
            end
        end
        return pm
    end
end


# Compute general pmf for birth-death process as given in eqn. (11) in 
#           "The linear birth-death process: an inferential retrospective" 
#           by Simon Tavare
# Input: variables for the birth death process: birth rate b, death rate d,
#           time t, and initial population size a
# Output: single value for the probability mass
function pmf_CritBranch_Tavare(n,b,t,a)
    ### Define auxiliary variables
    alpha = b*t/(1+b*t)
    beta = alpha
    
    ### Compute pstar
    pstar = Float64[]
    # For j=0
    pstar0 = 0.0
    if(n==0)
        pstar0 = 1.0
    else
        pstar0 = 0.0
    end
    # For j>0
    for j in 1:minimum([n,a])
        D = NegativeBinomial(j, 1-beta)
        push!(pstar, pdf(D,n-j))
    end
    
    ### Compute p
    # First, for j=0
    p = alpha^a*pstar0
    # Then for j>0
    for j in 1:minimum([n,a])
        #part1 = log(binomial(a,j))
        part1 = loggamma(a+1)-loggamma(j+1)-loggamma(a-j+1)
        part2 = j*log(1-alpha) + (a-j)*log(alpha)
        p += exp(part1+part2+log(pstar[j]))
    end

    return p
end


######################################
### Exact expectation and variance ###
######################################

Exp_BirthDeath(b,d,t,a) = a*exp((b-d)*t)
Var_BirthDeath(b,d,t,a) = a*(b+d)/(b-d)*exp((b-d)*t)*(exp((b-d)*t)-1)

Exp_CritBranch(b,t,a) = a
Var_CritBranch(b,t,a) = 2*a*b*t


###############################
### Approximations for pmfs ###
###############################

# Compute pmf for critical birth-death process using the the saddlepoint
#           approximation method based on the MGF in eqn. (8.51) in 
#           "The Elements of Stochstic Processes" by Norman Bailey
# Input: birth rate b (that equals the death rate) and time t, and
#           initial population size a
# Output: single value for the probability mass to have size n
#=
function pmf_CritBranch_SPA(n,b,t,a)
    theta = log( (a+2*(b*t)^2*n-n-sqrt((a-n)^2+4*a*n*(b*t)^2)) / (2*b*t*(b*t-1)*n) )
    CGF = a*log( (1-(b*t-1)*(exp(theta)-1)) / (1-b*t*(exp(theta)-1)) )
    ddCGF = 2*a*b*t*exp(2*theta)*(cosh(theta)-b*t*sinh(theta)) / (b*t*(exp(theta)-1)-1)^2 / (exp(theta)*(1-b*t) + b*t)^2
    pm = exp(CGF-theta*n)/sqrt(2*pi*ddCGF)
    return pm
end
=#

# Compute pmf for critical birth-death process using the the saddlepoint
#           approximation method based on the MGF described in Davison et al.
#           "Parameter estimation for discretely observed linear birth-and-death processes"
# Input: birth rate b (that equals the death rate) and time t, and
#           initial population size a
# Output: single value for the probability mass to have size n
function pmf_CritBranch_SPA(n,b,t,a)
    # Compute the log saddlepoint
    A=b*t-(b*t)^2
    B=2*(b*t)^2+(a/n)-1
    C=-b*t-(b*t)^2
    s=(1/(2*A))*(-B+sqrt(B^2-4*A*C))
    # Return the pmf
    pmf = (1/(sqrt(2*pi*a)))*(1/s^n)*((b*t*(1-s)+s)/(1-b*t*(s-1)))^a*((b*t*s*(-b*t*s^2+b*t+s^2+1))/((b*t*(s-1)-1)^2*(-b*t*s+b*t+s)^2))^(-1/2)
    isnan(pmf) ? println("ERROR: Saddlepoint approximation failed. Division by zero may be the cause!!!") : nothing
    return pmf
end


# Compute pmf for birth-death process using the the saddlepoint
#           approximation method based on the MGF described in Davison et al.
#           "Parameter estimation for discretely observed linear birth-and-death processes"
# Input: birth rate b (that equals the death rate) and time t, and
#           initial population size a
# Output: single value for the probability mass to have size n
function pmf_BirthDeath_SPA(n,b,d,t,a)
    # auxiliary variables
    L1 = exp((b-d)*t)
    # Compute the log saddlepoint
    A=b*(L1-1)*(b-d*L1)
    B=2*b*d*(1+L1^2-L1-(a/n)*L1)+L1*(b^2+d^2)*((a/n)-1)
    C=d*(L1-1)*(d-b*L1)
    s=(1/(2*A))*(-B+sqrt(B^2-4*A*C))
    # Compute pmf
    #pmf = (1/(sqrt(2*pi*a)))*1/s^n*((d-b*s+d*(s-1)*L1)/(d-b*s+b*(s-1)*L1))^a*((-1)*( (L1-1)*L1*s*(b-d)^2*(-b^2*s^2+b*L1*d*(s^2-1)+d^2))/((b*(L1*(s-1)-s)+d)^2*(b*s+d*(-L1*s+L1-1))^2))^(-1/2)
    part1 = -log(sqrt(2*pi*a))
    part2 = -n*log(s)
    part3 = a*log((d-b*s+d*(s-1)*L1)/(d-b*s+b*(s-1)*L1))
    part4 = log(((-1)*( (L1-1)*L1*s*(b-d)^2*(-b^2*s^2+b*L1*d*(s^2-1)+d^2))/((b*(L1*(s-1)-s)+d)^2*(b*s+d*(-L1*s+L1-1))^2))^(-1/2))
    pmf = exp(part1+part2+part3+part4)
    isnan(pmf) ? println("ERROR: Saddlepoint approximation failed. Division by zero may be the cause!!!") : nothing
    return pmf
end


# Compute pmf for critical birth-death process using the the normalized
#           saddlepoint approximation method based on the MGF in eqn. (8.51) 
#           in "The Elements of Stochstic Processes" by Norman Bailey
# Input: birth rate b (that equals the death rate) and time t, and
#           initial population size a, and normalization boundary bound
# Output: single value for the probability mass to have size n
function pmf_CritBranch_NSPA(n,b,t,a,boundary)
    # Compute normalization constant
    min, max = boundary
    yaxis = [pmf_CritBranch_SPA(np,b,t,a) for np in min:max]
    C = sum(yaxis)
    # Compute the zero value
    p0_crit = (b*t/(1 + b*t))^a
    # Compute normalized saddlepoint approximation pmf
    pm = (1-p0_crit)/C*pmf_CritBranch_SPA(n,b,t,a)
    return pm
end


# Approximate pmf for critical birth-death process using the the normalized
#           saddlepoint approximation method based on the MGF in eqn. (8.51) 
#           in "The Elements of Stochstic Processes" by Norman Bailey
# Input: birth rate b (that equals the death rate) and time t, and
#           initial population size a
# Output: pmf distribution over given xaxis normalized over that xaxis
function dist_CritBranch_NSPA(xaxis,b,t,a)
    # Compute un-normalized saddlepoint approximation 
    yaxis = [pmf_CritBranch_SPA(np,b,t,a) for np in xaxis]
    # Compute normalization constant
    C_norm = sum(yaxis)
    # Compute the zero value
    p0_crit = (b*t/(1 + b*t))^a
    # Compute normalized saddlepoint approximation 
    yaxis = (1-p0_crit)/C_norm*yaxis
    return yaxis
end


# Approximate pmf with the Consul Poisson distribution for 
#           sub and supercritical branching process
# Input: birth rate b, death rate d, time t, initial size a
# Output: pmf for value n
function pmf_BirthDeath_Consul(n,b,d,t,a)
    # Compute mean and variance
    mean = Exp_BirthDeath(b,d,t,a)
    variance = Var_BirthDeath(b,d,t,a)
    # Transform mean and variance into appropriate parameters
    mu = mean^(3/2)/sqrt(variance)
    lambda = (sqrt(variance)-sqrt(mean))/sqrt(variance)
    # Make sure we are well defined space
    if(mu<=0 || lambda<0 || lambda >= 1)
        println("ERROR: Consul Poisson is not valid with given parameters!!")
    end
    # Compute probability mass
    pm = log(mu)-(n*lambda+mu)+(n-1)*log(n*lambda+mu)-loggamma(n+1)
    return exp(pm)
end


# Approximate pmf with the Consul Poisson distribution
#           for critical branching process
# Input: birth rate b, death rate d, time t, initial size a
# Output: pmf for value n
function pmf_CritBranch_Consul(n,b,t,a)
    # Compute mean and variance
    mean = Exp_CritBranch(b,t,a)
    variance = Var_CritBranch(b,t,a)
    # Transform mean and variance into appropriate parameters
    mu = mean^(3/2)/sqrt(variance)
    lambda = (sqrt(variance)-sqrt(mean))/sqrt(variance)
    # Make sure we are well defined space
    if(mu<=0 || lambda<0 || lambda >= 1)
        println("ERROR: Consul Poisson is not valid with given parameters!!")
    end
    # Compute probability mass
    pm = log(mu)-(n*lambda+mu)+(n-1)*log(n*lambda+mu)-loggamma(n+1)
    return exp(pm)
end

#=
function pmf_BirthDeath_SPA_Davison(n,b,d,t,a)
    # Define auxiliary variables
    L = a*exp((b-d)*t)
    A = b*(L-1)*(b-d*L)
    B = 2*b*d*(1+L^2-L-a/n*L) + L*(b^2+d^2)*(a/n-1)
    C = d*(L-1)*(d-b*L)
    # Compute exp saddlepoint s = exp(x)
    s = (-B+sqrt(B^2-4*A*C))/(2*A)
    # Compute pmf
    part1 = 1/sqrt(2*pi*a)*1/s^n
    part2 = ( ( d-b*s+d*(s-1)*L )/( d-b*s+b*(s-1)*L ) )^a
    part3 = -(L-1)*L*s*(b-d)^2*(-b^2*s^2+b*L*d*(s^2-1)+d^2)
    part4 = b*(L*(s-1)-s)+d
    part5 = b*s+d*(-L*s+L-1)
    return part1*part2*(part3/part4^2/part5^2)^(-1/2)
end
=#


##################
### Test space ###
##################
#=
pmf_CritBranch_Tavare(1,1.0,1.0,3)
pmf_CritBranch_MyAttempt(1,1.0,1.0,3)
pmf_CritBranch_Bailey(1,1.0,1.0,3)


pmf_BirthDeath_Tavare(n,2.0,1.0,1.0,10)
pmf_BirthDeath_Bailey(n,2.0,1.0,1.0,10)

println(pmf_BirthDeath_SPA_Davison(50,7.0,6.0,1.0,20))

sum1 = 0
sum2 = 0
for n in 1:30
    println(pmf_BirthDeath_Tavare(n,7.0,5.0,1.0,10))
    #println(pmf_BirthDeath_Bailey(n,7.0,5.0,1.0,10))
    println(pmf_BirthDeath_SPA_Davison(n,7.0,5.0,1.0,10))
    println(pmf_BirthDeath_SPA_Davison3(n,7.0,5.0,1.0,10))
    println("and")
end

pmf_BirthDeath_SPA_Davison(10,7.0,6.0,1.0,10)

sum = 0
for n in 1:100
    #println(pmf_CritBranch_MyAttempt(n,1.0,2.0,50))
    #println(pmf_CritBranch_Bailey(n,1.0,2.0,50))
    println(pmf_CritBranch_Tavare(n,1.0,2.0,50))
    #sum += pmf_CritBranch_Tavare(n,1.0,2.0,50)
    println(pmf_CritBranch_SPA_Davison(n,1.0,2.0,50))
    println(pmf_CritBranch_SPA(n,1.0,2.0,50))
    println("and")
end
=#

