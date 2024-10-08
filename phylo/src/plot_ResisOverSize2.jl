# title: Plotting resistant subclone size over time
# author: alexander stein

using CairoMakie
using Distributions
using LandauDistribution

########################
### Parameterization ###
########################

nu = 1e-7
xaxis_popsizes = exp10.(range(start=8, stop=11, length=100))

b=0.14
d=0.13
sigma = (b-d)/b
sigma_inv = b/(b-d)
#sigma = 1

# Parameters of Kessler and levine as test for cdf
nu_test = 1e-7
N_test = Int(1e8)

###################################
### Probabilities and quantiles ###
###################################

numres(nu,N) = 2*nu*N*sigma_inv
numres0(nu,N) = 2*nu*N

# Define functions for the biggest clone size
kappa_pdf(kappa, nu, N) = numres(nu,N)/kappa^2*(1-1/kappa)^(numres(nu,N)-1)
kappa_cdf(kappa, nu, N) = (1 - 1/kappa)^(numres(nu,N))

freq_pdf(freq, nu, N) = kappa_pdf(freq*N, nu, N)*N  #numres/(freq*N)^2*(1-1/(freq*N)))^(nu*N-1)*N
freq_cdf(freq, nu, N) = kappa_cdf(freq*N, nu, N)    #(1 - 1/(freq*N))^(nu*N)

cumprobsize(percentage, nu, N) = -(1/((-1 + percentage^(1/numres(nu,N)))))  # returns median size
cumprobfreq(percentage, nu, N) = -(1/((-1 + percentage^(1/numres(nu,N)))*N))  # returns median frequency
#cumprobsize(0.5, nu_test, N_test)
#cumprobfreq(0.5, nu_test, N_test)

# Define functions for total number of resistant cells
Lan = Landau(1.0)
ResPop_pdf(m, nu, N) = 1/numres(nu,N)*pdf(Lan, m/numres(nu,N) - log(numres0(nu,N)) )
ResPop_cdf(x, nu, N) = cdf(Lan, x/numres(nu,N) - log(numres0(nu,N)) )
#ResPop_cdf(20, nu_test, N_test)

ResFreq_pdf(freq, nu, N) = ResPop_pdf(freq*N, nu, N)*N
ResFreq_cdf(freq, nu, N) = ResPop_cdf(freq*N, nu, N)

# Compute the 10%, 50% and 90% quantiles
function quantile_freq(percentage, nu, N)
    result = 0
    for freq in range(start=0.0, stop=1.0, length=100000000)
        CumFunc = ResFreq_cdf(freq, nu, N)
        if(CumFunc>percentage)
            result = freq
            break
        end
    end
    return result
end

#quantile_freq(0.5, nu_test, N_test)

function quantile_size(percentage, nu, N)
    result = 0
    for freq in range(start=0, stop=N, length=100000000)
        CumFunc = ResPop_cdf(freq, nu, N)
        if(CumFunc>percentage)
            result = freq
            break
        end
    end
    return result
end

#quantile_size(0.5, nu_test, N_test)

########################
### Prepare the plot ###
########################

function plot_quantilesoversize(xaxis, yaxis1, yaxis2)
    f = Figure()
    ax = Axis(f[1,1],
    xlabel = rich("Population size at detection, 𝘕", subscript("d")), xlabelsize = 30, xticklabelsize = 25,
    ylabel = "Frequency", ylabelsize = 30, yticklabelsize = 25,
    xticks = [3e10, 6e10, 9e10]
    )

    lines!(ax, xaxis, yaxis2[2], color=:blue, linewidth=4.0, label=L"f_R = R/N_d")
    band!(ax, xaxis, yaxis2[1], yaxis2[3], color = (:blue, 0.2))

    lines!(ax, xaxis, yaxis1[2], color=:green, linewidth=4.0, label=L"f_\kappa = \kappa/N_d")
    band!(ax, xaxis, yaxis1[1], yaxis1[3], color = (:green, 0.2))

    xlims!(ax, 0.95*minimum(xaxis), 1.0*maximum(xaxis))
    ylims!(ax, 0.0, 1.3*maximum(yaxis2[3]))
    axislegend(ax, labelsize=25, framecolor = :white, position = :lt)

    save("figures/Resis/quantilesoversize.png", f)
end

kappa_01 = [cumprobfreq(0.1, nu, x) for x in xaxis_popsizes]
kappa_05 = [cumprobfreq(0.5, nu, x) for x in xaxis_popsizes]
kappa_90 = [cumprobfreq(0.9, nu, x) for x in xaxis_popsizes]

kappa_all = [kappa_01, kappa_05, kappa_90]

m_01 = [quantile_freq(0.1, nu, x) for x in xaxis_popsizes]
m_05 = [quantile_freq(0.5, nu, x) for x in xaxis_popsizes]
m_90 = [quantile_freq(0.9, nu, x) for x in xaxis_popsizes]

m_all = [m_01, m_05, m_90]


plot_quantilesoversize(xaxis_popsizes, kappa_all, m_all)









