
using SomaticEvolution, Random, AbstractTrees
using Distributions

# Parameters before treatment
b_s1 = 1.0
d_s1 = 0.0
b_r1 = 1.0
d_r1 = 0.0

# Parameters during treatment
b_s2 = 1.0
d_s2 = 3.0
#b_r2 = 1.0
#d_r2 = 0.0

# rate of neutral mutations
mu = 2.0
omega = 2*mu
# rate of resistance mutations
nu = 1e-2
# resistance mutation distribution
dist1 = DiscreteNonParametric([0.0],[1.0])  # death rate 0.0 with probability 1
dist2 = DiscreteNonParametric([0.0],[1.0])  # death rate 0.0 with probability 1
# Maximum number of subclones
maxclones = 1000

# other parameters
Nd = 1000

# Grow to maximal population size
input = BranchingInput(Nmax=Nd, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
selection = SelectionDistribution(dist1, nu, maxclones)
simulation = runsimulation(SimpleTreeCell, input, selection)

input2 = BranchingInput(Nmax=2.0*Nd, birthrate=b_s2, deathrate=d_s2, μ=mu, ploidy=1)
#simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input2)
selection = SelectionDistribution(dist2, nu, maxclones)
simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input2, selection)

for subclone in simulation.output.subclones
    println("Birth rate = ", subclone.birthrate)
    println("Death rate = ", subclone.deathrate)
end

#simulation.output.subclones[1].birthrate
#simulation.output.subclones[1].deathrate
#simulation.output.subclones[2].birthrate
#simulation.output.subclones[2].deathrate
#simulation.output.subclones[3].birthrate
#simulation.output.subclones[3].deathrate