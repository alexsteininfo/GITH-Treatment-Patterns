### Testing Jessie's new simulation

using SomaticEvolution

b_1 = 1.0
d_1 = 0.0

b_2 = 1.0
d_2 = 2.0

mu = 2.0
omega = 2*mu

Nd =  1000

input11 = BranchingInput(Nmax=Nd, birthrate=b_1, deathrate=d_1, μ=mu, ploidy=1)
input22 = MoranInput(N=Nd, tmax=20.0, moranrate=b_2, μ=mu, ploidy=1)
input33 = BranchingInput(Nmax=10000, birthrate=b_1, deathrate=d_1, μ=mu, ploidy=1)
input44 = BranchingInput(Nmax=100, birthrate=b_2, deathrate=d_2, μ=mu, ploidy=1)

simulation1 = runsimulation(SimpleTreeCell, input11)

simulation2 = runsimulation(simulation1.output, SimpleTreeCell, WellMixed, input44)
simulation2 = runsimulation(simulation1, SimpleTreeCell, WellMixed, input22)

simulation3 = runsimulation(SimpleTreeCell, WellMixed, input22)

timesteps = [5.0,10.0]

function func(population)
    return 0
end

simulation3 = runsimulation_timeseries_returnfinalpop(SimpleTreeCell, input22, timesteps, func)


