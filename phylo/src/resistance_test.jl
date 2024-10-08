
using SomaticEvolution, Random, AbstractTrees

# Parameters before treatment
b_s1 = 1.0
d_s1 = 0.0
b_r1 = 1.0
d_r1 = 0.0

# Parameters during treatment
b_s2 = 1.0
d_s2 = 3.0
b_r2 = 1.0
d_r2 = 0.0

# rate of neutral mutations
mu = 2.0
omega = 2*mu

# Define input variables for simulations
mutant_time = [1.0, 2.0, 3.0]
mutant_selection = [0.0, 0.0, 0.0]

mutant_time2 = [Inf, Inf, Inf]
mutant_selection2 = [0.0, 0.0, 0.0]

# other parameters
Nd = 1000
Nmax_vector = Nd*[1.1, 1.5]
measurements = length(Nmax_vector) + 1

realizations = 20

# Grow to maximal population size
input = BranchingInput(Nmax=Nd, birthrate=b_s1, deathrate=d_s1, μ=mu, ploidy=1)
selection = SelectionPredefined(mutant_selection, mutant_time)
simulation = runsimulation(SimpleTreeCell, input, selection)

input2 = BranchingInput(Nmax=2.0*Nd, birthrate=b_s2, deathrate=d_s2, μ=mu, ploidy=1)
selection2 = SelectionPredefined(mutant_selection2, mutant_time2)
simulation = runsimulation(simulation.output, SimpleTreeCell, WellMixed, input2)


simulation.output.subclones[1].birthrate
simulation.output.subclones[1].deathrate
simulation.output.subclones[2].birthrate
simulation.output.subclones[2].deathrate
simulation.output.subclones[3].birthrate
simulation.output.subclones[3].deathrate