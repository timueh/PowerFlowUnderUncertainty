using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt
include("powersystem.jl")
sys = setupPowerSystem()

μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."
deg = 4
unc = setupUncertainty(μ,σ,w,sys[:Nd],deg)
ξ = sampleFromGaussianMixture(5000,μ,σ,w)

pf = Model(with_optimizer(Ipopt.Optimizer))
addCore!(pf,sys,unc)
addPVBus!(pf,sys,unc)
optimize!(pf)

pf_state = getGridState(pf,sys,unc)
pf_samples = generateSamples(ξ,pf_state,unc)

plotHistogram_gen(pf_samples[:pg], "pg"; fignum = 1, color = "green")
plotHistogram_gen(pf_samples[:qg], "qg"; fignum = 2)
plotHistogram_v(pf_samples[:v], "v"; fignum = 3, color = "cyan")
plotHistogram_i(pf_samples[:i], "i"; fignum = 4)
