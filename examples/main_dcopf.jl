using PowerFlowUnderUncertainty, LinearAlgebra, MosekTools, JuMP
include("powersystem.jl")
sys = setupPowerSystem()

μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."
deg = 4
unc = setupUncertainty(μ,σ,w,sys[:Nd],deg)
ξ = sampleFromGaussianMixture(5000,μ,σ,w)

opf = Model(with_optimizer(Mosek.Optimizer))
addCoreDC!(opf,sys,unc)
addConstraintsGenDC!(opf,sys,unc)
addConstraintsLineFlowDC!(opf,sys,unc)
addCostDC!(opf,sys,unc)
optimize!(opf)

# opf_state = getGridState(opf,sys,unc)
# opf_samples = generateSamples(ξ,opf_state,unc)
#
# plotHistogram_gen(opf_samples[:pg], "pg"; fignum = 1, color = "green")
# plotHistogram_gen(opf_samples[:qg], "qg"; fignum = 2)
# plotHistogram_v(opf_samples[:v], "v"; fignum = 3, color = "cyan")
# plotHistogram_i(opf_samples[:i], "i"; fignum = 4)
