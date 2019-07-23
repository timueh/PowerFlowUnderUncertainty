using PowerFlowUnderUncertainty, LinearAlgebra, MosekTools, JuMP, PolyChaos
include("powersystem.jl")
sys = setupPowerSystem()

μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."
deg = 4
unc = setupUncertainty(μ,σ,w,sys[:Nd],deg)
ξ = sampleFromGaussianMixture(5000,μ,σ,w)

dcopf = Model(with_optimizer(Mosek.Optimizer))
addCoreDC!(dcopf,sys,unc)
addConstraintsGenDC!(dcopf,sys,unc)
addConstraintsLineFlowDC!(dcopf,sys,unc)
addCostDC!(dcopf,sys,unc)
optimize!(dcopf)

dcopf_state = getGridStateDC(dcopf,sys,unc)
dcopf_samples = generateSamples(ξ,dcopf_state,unc)

plotHistogram_gen(dcopf_samples[:pg], "pg"; fignum = 1, color="red", alpha=0.3)
plotHistogram_nodal(dcopf_samples[:θ], "θ"; fignum = 4, color="red", alpha=0.3)
plotHistogram_branch(dcopf_samples[:pl], "pl"; fignum = 6, color="red", alpha=0.3)

T2 = Tensor(2, unc[:opq])
psol = value.(dcopf[:pg])
cost1 = sum([ sys[:costquad][i]*sum(T2.get([k-1,k-1])*psol[i,k]^2 for k in 1:unc[:dim]) + sys[:costlin][i]*psol[i,1] for i in 1:sys[:Ng] ])
lin, quad = buildCostSOC(sys[:costlin],sys[:costquad],T2,deg+1)
cost2 = objective_value(dcopf)^2 - sum( lin[i]^2/quad[i] for i in 1:length(lin) )
