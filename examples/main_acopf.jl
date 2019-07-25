using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt
include("powersystem.jl")
include("init.jl")

opf = Model(with_optimizer(Ipopt.Optimizer, max_iter=1000))
addCore!(opf,sys,unc)
addConsistencyConstraints!(opf)
# addGenerationConstraint!(1,:pg,opf,sys,unc)
# addGenerationConstraint!(2,:qg,opf,sys,unc)
# addVoltageConstraint!(2,opf,sys,unc)
addCurrentConstraint!(5,opf,sys,unc)
addCost!(opf,sys,unc)
optimize!(opf)

opf_state = getGridState(opf,sys,unc)
opf_samples = generateSamples(ξ,opf_state,unc)

mycolor = "red"
plotHistogram_gen(opf_samples[:pg], "pg"; fignum = 1, color = mycolor)
plotHistogram_gen(opf_samples[:qg], "qg"; fignum = 2, color = mycolor)
plotHistogram_nodal(opf_samples[:v], "v"; fignum = 3, color = mycolor)
plotHistogram_nodal(opf_samples[:θ], "θ"; fignum = 4, color = mycolor)
plotHistogram_branch(opf_samples[:i], "i"; fignum = 5, color = mycolor)
plotHistogram_branch(opf_samples[:pl_t], "pl_t"; fignum = 6, color = mycolor)
