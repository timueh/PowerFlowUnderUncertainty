using PowerFlowUnderUncertainty, LinearAlgebra, MosekTools, JuMP, PolyChaos, DelimitedFiles
include("powersystem.jl")
include("init.jl")

dcopf = Model(with_optimizer(Mosek.Optimizer))
addCoreDC!(dcopf,sys,unc)
addConstraintsGenDC!(dcopf,sys,unc)
# addConstraintsLineFlowDC!(dcopf,sys,unc)
addCostDC!(dcopf,sys,unc)
optimize!(dcopf)

dcopf_state = getGridStateDC(dcopf,sys,unc)
dcopf_samples = generateSamples(ξ,dcopf_state,sys,unc)

########################################################################
##### POST PROCESSING #####
########################################################################

mycolor = "green"
plotHistogram_gen(dcopf_samples[:pg], "pg"; fignum = 1, color=mycolor, alpha=0.3)
plotHistogram_nodal(dcopf_samples[:θ], "θ"; fignum = 4, color=mycolor, alpha=0.3)
plotHistogram_branch(dcopf_samples[:pl_t], "pl"; fignum = 6, color=mycolor, alpha=0.3)

width, height, color, compcolor = "3.9cm", "2.75cm", "black!60!green", "red!20"
files_to_save = Dict(:θ => Dict("name"=>"voltage_angle",
								"color"=>color,
								"compcolor"=>compcolor,
								"width"=>"2.95cm",
								"height"=>height),
					 :pg => Dict("name"=>"active_power",
								"color"=>color,
								"compcolor"=>compcolor,
								"width"=>"4.7cm",
								"height"=>height),
					 :pl_t => Dict("name"=>"active_power_to",
								"color"=>color,
								"compcolor"=>compcolor,
								"width"=>width,
								"height"=>height))

createCSV(files_to_save,dcopf_samples)
createTikz(files_to_save,dcopf_samples,"","../acopf_con/")
