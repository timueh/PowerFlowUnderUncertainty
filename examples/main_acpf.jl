using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles
include("powersystem.jl")
include("init.jl")

pf = Model(with_optimizer(Ipopt.Optimizer))
addCore!(pf,sys,unc)
addPVBus!(pf,sys,unc)
optimize!(pf)

pf_state = getGridState(pf,sys,unc)
pf_samples = generateSamples(ξ,pf_state,unc)

mycolor = "red"
plotHistogram_gen(pf_samples[:pg], "pg"; fignum = 1, color = mycolor)
plotHistogram_gen(pf_samples[:qg], "qg"; fignum = 2, color = mycolor)
plotHistogram_nodal(pf_samples[:v], "v"; fignum = 3, color = mycolor)
plotHistogram_nodal(pf_samples[:θ], "θ"; fignum = 4, color = mycolor)
plotHistogram_branch(pf_samples[:i], "i"; fignum = 5, color = mycolor)
plotHistogram_branch(pf_samples[:pl_t], "pl_t"; fignum = 6, color = mycolor)
plotHistogram_gen(pf_samples[:pd], "pd"; fignum = 1+10, color = mycolor)
plotHistogram_gen(pf_samples[:qd], "qd"; fignum = 2+10, color = mycolor)

width, height = "3.9cm", "2.75cm" 
files_to_save = Dict(:v => Dict("name"=>"voltage_magnitude",
								"color"=>"blue",
								"width"=>"2.95cm",
								"height"=>height),
					 :θ => Dict("name"=>"voltage_angle",
								"color"=>"blue",
								"width"=>"2.95cm",
								"height"=>height),
					 :pg => Dict("name"=>"active_power",
								"color"=>"blue",
								"width"=>width,
								"height"=>height),
					 :qg => Dict("name"=>"reactive_power",
								"color"=>"blue",
								"width"=>width,
								"height"=>height),
					 :pd => Dict("name"=>"active_power_demand",
								"color"=>"blue",
								"width"=>"4.7cm",
								"height"=>height),
					 :qd => Dict("name"=>"reactive_power_demand",
								"color"=>"blue",
								"width"=>"4.7cm",
								"height"=>height),
					 :i => Dict("name"=>"current_magnitude",
								"color"=>"black!20",
								"width"=>width,
								"height"=>height))

createCSV(files_to_save,pf_samples)
createTikz(files_to_save,pf_samples,"figures/CaseStudy/ppf/")




