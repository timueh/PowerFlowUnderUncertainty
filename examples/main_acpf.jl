using PowerFlowUnderUncertainty, LinearAlgebra, JuMP, Ipopt, DelimitedFiles
include("powersystem.jl")
include("init.jl")

pf = Model(with_optimizer(Ipopt.Optimizer))
addCore!(pf,sys,unc)
addPVBus!(pf,sys,unc)
optimize!(pf)

pf_state = getGridState(pf,sys,unc)
pf_samples = generateSamples(ξ,pf_state,sys,unc)

mycolor = "red"

# plotHistogram_gen(pf_samples[:pd], "pd"; fignum = 1+10, color = mycolor)
# plotHistogram_gen(pf_samples[:qd], "qd"; fignum = 2+10, color = mycolor)

width, height, color = "3.9cm", "2.75cm", "black!40"
files_to_save = Dict(:v => Dict("name"=>"voltage_magnitude",
								"color"=>color,
								"width"=>"2.95cm",
								"height"=>height),
					 :θ => Dict("name"=>"voltage_angle",
								"color"=>color,
								"width"=>"2.95cm",
								"height"=>height),
					 :pg => Dict("name"=>"active_power",
								"color"=>color,
								"width"=>width,
								"height"=>height),
					 :qg => Dict("name"=>"reactive_power",
								"color"=>color,
								"width"=>width,
								"height"=>height),
					 :pd => Dict("name"=>"active_power_demand",
								"color"=>color,
								"width"=>"4.7cm",
								"height"=>height),
					 :qd => Dict("name"=>"reactive_power_demand",
								"color"=>color,
								"width"=>"4.7cm",
								"height"=>height),
					 :i => Dict("name"=>"current_magnitude",
								"color"=>color,
								"width"=>width,
								"height"=>height),
					 :pl_t => Dict("name"=>"active_power_to",
								"color"=>color,
								"width"=>width,
								"height"=>height),
					 :ql_t => Dict("name"=>"reactive_power_to",
								"color"=>color,
								"width"=>width,
								"height"=>height))

createCSV(files_to_save,pf_samples)
createTikz(files_to_save,pf_samples,"")
# mv ~/Documents/JuliaDev/PowerFlowUnderUncertainty/examples/*.tex ./
# mv ~/Documents/JuliaDev/PowerFlowUnderUncertainty/examples/*.csv ./



