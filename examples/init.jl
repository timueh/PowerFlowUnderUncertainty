function createCSV(fileDict,data::Dict)
	mynames = String[]
	for (key, val) in fileDict
		for (i, samples) in enumerate(eachrow(data[key]))
			push!(mynames, val["name"]*"_"*string(i))
			fname = val["name"]*"_"*string(i)*".csv"
			open(fname, "w") do io
				writedlm(fname, samples, ',')
			end
		end
	end
end

function createTikz(fileDict,data,foldername="figures/CaseStudy/")
	for (key, val) in fileDict
		display(val)
		for (i, samples) in enumerate(eachrow(data[key]))
			fname = val["name"]*"_"*string(i)*".csv"
			# x, xsub = string(key), string(key)*"_"*string(i)
			# d = createDictForTikz(xsub,"\\density_{\\rv{"*x*"}_"*string(i)*"}("*xsub*")",foldername*fname; color=val["color"],width=val["width"],height=val["height"])
			xlab, ylab = createLabels(key,i)
			d = createDictForTikz(xlab,ylab,foldername*fname; color=val["color"],width=val["width"],height=val["height"])
			# *x*"\}_"*string(i)*"\}("*xsub*")"
			fname = val["name"]*"_"*string(i)*".tex"
			createTikzFigure(fname,d)
		end
	end
end

function createLabels(x::Symbol,i::Int)
	xlab, ylab = [], []
	if x in [:i]
		xlab = string("\\ibus_{",i,"}")
		ylab = string("\\density_{\\iRV_{",i,"}}(",xlab,")")
	elseif x == :v
		xlab = string("\\vmag_{",i,"}")
		ylab = string("\\density_{\\vmagRV_{",i,"}}(",xlab,")")
	elseif x in [:θ]
		xlab = string("\\phase_{",i,"}")
		ylab = string("\\density_{\\phaseRV_{",i,"}}(",xlab,")")
	elseif x in [:pg]
		xlab = string("\\p_{",i,"}")
		ylab = string("\\density_{\\pRV_{",i,"}}(",xlab,")")
		# xlab = string("p_{",i,"}^{\\text{g}}")
		# ylab = string("\\density_{\\rv{p}_{",i,"}}^{\\text{g}}(",xlab,")")
	elseif x in [:qg]
		xlab = string("\\q_{",i,"}")
		ylab = string("\\density_{\\qRV_{",i,"}}(",xlab,")")
		# xlab = string("q_{",i,"}^{\\text{g}}")
		# ylab = string("\\density_{\\rv{q}_{",i,"}}^{\\text{g}}(",xlab,")")
	else
		throw(error("Symbol $x is not known."))
	end
	xlab, ylab
end

function createDictForTikz(x::String,y::String,fname::String; color="black!20", width="5cm", height="3cm")
	# keys = ["xlab", "ylab", "filename", "color", "width", "height"]
	D = Dict("xlab"=>	"xlabel = \$"*x*"\$,",
			"ylab"=>	"ylabel = \$"*y*"\$,",
			"filename"=> "file {"*fname*"};",
			"color" => "color = "*color,
			"width" => "width = "*width*",",
			"height" => "height = "*height*",")
end




function createTikzFigure(fname::String,d::Dict)
	newfile = open(fname,"w")
	open("output/bareTikzFile.txt") do file
	           for ln in eachline(file)
	           		write(newfile, haskey(d, strip(ln)) ? d[strip(ln)]*"\n" : ln*"\n")
	           end
	       end
	close(newfile)
end

width, height = "3.9cm", "2.75cm" 

sys = setupPowerSystem()

μ, σ, w = [2.1, 3.2], [0.3, 0.4], [0.3, 0.7]
@assert sum(w) == 1 "The weights do not sum to one."
deg = 4
unc = setupUncertainty(μ,σ,w,sys[:Nd],deg)
ξ = sampleFromGaussianMixture(5000,μ,σ,w)

