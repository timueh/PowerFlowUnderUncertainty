export 	createCSV,
		createTikz,
		createLabels,
		createDictForTikz,
		createTikzFigure


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

function createTikz(fileDict,data,foldername)
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

function createTikz(fileDict,data,foldername,compfoldername)
	for (key, val) in fileDict
		display(val)
		for (i, samples) in enumerate(eachrow(data[key]))
			fname = val["name"]*"_"*string(i)*".csv"
			xlab, ylab = createLabels(key,i)
			d = createDictForTikz(xlab,ylab,foldername*fname,compfoldername*fname; compcolor=val["compcolor"], color=val["color"],width=val["width"],height=val["height"])
			fname = val["name"]*"_"*string(i)*".tex"
			createTikzFigure(fname,d,source="output/bareTikzFile_comp.tex")
		end
	end
end

function createLabels(x::Symbol,i::Int)
	xlab, ylab = [], []
	if x == :i
		xlab = string("\\ibus_{",i,"}")
		ylab = string("\\density_{\\iRV_{",i,"}}(",xlab,")")
	elseif x == :v
		xlab = string("\\vmag_{",i,"}")
		ylab = string("\\density_{\\vmagRV_{",i,"}}(",xlab,")")
	elseif x == :θ
		xlab = string("\\phase_{",i,"}")
		ylab = string("\\density_{\\phaseRV_{",i,"}}(",xlab,")")
	elseif x == :pg
		# assign generator number to bus number!
		i == 1 ? i = 1 : i = 3
		xlab = string("\\pbus_{",i,"}")
		ylab = string("\\density_{\\pbusRV_{",i,"}}(",xlab,")")
	elseif x == :qg
		# assign generator number to bus number!
		i == 1 ? i = 1 : i = 3
		xlab = string("\\qbus_{",i,"}")
		ylab = string("\\density_{\\qbusRV_{",i,"}}(",xlab,")")
	elseif x == :pd
		i == 1 ? i = 2 : i = 4
		xlab = string("\\pbus_{",i,"}")
		ylab = string("\\density_{-\\pbusRV_{",i,"}}(",xlab,")")
	elseif x == :qd
		i == 1 ? i = 2 : i = 4
		xlab = string("\\qbus_{",i,"}")
		ylab = string("\\density_{-\\qbusRV_{",i,"}}(",xlab,")")
	elseif x == :pl_t
		xlab = string("\\pbranch_{",i,"}")
		ylab = string("\\density_{\\pbranchRV_{",i,"}}(",xlab,")")
	elseif x == :ql_t
		xlab = string("\\qbranch_{",i,"}")
		ylab = string("\\density_{\\qbranchRV_{",i,"}}(",xlab,")")
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

function createDictForTikz(x::String,y::String,fname::String,compfname::String; compcolor="gray!20", color="black!20", width="5cm", height="3cm")
	# keys = ["xlab", "ylab", "filename", "color", "width", "height"]
	D = Dict("xlab"=>	"xlabel = \$"*x*"\$,",
			"ylab"=>	"ylabel = \$"*y*"\$,",
			"filename"=> "file {"*fname*"};",
			"filename_comp"=> "file {"*compfname*"};",
			"color" => "color = "*color,
			"color_comp" => "color = "*compcolor,
			"width" => "width = "*width*",",
			"height" => "height = "*height*",")
end

function createTikzFigure(fname::String,d::Dict; source="output/bareTikzFile.tex")
	newfile = open(fname,"w")
	open(source) do file
	           for ln in eachline(file)
	           		write(newfile, haskey(d, strip(ln)) ? d[strip(ln)]*"\n" : ln*"\n")
	           end
	       end
	close(newfile)
end