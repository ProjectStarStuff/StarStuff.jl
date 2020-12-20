using DrWatson
@quickactivate "crprop"

module PlotCSV

using DrWatson
using AbstractPlotting
using AbstractPlotting.MakieLayout
using CSV
using DataFrames
using Statistics
using Unitful
using UnitfulAstro

function plot_csv(filename::String)
	df = DataFrame(CSV.File(filename))

	coordnames = [:x,:y,:z]
	centernames = []

	for name ∈ coordnames
		centername = Symbol(string(name,"center"))
		push!(centernames,centername)
		df[centername] = ustrip.(u"pc",(df[!,name] .- mean(df[!,name]))*u"m")
		df[name] = ustrip.(u"pc",df[!,name]*u"m")
	end

	# lines(df[:xcenter],df[:ycenter],df[:zcenter])
	lines(df[:x],df[:y],df[:z])
end
export plot_csv

function plot_csv!(scene::Scene,filename::String)
	df = DataFrame(CSV.File(filename))

	coordnames = [:x,:y,:z]
	centernames = []

	for name ∈ coordnames
		centername = Symbol(string(name,"center"))
		push!(centernames,centername)
		df[centername] = ustrip.(u"pc",(df[!,name] .- mean(df[!,name]))*u"m")
		df[name] = ustrip.(u"pc",df[!,name]*u"m")
	end

	# lines(df[:xcenter],df[:ycenter],df[:zcenter])
	lines!(scene,df[:x],df[:y],df[:z])
end
export plot_csv!


end