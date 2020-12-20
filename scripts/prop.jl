using DrWatson
@quickactivate "crprop"

module Prop

using DrWatson
using Reexport
using Statistics
#using FITSIO

include(srcdir("CrProp.jl"))
@reexport using .CrProp

"""
	propRun()

Main function to propagate cosmic rays
"""
function propRun()
	# make an electron and a proton
	electron = ParticleID(0,0,1)
	proton = ParticleID(1,0,0)
	particleTypes = [electron,proton]

	# energies at which to generate particles
	energies = [1.0,10.0,100.0,1000.0]*u"GeV"
	energies = upreferred.(energies)

	# fire particles from different galactic radii
	rStart = [1.0,5.0,8.5,10.0,15.0]*u"kpc"
	rStart = upreferred.(rStart)

	# directions to fire particles (these will be normalized)
	xGal = [1.0,0.0,0.0]
	yGal = [0.0,1.0,0.0]
	zGal = [0.0,0.0,1.0] 
	outAlongDisk = [1.0,1.0,0.0]
	inAlongDisk  = [-1.0,-1.0,0.0]
	outRadially  = [1.0,1.0,1.0]
	particleAim = [xGal,yGal,zGal,outAlongDisk,inAlongDisk,outRadially]

	gpFrac = 0.05

	# define electric and magnetic fields
	efield(x::Array{T}) where {T<:Number} = vec([0.0, 0.0, 0.0])
	bfield = Bfield.model[:galacticDipole]



	# for construction ParticleTree and input to integrator
	params = Dict(
			:t => 0.0,
			:efield => efield,
			:bfield => bfield
		)


	initialPT = initParticles(
							 particleTypes,
							 energies,
							 rStart,
							 particleAim,
							 gpFrac,
							 bfield
							)

	
	solPT = nothing
	# iterate over Timescales
	for i = 1:length(initialPT.nodes)
		# Build problem using boris algorithm integrator and initial Snapshot
		# of each Timescale
		prob = DiscreteProblem(rintegrator[:boris], 
			initialPT[i,1], (0.0,initialPT.nodes[i].dt*1000), params)

		# Solve problem using time step specified in Timescale node
		solution = solve(prob, FunctionMap(), dt = initialPT.nodes[i].dt)
		# Construct a new Timescale out of the solution and add to solution particle tree
		solTS = construct(Timescale,solution.u,Float64[],solution.t,initialPT.nodes[i].dt)
		if solPT == nothing
			solPT = construct(ParticleTree, Vector([solTS]))
		else
			add_node!(solPT,solTS)
		end
	end

	tVec = Vector()

	for i=1:6
		
		times, traj = getParticleTrajectory(solPT,20,1,i)
		saveMe = Dict(
			:x => [coords.x[1][1] for coords in traj],
			:y => [coords.x[1][2] for coords in traj],
			:z => [coords.x[1][3] for coords in traj]
			)
		push!(tVec,saveMe)
	end

	"""
		scaleShift(a::T) where {T<:AbstractVector}

	Shifts values of array so it is centered on its mean and range is normalized to 1
	"""
	function scaleShift(a::T) where {T<:AbstractVector}
		aRange = max(a...) - min(a...)
		aMean  = mean(a)
		out = a .- aMean
		if aRange == 0.0
			return out
		else
			return out ./ aRange
		end
	end

	# Plotting stuff

	# Parent scene
	outerPadding = 30
	parentScene, parentLayout = layoutscene(outerPadding, resolution = (1600,900), backgroundcolor = RGBf0(0.98, 0.98, 0.98))

	# UI Setup
	#	ui grid
	uiGrid = GridLayout()
	#	put in first 2 rows of first column
	parentLayout[1:2,1] = uiGrid

	colsize!(parentLayout,1,Fixed(150))

	#	first slider at top of uiGrid
	startVal = 2
	ls = labelslider!(parentScene,"# helices",1:10;format = x->"count = $(x)",startvalue = startVal)
	uiGrid[1,1] = ls.layout
	sl1 = uiGrid[2,1] = ls.slider
	sl1node = sl1.value

	# Run button on the bottom
	runButton = LButton(parentScene,label="RUN")

	# 3D Scene
	#	Place scene to span all rows of the second column
	plot3dScene = parentLayout[:, 2] = LScene(parentScene, scenekw = (camera = cam3d!, raw = false))

	#	Get axes from the scene
	plot3dAxis  = plot3dScene.scene[Axis]

	#	helix function for testing
	function helix(c::Number,start::Number,step::Number,stop::Number)
		zz = collect(start:step:stop)
		ϕ  = zz .+ c
		xx = cos.(ϕ)
		yy = sin.(ϕ)
		return [xx,yy,zz]
	end


	#	Initial plot
	tIndex = startVal
	c = 2*π/tIndex
	for i =0:tIndex
		x,y,z = helix(i*c,0.0,π/8,8*π)
		lines!(plot3dScene, x,y,z)#scaleShift.(collect(values(tVec[i])))...)
	end



	on(runButton.clicks) do clicks
		clear!(plot3dScene.scene)
		tIndex = to_value(sl1node)

		c = 2*π/tIndex
		for i =0:tIndex
			x,y,z = helix(i*c,0.0,π/8,8*π)
			lines!(plot3dScene, x,y,z)#scaleShift.(collect(values(tVec[i])))...)
		end
	end

	return parentScene

end

export propRun
end # Prop

import .Prop

Prop.propRun()
