using DrWatson
@quickactivate "crprop"

module GUItest

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout

function test()
	# Parent scene
	outerPadding = 30
	parentScene, parentLayout = layoutscene(outerPadding, resolution = (1600,900), backgroundcolor = RGBf0(0.98, 0.98, 0.98))

	# UI Setup
	#	ui grid
	uiGrid = GridLayout()
	#	put in first 2 rows of first column
	parentLayout[1:2,1] = uiGrid


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

	function clear!(scene::Scene) 
		if length(scene.plots) > 2
			delete!(scene,scene.plots[end])		
			clear!(scene)
		elseif length(scene.plots) == 2
			delete!(scene,scene.plots[end])		
		end
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

export test

end # GUItest