
module PointNshoot
include("StarStuff.jl")

using .StarStuff
using GLMakie

"""
        singleParticleTree(particle::ParticleID,energy::Quantity,r0::Quantity,direction::Vector,gpFrac::Float64)

Build particle tree for a single particle
"""
function singleParticleTree(particle::ParticleID,energy::Quantity,r0::Quantity,direction::Vector,gpFrac::Float64)
    # TODO: allow fields as input
    # define electric and magnetic fields
    efield(x::Array{T}) where {T<:Number} = vec([0.0, 0.0, 0.0])
    bfield = Bfield.model[:jansson2012m]

    # generate initial particle tree
    initialPT = StarStuff.initParticles(
        Vector([particle]),
        Vector([upreferred(energy)]),
        Vector([upreferred(r0)]),
        Vector([direction]),
        gpFrac,
        bfield
    )

    # for construction ParticleTree and input to integrator
    params = Dict(
                  :t => 0.0,
                  :efield => efield,
                  :bfield => bfield
                 )

    # INTEGRATION
    # declare variable name for solution particle tree
    solPT = nothing
    # iterate over Timescales
    for i = 1:length(initialPT.nodes)
        # Build problem using boris algorithm integrator and initial Snapshot
        # of each Timescale
        prob = DiscreteProblem(rintegrator[:boris], 
                               initialPT[i,1], (0.0,initialPT.nodes[i].dt*100000), params)
        println("test $i")
        println(typeof(initialPT[i,1]))
        println(typeof(prob))
        # Solve problem using time step specified in Timescale node
        solution = solve(prob, FunctionMap(), dt = initialPT.nodes[i].dt, progress = false)
        println("test $(i) end")
        # Construct a new Timescale out of the solution and add to solution particle tree
        solTS = construct(Timescale,solution.u,Float64[],solution.t,initialPT.nodes[i].dt)
        if solPT === nothing
            solPT = construct(ParticleTree, Vector([solTS]))
        else
            add_node!(solPT,solTS)
        end
    end
    return solPT
end


function best_units(x::Float64)
    # order of magnitude
    oom = mag10(x)
    units = nothing
    # TODO: find method to get strings from units
    if oom < 3
        units = (u"m","m")
    elseif oom < 6
        units = (u"km","km")
    elseif oom < 9
        units = (u"Mm","Mm")
    elseif oom < 14
        units = (u"AU","AU")
    elseif oom < 20
        units = (u"pc","pc")
    else
        units = (u"kpc","kpc")
    end
    return units
end

function show_bfield!(scene::Scene)
    println("Entered show_bfield!")
    # units for all distances
    units_distance = u"kpc"

    # number of points along radial coordinate
    nr = 10
    # innermost radius is 1 kpc IAW Jansson2012
    rmin = 1.0u"kpc"
    # take outermost radius from disk limits in Jansson field
    rmax = Bfield.janssonDict[:rdisklims][end]

    # number of points along z-coordinate
    nz = 20
    # set z-coordinate boundaries to be ±10 kpc for now to bring in some of the halo field
    zmin = -10.0u"kpc"
    zmax = 10.0u"kpc"

    # generate r and z coordinate arrays
    rr = LinRange(rmin,rmax,nr)
    zz = LinRange(zmin,zmax,nz)

    # array to store coordinates
    coords = nothing

    # number of points around innermost radius
    nmin = 5
    # number of points around outermost radius
    nmax = 20    
    # generate cartesian coordinates and store in coords
    for (index,r) in enumerate(rr)
        # calculate number of points along this radius
        frac = (nmax-nmin)*index/nr
        n = nmin + Int(round(frac))
        println("Number of points around radius $index = ", n)
        # generate polar angle coordinates
        θ = LinRange(0.0,2*π,n)
        # transform to cartesian
        yy = r*sin.(θ)
        xx = r*cos.(θ)
        # make array for storing coordinates
        nxnynz = length(xx)*length(yy)*length(zz)
        tempcoords = zeros(Float64,(nxnynz,3))*units_distance
        # iterate over all combinations of generated coordinates to form grid
        ii = 1
        for xi in xx
            for yi in yy
                for zi in zz
                    tempcoords[ii,:] .= xi,yi,zi
                    ii+=1
                end
            end
        end
        # concatenate new coordinates 
        if index == 1
            coords = tempcoords
        else
            coords = vcat(coords,tempcoords)
        end
    end
    # xx = nothing
    # yy = nothing
    # for (index,r) in enumerate(rr)
    #     # calculate number of points along this radius
    #     frac = (nmax-nmin)*index/nr
    #     n = nmin + Int(round(frac))
    #     println("Number of points around radius $index = ", n)
    #     # generate polar angle coordinates
    #     θ = LinRange(0.0,2*π,n)
    #     # transform to cartesian
    #     y = r*sin.(θ)
    #     x = r*cos.(θ)
    #     if index == 1
    #         yy = y
    #         xx = x
    #     else
    #         append!(yy,y)
    #         append!(xx,x)
    #     end
    # end

    println("Made all spatial coordinates. Calculating magnetic field.")
    # generate array to store magnetic field values
    # bunits = Unitful.unit(Bfield.jansson2012(coords[1,:])[1])
    # bb = zeros(Float64,size(coords))*bunits
    # calculate magnetic field values
    # for i = 1:size(coords)[1]
    #     bb[i,:] .= Bfield.jansson2012(coords[i,:])
    # end
    # println("Finished calculating magnetic field")
    janssonx(x,y,z) = Bfield.jansson2012kpc([x,y,z])[1]*1.0e10
    janssony(x,y,z) = Bfield.jansson2012kpc([x,y,z])[1]*1.0e10
    janssonz(x,y,z) = Bfield.jansson2012kpc([x,y,z])[1]*1.0e10

    coords = ustrip.(units_distance, coords)

    bz = [janssonz(v...) for v ∈ eachrow(coords)]
    println(bz[1:10])
    println(min(bz...))
    println(max(bz...))

    scatter!(scene,coords[:,1],coords[:,2],coords[:,3],color=bz*10,markersize=100)
    # volume!(scene,coords[:,1],coords[:,2],coords[:,3],jansson1,algorithm=:mip)
    # volume!(scene,rand(32,32,32), algorithm=:mip)
end

function pnsScene!(scene::Scene,pt::ParticleTree;filename::Union{String,Nothing}=nothing)
    # Extract coordinates from particle tree
    times, traj = getParticleTrajectory(pt,1,1,1)
    # println("dt = ", times[51]-times[50])
    # println("steps = ",length(times))
    # println("max time = ", times[end])

    xx = [c.x[1][1] for c in traj]
    yy = [c.x[1][2] for c in traj]
    zz = [c.x[1][3] for c in traj]
    coords = [xx,yy,zz]

    #out_df = DataFrame(t=times,x=xx,y=yy,z=zz)
    #CSV.write()

    # Shift coordinates so they're ceneterd at zero
    μcoords = mean.(coords)
    shiftedCoords = [coords[i] .- μcoords[i] for i=1:3]
    localOrigin = μcoords*u"m"
    localOrigin = uconvert.(u"kpc",localOrigin)

    # Calculate bounds for plotting
    xtreme,ystreme,zstreme = extrema.(shiftedCoords)
    mincoord = min(xtreme[1],ystreme[1],zstreme[1])
    maxcoord = max(xtreme[2],ystreme[2],zstreme[2])

    bounds = FRect3D([mincoord,mincoord,mincoord],[maxcoord,maxcoord,maxcoord])

    # 3D SCENE
    makelabels(xs,tickunits) = ["$(@sprintf("%.2e",ustrip(tickunits,x*u"m")))" for x in xs]	
    maxticks = collect(LinRange(mincoord,0.0,5))
    append!(maxticks,collect(LinRange(0.0,maxcoord,5)[2:end]))
    tickunits,axesunits = best_units(maxcoord-mincoord)
    ticklabels = makelabels(maxticks,tickunits)
    #	Plot coordinates and labels
    lines!(scene, shiftedCoords..., scale_plot=false,limits = bounds)
    xticks!(scene; xtickrange=maxticks, xticklabels=ticklabels)
    yticks!(scene; ytickrange=maxticks, yticklabels=ticklabels)
    zticks!(scene; ztickrange=maxticks, zticklabels=ticklabels)
    xlabel!(scene, "X ($axesunits)")
    ylabel!(scene, "Y ($axesunits)")
    zlabel!(scene, "Z ($axesunits)")
    AbstractPlotting.update!(scene)

    return scene
end



"""
        run()

Run PointNshoot script
"""
function run()
    GLMakie.activate!()
    scene = Scene()
    show_bfield!(scene)
    return scene
end

export run

end # PointNshoot

import .PointNshoot
PointNshoot.run()