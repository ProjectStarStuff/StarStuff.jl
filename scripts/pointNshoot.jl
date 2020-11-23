module PointNshoot

using Blink
using Interact

using StarStuff
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

    return scene
end

"""
    function particlepicker()

Widget to select a particle. Returns a ParticleID
"""
function particlepicker()
    el = spinbox(label="Electrons"; value=0)
    pr = spinbox(label="Protons"; value=0)
    nu = spinbox(label="Neutrons"; value=0)
    output = Interact.@map(ParticleID(&el,&pr,&nu))
    wdg = Widget(["e"=>el,"p"=>pr,"n"=>nu],output=output)
    @layout! wdg Interact.hbox(el,pr,nu)
end


"""
        run()

Run PointNshoot script
"""
function run()
    GLMakie.activate!()

    s = Observable{Any}(Scene())

    sine_button = Interact.button("Sine plot")
    cos_button = Interact.button("Cos plot")
    spi_button = Interact.button("Spiral plot")
    clear_button = Interact.button("Clear plots")

    buttons = Interact.vbox(sine_button,cos_button,spi_button,clear_button)


    function sineplot!(a::Int64,s::Scene)
        t = collect(0:0.1:10)
        # x = cos.(2*π*t)
        y = sin.(2*π*t)
        # z = t
        println("Sin plot")
        lines!(s,t,y)
        update!(s)
    end

    function cosplot!(a::Int64,s::Scene)
        t = collect(0:0.1:10)
        x = cos.(2*π*t)
        # y = sin.(2*π*t)
        # z = t
        println("Cos plot")
        lines!(s,t,x)
        update!(s)
    end

    function spiplot!(a::Int64,s::Scene)
        t = collect(0:0.1:10)
        x = cos.(2*π*t)
        y = sin.(2*π*t)
        z = t
        println("Spiral plot")
        lines!(s,x,y,z)
        update!(s)
    end

    function clearplot!(a::Int64,s::Scene)
        clear!(s)
    end



    map!(sineplot!,s,sine_button,s[])
    map!(cosplot!,s,cos_button,s[])
    map!(spiplot!,s,spi_button,s[])
    map!(clearplot!,s,clear_button,s[])

    ui = dom"div"(buttons)
    w  = Window()
    body!(w, ui)

    return s[]
    # # REQUEST SIMULATION PARAMETERS
    # # get number of particles
    # nparticles = get_until_correct("Number of particles: ",Int)

    # # immediately get parameters if nparticles=1
    # if nparticles == 1
    #     # particle type
    #     particles = Dict(
    #                      "e-" => ParticleID(0,0,1),
    #                      "p+"   => ParticleID(1,0,0),
    #                      "He+2" => ParticleID(2,2,0)
    #                     )
    #     particlekeys = Vector(collect(keys(particles)))
    #     particlemenu = request("Select a particle:",RadioMenu(particlekeys))
    #     particle = particles[particlekeys[particlemenu]]

    #     println("particle = ", particle)
    #     # energy
    #     eunits = ["MeV","GeV","TeV"]
    #     eunitsmenu = request("Select energy units:",RadioMenu(eunits))
    #     eunit = eunits[eunitsmenu]
    #     eprompt = "Enter energy value ($eunit):"
    #     eval = get_until_correct(eprompt,Float64)
    #     eunit = uparse(eunit)
    #     energy = eval*eunit

    #     println("energy = ", energy)

    #     # position

    #     # direction
    # end
    # println("Number of particles = ", nparticles)
    # #    # energies at which to generate particles
    # #    energy = [0.01,0.1,1.0,10.0]*u"TeV"
    # #    # fire particles from different galactic radii
    # #    rStart = 8.5u"kpc"
    # #    # directions to fire particles (these will be normalized)
    # #    outRadially  = [-1.0,0.0,0.0]
    # #    particleAim = Vector([outRadially])
    # #
    # #    gpFrac = 0.05
    # #
    # #    filename = ""
    # #    for en in energy
    # #        # Generate particle tree
    # #        pnsPT = singleParticleTree(proton,en,rStart,outRadially,gpFrac)
    # #
    # #        # Plot trajectory
    # #        pnsScene!(scene,pnsPT,filename)
    # #    end

    # # return scene
    # return particle

end

export run

end # PointNshoot

import .PointNshoot
PointNshoot.run()