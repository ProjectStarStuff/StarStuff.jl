
module PointNshoot
include("StarStuff.jl")

using .StarStuff


using Blink
using Interact

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
    update!(scene)

    return scene
end

function sineplot!(a::Int64,s::Scene)
    cam2d!(s)
    t = collect(0:0.1:10)
    # x = cos.(2*π*t)
    y = sin.(2*π*t)
    # z = t
    println("Sin plot")
    lines!(s,t,y)
    update!(s)
end

function cosplot!(a::Int64,s::Scene)
    cam2d!(s)
    t = collect(0:0.1:10)
    x = cos.(2*π*t)
    # y = sin.(2*π*t)
    # z = t
    println("Cos plot")
    lines!(s,t,x)
    update!(s)
end

function spiplot!(a::Int64,s::Scene)
    cam3d!(s)
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

"""
    function particlepicker()

Widget to select a particle. Returns a ParticleID
"""
function particlepicker()
    el = spinbox(label="Electrons"; value=0)
    pr = spinbox(label="Protons"; value=0)
    nu = spinbox(label="Neutrons"; value=0)
    output = Interact.@map(ParticleID(&pr,&nu,&el))
    wdg = Widget(["e"=>el,"p"=>pr,"n"=>nu],output=output)
    @layout! wdg Interact.hbox(pr,nu,el)
end

"""
        run()

Run PointNshoot script
"""
function run()
    GLMakie.activate!()

    s = Observable{Any}(Scene()) # scene output
    pt = Observable{Any}(nothing) # stored particles
    pid = Observable{Any}(nothing)

    # sine_button = Interact.button("Sine plot")
    # cos_button = Interact.button("Cos plot")
    # spi_button = Interact.button("Spiral plot")
    # clear_button = Interact.button("Clear plots")

    # map!(sineplot!,s,sine_button,s[])
    # map!(cosplot!,s,cos_button,s[])
    # map!(spiplot!,s,spi_button,s[])
    # map!(clearplot!,s,clear_button,s[])

    # *** PARTICLE ID INPUT ***
    part_text = Interact.latex("\\text{\\textbf{Particle components}}")
    el = Interact.spinbox(label="Electrons"; value=0,step="any")
    pr = Interact.spinbox(label="Protons"; value=0)
    nu = Interact.spinbox(label="Neutrons"; value=0)
    part_components = Interact.hbox(pr,nu,el)

    part_display = Interact.vbox(part_text,part_components)

    # *** POSITION INPUT ***
    dist_names = ["kpc","pc","km","m","cm"]
    dist_objects = [1.0u"kpc",1.0u"pc",1.0u"km",1.0u"m",1.0u"cm",]
    dist_dict = OrderedDict(zip(dist_names,dist_objects))
    dist_wdg = Interact.dropdown(dist_names)
    # Interact.@on println(&unit_wdg)
    dist_disp = Interact.@map Interact.latex("\\text{$(&dist_wdg)}")

    xstart = Interact.spinbox(label="X: "; value = -8.5)
    xval = Interact.@map &xstart*dist_dict[&dist_wdg]
    xdisp = Interact.hbox(xstart)
    ystart = Interact.spinbox(label="Y: "; value = 0.0)
    yval = Interact.@map &ystart*dist_dict[&dist_wdg]
    ydisp = Interact.hbox(ystart)
    zstart = Interact.spinbox(label="Z: "; value = 0.0)
    zval = Interact.@map &zstart*dist_dict[&dist_wdg]
    zdisp = Interact.hbox(zstart)

    dist_text = Interact.hbox(latex("\\text{\\textbf{Initial Position}(}"),dist_wdg,latex("\\text{)}"))
    dist_display = Interact.vbox(dist_text,xdisp,ydisp,zdisp)

    # *** ENERGY INPUT ***
    enu_names = ["GeV","TeV","MeV","keV"]
    enu_objects = [1.0u"GeV",1.0u"TeV",1.0u"MeV",1.0u"keV"]
    enu_dict = OrderedDict(zip(enu_names,enu_objects))
    enu_wdg = Interact.dropdown(enu_names)
    enu_disp = Interact.@map Interact.latex("\\text{$(&enu_wdg)}")

    enstart = Interact.spinbox(value = 0.0)
    enval = Interact.@map &enstart*enu_dict[&enu_wdg]
    en_disp = Interact.hbox(enstart,enu_disp)

    energy_display = Interact.hbox(latex("\\text{\\textbf{Energy = }}"),enstart,enu_wdg)

    # *** DIRECTION INPUT ***
    aim_text = Interact.latex("\\text{\\textbf{Initial Direction}}")
    xaim = Interact.spinbox(label="X: ",value = 1.0)
    yaim = Interact.spinbox(label="Y: ",value = 0.0)
    zaim = Interact.spinbox(label="Z: ",value = 0.0)
    aim_val = Interact.@map normalize(Vector([&xaim,&yaim,&zaim]))

    aim_display = Interact.vbox(aim_text,xaim,yaim,zaim)

    pos_aim = Interact.hbox(dist_display,aim_display)
 
    # *** TIME INPUT***

    ## Manual timestep
    # make dictionary for time and dt units
    timeu_names   = ["sec","min","hr","day","yr","kyr","Myr"]
    timeu_objects = [1.0u"s",1.0u"minute",1.0u"hr",1.0u"d",1.0u"yr",1.0u"kyr",1.0u"Myr"]
    timeu_dict    = OrderedDict(zip(timeu_names,timeu_objects))

    # make time and dt unit widgets
    timeu_wdg = Interact.dropdown(timeu_names)
    dtu_wdg   = Interact.dropdown(timeu_names)

    # make numerical input widgets for time and dt
    time_num = Interact.spinbox(label = "Duration:",value = 0.0)
    dt_num   = Interact.spinbox(label = "Time step:",value = 0.0)

    # bind time and dt widget values into quantities
    time_val = Interact.@map &time_num*timeu_dict[&timeu_wdg]
    dt_val   = Interact.@map &dt_num*timeu_dict[&dtu_wdg]

    # join number and unit widgets together for display
    time_disp = Interact.hbox(time_num,timeu_wdg)
    dt_disp = Interact.hbox(dt_num,dtu_wdg)

    # join time and dt widgets together
    mantime_disp = Interact.vbox(time_disp,dt_disp)

    ## Automatic timestep
    # get the field value at the current location

    autotime_disp = Interact.latex("\\text{Not yet implemented...}")


    ## Choose between time modes
    # make dictionary of different time modes
    timemodes_list = ["MANUAL","AUTOMATIC"]
    timemodes_objects = [mantime_disp,autotime_disp]
    timemodes_dict = OrderedDict(zip(timemodes_list,timemodes_objects))

    # make toggle widget for different time modes
    timemodes_wdg = Widgets.togglebuttons(timemodes_list)

    # make observable for proper time display
    timemodes_disp = Interact.@map timemodes_dict[&timemodes_wdg]

    # join toggle widget and time input
    time_text = Interact.latex("\\text{\\textbf{Time settings}}")
    time_section = Interact.@map Interact.vbox(time_text,timemodes_wdg,&timemodes_disp)

    ## *** Construct particle
    pt_obj = Observable(ParticleTree())
    clearpt_button = Interact.button("Clear particles")
    addpt_button   = Interact.button("Add particle")
    pt_buttons = Interact.hbox(pad(1em,addpt_button),pad(1em,clearpt_button))

    # TODO: add particle whenever "addpt_button" is pressed

    # *** MAKE DISPLAY ***
    ui = dom"div"(
        part_display,
        Interact.hline(),
        energy_display,
        Interact.hline(),
        pos_aim,
        Interact.hline(),
        time_section,
        Interact.hline(),
        pt_buttons
        )
    w  = Window()
    body!(w,ui)

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