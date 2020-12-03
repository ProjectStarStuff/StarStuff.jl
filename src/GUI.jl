
# DEFINE WIDGETS
"""
particle_wdg()

Widget to get particle type. Returns Observable{ParticleID}
"""
function particle_wdg()
    part_text = Interact.latex("\\text{\\textbf{Particle components}}")
    el = Interact.spinbox(label="Electrons"; value=1)
    pr = Interact.spinbox(label="Protons"; value=0)
    nu = Interact.spinbox(label="Neutrons"; value=0)

    output = Interact.@map ParticleID(&pr,&nu,&el)

    wdg = Widget(["p"=>pr,"n"=>nu,"e"=>el],output = output)
    part_components = Interact.hbox(pr,nu,el)
    Interact.@layout! wdg Interact.vbox(part_text,part_components)
end
export particle_wdg

"""
position_wdg()

Widget to get position. Returns Observable{Array{Quantity{Length},1}}   
"""
function position_wdg()
    # *** POSITION INPUT ***
    dist_names = ["kpc","pc","km","m","cm"]
    dist_objects = [1.0u"kpc",1.0u"pc",1.0u"km",1.0u"m",1.0u"cm",]
    # use an OrderedDict to display desired units first
    dist_dict = OrderedDict(zip(dist_names,dist_objects))
    dist_wdg = Interact.dropdown(dist_names)
    # Interact.@on println(&unit_wdg)

    xstart = Interact.spinbox(label="X: "; value = -8.5)
    xval = Interact.@map &xstart*dist_dict[&dist_wdg]
    xdisp = Interact.hbox(xstart)
    ystart = Interact.spinbox(label="Y: "; value = 0.0)
    yval = Interact.@map &ystart*dist_dict[&dist_wdg]
    ydisp = Interact.hbox(ystart)
    zstart = Interact.spinbox(label="Z: "; value = 0.0)
    zval = Interact.@map &zstart*dist_dict[&dist_wdg]
    zdisp = Interact.hbox(zstart)

    output = Interact.@map Vector([&xval,&yval,&zval])

    wdg = Widget(["x"=>xstart,"y"=>ystart,"z"=>zstart,"u"=>dist_wdg],output=output)

    dist_text = Interact.latex("\\text{\\textbf{Initial Position}}")
    dist_utext = Interact.hbox(Interact.latex("\\text{Units: }"),dist_wdg)
    Interact.@layout! wdg Interact.vbox(dist_text,dist_utext,xdisp,ydisp,zdisp)
end
export position_wdg

"""
energy_wdg()

Widget to get energy. Returns Observable{Quantity{Energy}}
"""
function energy_wdg()
    # unit input
    enu_names = ["GeV","TeV","MeV","keV"]
    enu_objects = [1.0u"GeV",1.0u"TeV",1.0u"MeV",1.0u"keV"]
    enu_dict = OrderedDict(zip(enu_names,enu_objects))
    enu_wdg = Interact.dropdown(enu_names)
    enu_disp = Interact.@map Interact.latex("\\text{$(&enu_wdg)}")
    # value input
    enstart = Interact.spinbox(value = 1.0)
    # combine unit and value
    output = Interact.@map &enstart*enu_dict[&enu_wdg]
    # define widget
    wdg = Widget(["e"=>enstart,"u"=>enu_wdg],output=output)
    # return widget layout
    Interact.@layout! wdg Interact.hbox(latex("\\text{\\textbf{Energy = }}"),enstart,enu_wdg)
end
export energy_wdg

"""
    direction_wdg()

Widget to get initial particle direction. Returns Observable{Array{Float64,1}}
"""
function direction_wdg()
    text = Interact.latex("\\text{\\textbf{Initial Direction}}")
    x = Interact.spinbox(label="X: ",value = 1.0)
    y = Interact.spinbox(label="Y: ",value = 0.0)
    z = Interact.spinbox(label="Z: ",value = 0.0)
    output = Interact.@map Vector([&x,&y,&z])

    wdg = Widget(["x"=>x,"y"=>y,"z"=>z], output = output)

    Interact.@layout! wdg Interact.vbox(text,x,y,z)    
end
export direction_wdg

"""
    manualdt_wdg()

Widget to manually set particle timestep. Returns Observable{Quantity{Time}}
"""
function manualdt_wdg()
    # make dictionary for time and dt units
    tu_names   = ["sec","min","hr","day","yr","kyr","Myr"]
    tu_objects = [1.0u"s",1.0u"minute",1.0u"hr",1.0u"d",1.0u"yr",1.0u"kyr",1.0u"Myr"]
    tu_dict    = OrderedDict(zip(tu_names,tu_objects))

    # make time and dt unit widgets
    dtu = Interact.dropdown(tu_names)

    # make numerical input widgets for time and dt
    dt = Interact.spinbox(value = 1.0)

    # bind units and values into quantities
    output = Interact.@map &dt*tu_dict[&dtu]

    # construct widget
    wdg = Widget(["dt"=>dt,"u"=>dtu], output = output)

    # join number and unit widgets together for display
    dt_text = Interact.latex("\\text{\\textbf{Time Step}}")
    dt_val  = Interact.hbox(dt,dtu)
    Interact.@layout! wdg Interact.vbox(dt_text,dt_val)
end
export manualdt_wdg

"""
    pt_manualdt_wdg()

Widget to generate particles manual timestep. Returns Observable{ParticleTree}
"""
function pt_manualdt_wdg()
    particle  = particle_wdg()
    energy    = energy_wdg()
    position  = position_wdg()
    direction = direction_wdg()
    dt        = manualdt_wdg()

    output = Observable{Any}(ParticleTree())
    # ptnew  = Observable(ParticleTree())

    button_addpt = Interact.button("Add particle")

    ptnew = Interact.@map (&button_addpt; ParticleTree(particle[],energy[],position[],direction[],dt[]))
    # output = Interact.@map (&button_addpt; ParticleTree(particle[],energy[],position[],direction[],dt[]))

    # Interact.@on (&button_addpt;println("inside function: output length = ", length(getproperty(output[],:nodes))))

    Interact.@on push!(output[],&ptnew)
    connect!(ptnew,output)

    wdg = Widget(["particle"=>particle,"energy"=>energy,"position"=>position,"direction"=>direction,"dt"=>dt], output = output)

    pos_dir = Interact.hbox(position,direction)

    Interact.@layout! wdg Interact.vbox(particle,
                                        hline(),
                                        energy,
                                        hline(),
                                        pos_dir,
                                        hline(),
                                        dt,
                                        hline(),
                                        button_addpt
                                        )
end
export pt_manualdt_wdg
# END DEFINE WIDGETS
###############################################################################

"""
    get_until_correct(prompt::String, dtype::Type)

Prompts user until a valid input is given
"""
function get_until_correct(prompt::String, dtype::Type; parser::Function = parse)
    out = nothing
    try
        print(prompt)
        out = readline()
        out = parse(dtype,out)
    catch
        println("Invalid entry. Try again...")
        out = get_until_correct(prompt,dtype)
    end
    return out
end
export get_until_correct

function parse_array(a::String, dtype::Type; token::String = ",")
    out_str = strip.(split(a,token))
    out = parse.(dtype,out_str)
    return out
end

"""
    request_particle_type()

Prompts user to enter a particle type. Returns a ParticleID.

"""
function request_particleid()
    particle_menu = request("Select a particle:",RadioMenu(particle_keys))
    particle = particle_dict[particle_keys[particle_menu]]
    return particle
end


"""
    request_energy()

Prompts user to enter energy. Returns Quantity{Float64,Energy}
"""
function request_energy()
    eunits = ["MeV","GeV","TeV"]
    eunits_menu = request("Select energy units:",RadioMenu(eunits))
    eunit = eunits[eunits_menu]
    eprompt = "Enter energy value ($units):"
    eval = get_until_correct(eprompt,Float64)
    eunit = uparse(eunit)
    energy = eval*eunit
    return energy
end


"""
    request_position()

Prompts user to enter a position. Returns Array{Quantity{Float64,Length},1}
"""
function request_position()
    posprompt = "Enter initial position in kpc as `XX.XX, YY.YY, ZZ.ZZ`: "
    pos = get_until_correct(posprompt,Float64,parse_array)
    return pos
end


"""
    request_direction()

Prompts user to enter a direction. Returns Array{Float64,1}
"""
function request_direction()
    dirprompt = "Enter initial direction as `XX.XX, YY.YY, ZZ.ZZ`: "
    dir = get_until_correct(dirprompt,Float64,parse_array)
    return dir
end


"""
    request_particle_ic()

Prompts user to enter initial conditions for a single particle
"""
function request_particle_ic()
    particle = request_particleid()
    energy   = request_energy()
    position = request_position()
    direction = request_direction()
    dtmenu = ["Automatic", "Manual"]
    dtselect = request("Time stepping: ",RadioMenu(dtmenu)) 
    dt = nothing
    if dtselect == 1
        dtprompt = "Enter dt (s): "
        dt       = get_until_correct(dtprompt,Float64)    
    else
        timeprompt = "Enter duration (s): "
        time       = get_until_correct(timeprompt, Float64)
    end

end


"""
    particle_prompt()

Gets input from user to construct initial particle tree
"""
function particle_prompt()


end

"""
    clear!(scene::Scene)

Removes everything from scene except axes
"""
function clear!(scene::Scene) 
    if length(scene.plots) > 1
        delete!(scene,scene.plots[end])		
        clear!(scene)
    elseif length(scene.plots) == 1
        delete!(scene,scene.plots[end])		
    end
end

export clear!

"""
    scaleShift(a::T) where {T<:AbstractVector}

Shifts values of array so it is centered on its mean and range is normalized to 1
"""
function scaleshift(a::T) where {T<:AbstractVector}
    aRange = max(a...) - min(a...)
    aMean  = mean(a)
    out = a .- aMean
    if aRange == 0.0
        return out
    else
        return aMean, aRange, out ./ aRange
    end
end

export scaleshift
