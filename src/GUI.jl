
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
    if length(scene.plots) > 2
        delete!(scene,scene.plots[end])		
        clear!(scene)
    elseif length(scene.plots) == 2
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
