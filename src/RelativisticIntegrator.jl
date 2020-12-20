"""
    function stepBoris(x,p,efield,bfield,q,m,dt)

Performs one complete time step of the Boris algorithm for relativistic charged particles in an electromagnetic field. This version is from Ripperda2017
"""
function stepBoris(x, p, efield, bfield, q, m, dt)
    qdt2m = q * dt / (2.0 * m)

    un = p / m
    γn = γFromP(p, m)
    xnph = x + un * dt / (2.0 * γn)

    eph = efield(xnph)
    um = un + qdt2m * eph

    γm = γFromP(m * um, m)
    t = (qdt2m / γm) * bfield(xnph)
    s = (2.0 / (1 + t ⋅ t)) * t

    up = um + cross((um + cross(um, t)), s)

    unp1 = up + eph * q * dt / (2.0 * m)
    pnew = m * unp1

    γnp1 = γFromP(pnew, m)
    xnew = xnph + (dt / (2.0 * γnp1)) * unp1
    return ArrayPartition(xnew, pnew)
end
export stepBoris


"""
    fBoris!(u::ArrayPartition, params::Dict, t)

Boris stepper wrapper for the discrete solver.
Use with single particle coordinates.

Input:
    - u{ArrayPartition{Number}} ≡ ArrayPartition([position],[momentum])
    - params{Dict} ≡ Dict[
            :t => {Number} current time value,
            :efield => {Function} returns electric field at [x,y,z]
            :bfield => {Function} returns magnetic flux density at [x,y,z]
            :q      => {Number} charge of particle in SI
            :m      => {Number} mass of particle in SI
        ]
    - t{Number} ≡ time at next time step

Returns:
    Particle coordinates at time t in the form 
    {ArrayPartition{Number}}([[x,y,z],[px,py,pz]]).
"""
function fBoris!(u::ArrayPartition, params::Dict, t)
    dt = t - params[:t]
    params[:t] = t
    stepBoris(
        u.x[1],
        u.x[2],
        params[:efield],
        params[:bfield],
        params[:q],
        params[:m],
        dt
    )
end

"""
    fBoris!(u::Snapshot, params::Dict, t)

Boris stepper wrapper for the discrete solver.
Use with Snapshot branch of ParticleTree.

Input:
    - u{Snapshot}  ≡ level 2 ParticleTree branch
    - params{Dict} ≡ Dict[
            :t      => {Number} current time value,
            :efield => {Function} returns electric field at [x,y,z]
            :bfield => {Function} returns magnetic flux density at [x,y,z]
        ]
    - t{Number}    ≡ time at next time step

Returns:
    Snapshot at time t
"""
function fBoris!(u::Snapshot, params::Dict, t)
    dt = t - params[:t]
    params[:t] = t
    snapOut = deepcopy(u)
    for species in snapOut.nodes
        q = usistrip(charge(species.id))
        m = usistrip(mass(species.id))
        map(species.nodes) do coords
            stepCoords = stepBoris(
                    coords.values.x[1],
                    coords.values.x[2],
                    params[:efield],
                    params[:bfield],
                    q,
                    m,
                    dt
                )
            coords.values .= stepCoords
        end
    end
    return snapOut
end

export fBoris!


# dictionary of all relativistic integrators availble
rintegrator = Dict(
                :boris => fBoris!
)

export rintegrator
