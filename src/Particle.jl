"""
Associated with various objects that contain constant particle properties
"""
abstract type Particle end


"""
ID cosmic ray by number of protons, neutrons, and electrons
"""
struct ParticleID <: Particle
    p::Int64
    n::Int64
    e::Int64
end # ParticleID

fundamentalRestMassDict = Dict(
    :e => upreferred(511.0u"keV/c^2"),
    :p => upreferred(938.0u"MeV/c^2")
    )

"""
    mass(r::ParticleID; m::Dict = fundamentalRestMass)

Simple particle mass.

TODO: Account for binding energy
"""
mass(r::ParticleID; m::Dict = fundamentalRestMassDict) = r.p*m[:p]+r.e*m[:e] + 1.008665*m[:p]*r.n

"""
    charge(r::ParticleID)

Returns particle charge
"""
charge(r::ParticleID) = (r.p-r.e)*Unitful.q

Base.:(==)(a::ParticleID, b::ParticleID) = (a.p == b.p && a.n == b.n && a.e == b.e) ? true : false
Base.length(a::ParticleID) = 1 #defined for broadcasting purposes
Base.iterate(a::ParticleID) = (a,nothing)
export ParticleID, mass, charge
