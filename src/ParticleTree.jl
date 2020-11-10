# TODO: Add empty tree constructors
# BEGIN PARTICLE TREE
"""
Leaf structure in galaxy tree consists of position and momentum coordinates of a specific particle
"""
struct Coords{B} <: AbstractMultiScaleArrayLeaf{B}
    values::ArrayPartition{B}

    function Coords(values::ArrayPartition{B}) where B <: Number
        @assert length(values.x[1]) == 3 "Position must have length == 3"
        @assert length(values.x[2]) == 3 "Momentum must have length == 3"
        new{B}(values)
    end

end # Coords

"""
Internally consists of a particle ID and an array of particle coordinates
"""
struct Species{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
    id::ParticleID
end # Species

"""
Internally consists of a collection of all species at one particular time
"""
struct Snapshot{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

"""
Internally a collection of all Snapshots with the same time stepping
"""
struct Timescale{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
    times::Vector{Float64}
    dt::Float64

    function Timescale(nodes::Vector{T},values::Vector{B},end_idxs::Vector{Int},times::Vector{Float64},dt::Float64 = -1.0) where {T,B}
        @assert length(nodes) == length(times) "Number of nodes must equal the number of times"
        new{T,B}(nodes,values,end_idxs,times,dt)
    end
end

"""
Collection of all timescales
"""
struct ParticleTree{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end	

export Coords, Species, Snapshot, Timescale, ParticleTree
# END PARTICLE TREE

###############################################################################
# FUNCTIONS TO HELP GENERATE PARTICLE TREES

"""
    single_particletree(particletype::T, energy::U, xyz::V, aim::W, dt::X) where {T<:ParticleID,U<:typeof(1.0u"J"),V<:typeof(Vector([1.0,1.0,1.0]*u"m")),W<:typeof(Vector([1.0,0.0,0.0])),X<:typeof(1.0u"s")}

Return a ParticleTree with a single particle in it.

Input:
- particletype: ParticleID
- energy: initial particle energies (must have units ∈ Unitful.Energy)
- xyz: starting Cartesian coordinates with Earth along the -x axis (must have units ∈ Unitful.Distance)
- aim: Vector with length=3 with initial direction of particle trajectory. 
- dt: initial time step for simulation
- bfield: Function that returns the value of the galactic magnetic field when given spatial coorinates

Returns:
        ParticleTree
"""
function single_particletree(particletype::T, energy::U, xyz::V, aim::W, dt::X) where {T<:ParticleID,U<:typeof(1.0u"J"),V<:typeof(Vector([1.0,1.0,1.0]*u"m")),W<:typeof(Vector([1.0,0.0,0.0])),X<:typeof(1.0u"s")}

    # Get ParticleTree leaf type
    leafType = typeof(Coords(ArrayPartition([1.,2.,3.],[4.,5.,6.])))

    # Ensure that the initial direction is normalized for subsequent
    # momentum scaling.
    particleAim = map(x->normalize(x),aim_vec)

    # starts at time 0.0
    timeInit = Vector([0.0])

    particleTree = nothing

    # first make Coords
    pos = usistrip.(xyz)
    dir = normalize(aim)
    pmag = usistrip(momentum(energy,particletype))
    mom = pmag*dir
    coords = Coords(ArrayPartition(pos,mom))

    # stick Coords in a Species
    sp = construct(Species,Vector([coords]),Float64[],particletype)

    # stick Species in a Snapshot
    snap = construct(Snapshot,Vector([sp]))

    # stick Snapshot in a Timescale
    time_init = Vector([0.0])
    dt_strip = usistrip(dt)
    scale = construct(Timescale,Vector([snap]),Float64[],time_init,dt_strip)

    # stick Timescale in a ParticleTree
    particletree = construct(ParticleTree,Vector([scale]))

    return particletree

end # single_particletree
export single_particletree
"""
    initParticles(particleTypes, energies, rStart, particleAim, gpFrac, bfield)

Return a ParticleTree consisting of all combinations of particles, energies, radii, and initial direction.
Timestep is calculated as a fraction of the gyroperiod at the particle point of insertion.

Input:
- particleTypes: Vector of ParticleID(s)
- energies: Vector of initial particle energies (must have units ∈ Unitful.Energy)
- rStart: starting radius from galactic center (must have units ∈ Unitful.Distance)
- particleAim: Vector with length=3 with initial direction of particle trajectory. 
- gpFrac: Float giving the fraction of the gyroperiod to use as a time step
- bfield: Function that returns the value of the galactic magnetic field when given spatial coorinates

Returns:
        ParticleTree
"""
function initParticles(particleTypes, energies, rStart, particleAim, gpFrac, bfield)

    # Get ParticleTree leaf type
    leafType = typeof(Coords(ArrayPartition([1.,2.,3.],[4.,5.,6.])))

    # Ensure that the initial direction is normalized for subsequent
    # momentum scaling.
    particleAim = map(x->normalize(x),particleAim)

    # starts at time 0.0
    timeInit = Vector([0.0])

    particleTree = nothing
    # loop over all arrays to get combination of all parameters 
    for r in rStart
        rVec = Vector([usistrip(r),0.0,0.0])	
        bMag = norm(bfield(rVec))	
        for en in energies
            for p in particleTypes
                pMag = momentum(en,p) 
                vMag = vFromP(pMag,mass(p)) 
                coordsVec = nothing
                for d in particleAim
                    pVec = usistrip(pMag)*d 
                    coords = Coords(ArrayPartition(rVec,pVec))	
                    if coordsVec == nothing
                        coordsVec = Vector([coords]) 
                    else
                        push!(coordsVec,coords) 
                    end
                end
                # each species will have a different gyroradius, so each
                # one also gets its own time step -> Timescale
                sp = construct(Species,coordsVec,Float64[],p) 
                spVec = Vector([sp])
                snap = construct(Snapshot,spVec)
                snapVec = Vector([snap])
                gCirc = 2*π*gyroradius(
                                       ustrip(u"GeV/c",pMag),
                                       abs(p.p-p.e),
                                       usistrip(bMag)*1.0e4
                                       )
                gPrd = gCirc/usistrip(vMag)
                dt = gPrd*gpFrac
                scale = construct(Timescale,snapVec,Float64[],timeInit,dt)
                if particleTree == nothing
                    particleTree = construct(ParticleTree,Vector([scale]))
                else
                    add_node!(particleTree,scale)
                end
            end
        end
    end
    return particleTree

end # initParticles

# TODO: generate readable lookup tree if inputs are given as dictionaries

# TODO: generate from a (.json?) file

###############################################################################
# FUNCTIONS TO HELP READ PARTICLE TREES

"""
        getCoordArray(pt::ParticleTree)

Returns all (position,momentum) coordinates in an array of size (6,nCoords)
"""
function getCoordArray(pt::ParticleTree)
    out = Vector(pt)
    out = transpose(reshape(pt,(6,length(pt)÷6)))
    return out
end

export getCoordArray

"""
        getTimescaleIndices(pt::ParticleTree)

Returns internal MultiScaleArray indices denoting the last index of each 
Timescale.
"""
function getTimescaleIndices(pt::ParticleTree)
    return pt.end_idxs
end

export getTimescaleIndices

"""
        getSnapshotIndices(pt::ParticleTree)

Returns internal MultiScaleArray indices denoting the last index of each
Snapshot.
"""
function getSnapshotIndices(pt::ParticleTree)
    out = nothing
    for snapshot in LevelIter(2,pt)
        if out == nothing
            out = Vector([snapshot[1].end_idxs])
        else
            push!(out,snapshot[1].end_idxs)
        end
    end
    return out
end

export getSnapshotIndices

"""
        getSpeciesIndices(pt::ParticleTree)

Returns internal MultiScaleArray indices denoting the last index of each
Species.
"""
function getSpeciesIndices(pt::ParticleTree)
    out = nothing
    for species in LevelIter(3,pt)
        if out == nothing
            out = Vector([species[1].end_idxs])
        else
            push!(out,species[1].end_idxs)
        end
    end
    return out
end

export getSpeciesIndices

"""
        getParticleTrajectory(pt::ParticleTree,iTs::Int64,iSs::UnitRange,iSp::Int64,iCd::Int64)

Returns particle trajectory over given range iSs
"""
function getParticleTrajectory(pt::ParticleTree,iTs::Int64,iSs::UnitRange,iSp::Int64,iCd::Int64)
    return pt.nodes[iTs].times, Vector([pt[iTs,i,iSp,iCd].values for i = iSs])
end

"""
        getParticleTrajectory(pt::ParticleTree,iTs::Int64,iSs::UnitRange,iSp::Int64,iCd::Int64)

Returns entire particle trajectory
"""
function getParticleTrajectory(pt::ParticleTree,iTs::Int64,iSp::Int64,iCd::Int64)
    nTime = length(pt.nodes[iTs].nodes)
    getParticleTrajectory(pt,iTs,1:nTime,iSp,iCd)
end

export getParticleTrajectory
###############################################################################
# SAVE PARTICLE TREES TO FILE
