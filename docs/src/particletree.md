# Particle Trees
StarStuff uses the [MultiScaleArrays.jl](https://github.com/SciML/MultiScaleArrays.jl) package to keep track of particle properties by constructing a `ParticleTree`.

## `ParticleTree` constructors

```@docs
Coords{B}
Species{T<:AbstractMultiScaleArray,B<:Number}
Snapshot{T<:AbstractMultiScaleArray,B<:Number}
Timescale{T<:AbstractMultiScaleArray,B<:Number}
ParticleTree{T<:AbstractMultiScaleArray,B<:Number}
```

## `ParticleTree` generating functions

```@docs
ParticleTree(particletype::ParticleID,energy::Quantity,xyz::Vector{T},aim::Vector,dt::Quantity) where {T<:Quantity}
ParticleTree(particletype::ParticleID, energy::Float64, xyz::Vector{Float64}, aim::Vector{Float64}, dt::Float64)
push!(pt::ParticleTree,val::Union{ParticleTree,Timescale,Snapshot,Species,Coords})
push(pt::ParticleTree,val::Union{ParticleTree,Timescale,Snapshot,Species,Coords})
initParticles(particleTypes, energies, rStart, particleAim, gpFrac, bfield)
```

## Functions to get properties from a `ParticleTree` 

```@docs
getCoordArray(pt::ParticleTree)
getTimescaleIndices(pt::ParticleTree)
getSnapshotIndices(pt::ParticleTree)
getSpeciesIndices(pt::ParticleTree)
getParticleTrajectory(pt::ParticleTree,iTs::Int64,iSs::UnitRange,iSp::Int64,iCd::Int64)
getParticleTrajectory(pt::ParticleTree,iTs::Int64,iSs::UnitRange,iSp::Int64,iCd::Int64)
```