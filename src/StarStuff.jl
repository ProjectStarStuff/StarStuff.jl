module StarStuff
    using Reexport

    # external dependencies
    @reexport using Printf
    @reexport using Unitful
    @reexport using UnitfulAstro
    @reexport using UnitfulAtomic
    @reexport using LinearAlgebra
    @reexport using Random
    @reexport using Distributions
    @reexport using DifferentialEquations
    @reexport using MultiScaleArrays
    @reexport using GeometryBasics
    @reexport using IterTools
    @reexport using Measurements

    @reexport using Makie
    @reexport using AbstractPlotting
    @reexport using AbstractPlotting.MakieLayout

    include("Particle.jl")
    include("Kinematics.jl")    
    include("Generator.jl")
    include("PropUtils.jl")
    include("ParticleTree.jl")
    include("Bfield.jl")
    include("RelativisticIntegrator.jl")
    include("GUI.jl")
end # StarStuff
