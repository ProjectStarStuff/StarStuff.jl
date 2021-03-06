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
    @reexport using CSV
    @reexport using Statistics
    @reexport using DataFrames
    @reexport using DataStructures

    @reexport using Makie
    @reexport using AbstractPlotting
    @reexport using AbstractPlotting.MakieLayout

    # @reexport using Interact
    # @reexport using Blink

    include("Particle.jl")
    include("Kinematics.jl")    
    include("Generator.jl")
    include("PropUtils.jl")
    include("ParticleTree.jl")
    include("Bfield.jl")
    include("RelativisticIntegrator.jl")
    include("GUI.jl")
    
    # scripts
    # include("PointNshoot.jl")
end # StarStuff
