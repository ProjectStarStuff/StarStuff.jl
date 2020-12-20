module ParticleTreeTest
	using DrWatson

    using Unitful
    using UnitfulAstro
    using UnitfulAtomic
    using LinearAlgebra
    using Random
    using DifferentialEquations
	using MultiScaleArrays
	using RecursiveArrayTools
	using Test

	# include(srcdir("Particle.jl"))
	# include(srcdir("Generator.jl"))
	# include(srcdir("ParticleTree.jl"))
	include(srcdir("CrProp.jl"))
	using .CrProp

	function test()

		function randSpotInGalaxy(rng = MersenneTwister())
			r  = 10.0
			ϕ  = 2*π*rand(rng)
			z  = randbtw(-1.0,1.0,rng)
			return [r*cos(ϕ),r*sin(ϕ),z]
		end

		# generate electrons
		eStart = 1.0
		eStop  = 10.0
		eα     = 3.0

		powerLawElectrons = PowerLaw(eα,(eStart,eStop))

		electron = ParticleID(0,0,1)

		nElectrons = 10
		electronCoords = Vector([Coords(ArrayPartition(randSpotInGalaxy(), rand(powerLawElectrons)*unitRandom3Cartesian())) for i=1:nElectrons])
		testElectrons = construct(Species,deepcopy(electronCoords),Float64[],electron)

		# generate protons
		pStart = 1.0
		pStop  = 10.0
		pα     = 2.7

		powerLawProtons = PowerLaw(pα,(pStart,pStop))

		proton = ParticleID(1,0,0)

		nProtons = 5
		protonCoords = Vector([Coords(ArrayPartition(randSpotInGalaxy(), rand(powerLawProtons)*unitRandom3Cartesian())) for i=1:nProtons]) 
		testProtons = construct(Species,deepcopy(protonCoords),Float64[],proton)

		# make snapshot(s)
		snap = construct(Snapshot,[deepcopy(testElectrons),deepcopy(testProtons)])
		times = [0.0,1.0]

		# make a timescale
		testTimescale = construct(Timescale,[deepcopy(snap),deepcopy(snap)],Float64[],times)

		# make MultiParticleHistory
		testPT = construct(ParticleTree,[deepcopy(testTimescale),deepcopy(testTimescale)])


		print_human_readable(testPT)
		println()
		println(length(testPT))
		println(testPT[1,1,1,1].values.x[2])
		

	end #test

	export test
end
