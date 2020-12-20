using DrWatson
@quickactivate "crprop"


module SpeciesTest
    using DrWatson
    using Test

    include(srcdir("CrProp.jl"))
    using .CrProp

    """
        test()

    Run tests on cosmic ray species structs
    """
    function test()
        electron = Cr(0,0,1)
        electron2 = Cr(0,0,1)

        @testset "CrSpecies tests" begin
            @test qcr(electron) == -1*Unitful.q
            @test mcr(electron) == upreferred(511.0u"keV/c^2")
            @test electron == electron2
        end
    end
end #SpeciesTest

SpeciesTest.test()
