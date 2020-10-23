using DrWatson
@quickactivate "crprop"


module KinematicsTest
    using DrWatson
    using Test

    include(srcdir("CrProp.jl"))
    using .CrProp

    """
        test()

    Run tests on kinematics functions
    """
    function test(tol = 1e-9)
        residual(val,expect) = abs((expect-val)/expect)

        electron = Cr(0,0,1)
        eElectron = 1.0u"GeV"
        mElectron = uconvert(u"GeV/c^2", mcr(electron))
        pElectron = momentum(eElectron, mElectron)

        eElectronSI = usistrip(eElectron)
        mElectronSI = usistrip(mElectron)
        pElectronSI = usistrip(pElectron)
        cSI         = float(Unitful.c0.val)
        @testset "Cosmic ray kinematics tests" begin
            @test β(1) == 0
            @test residual(eElectron^2-pElectron^2*Unitful.c^2, mElectron^2*Unitful.c^4) < tol
            @test residual(eElectronSI^2-pElectronSI^2*cSI^2, mElectronSI^2*cSI^4) < tol
            #TODO: Complete tests for all functions in Kinematics.jl

        end
    end
end #KinematicsTest

KinematicsTest.test()
