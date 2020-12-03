module TestStuff

include("StarStuff.jl")
using .StarStuff


function run()

    pt = pt_manualdt_wdg()

    # Interact.@on println("Number of particles = ",length(getproperty(&pt,:nodes)))
    Interact.@on println("*** New Particle ***\n",&pt)

    w = Window()
    ui = dom"div"(
        pt
    )
    body!(w,ui)


    # test = ParticleTree()
    # println(length(test.nodes))
    # push!(test,ParticleTree())
    # println("length(test.nodes) = ", length(test.nodes))
    # # parameters for push! test
    # el = ParticleID(0,0,1)
    # en = upreferred(1.0u"GeV")
    # pos = Vector([-8.5,0.0,0.0]*u"kpc")
    # dir = Vector([1.0,0.0,0.0])
    # dt = 1.0u"s"

    # elpt = ParticleTree(el,en,pos,dir,dt)
    # elcopy = deepcopy(elpt)
    # println("length(elpt.nodes) = ", length(elpt.nodes))
    # println("length(elcopy.nodes) = ", length(elcopy.nodes))

    # println(elpt)
    # # test = push(test,elpt)
    # push!(test,elpt)
    # println("length(test.nodes) = ", length(test.nodes))

    # println(test)
    # # test = push(test,elpt)
    # push!(test,elpt)
    # println("length(test.nodes) = ", length(test.nodes)) 
    # println(test)


    # aim2 = Vector([0.0,1.0,0.0])
    # aim2 = normalize(aim2)
    # pos2 = usistrip.(pos)
    # mom2 = momentum_si(usistrip(en),el)
    # mom2 = aim2*mom2
    # coords2 = Coords(ArrayPartition(pos2,mom2))

    # println("*** debug 1 ***")
    # function make_pt()
    #     # parameters to test build 1
    #     el = ParticleID(0,0,1)
    #     en = upreferred(1.0u"GeV")
    #     pos = Vector([-8.5,0.0,0.0]*u"kpc")
    #     aim = Vector([1.0,0.0,0.0])
    #     dt = 3.6e3u"s"

    #     ### build 1 ###
    #     return ParticleTree(el,en,pos,aim,dt)
    # end

    # pt = make_pt()
    # println("*** debug 2 ***")

    # # push! test
    # push!(pt,coords2)

    # println(typeof(pt[1,1,1,1]) == typeof(coords2))

    # println("*** debug 3 ***")

    # return pt
    return nothing
end

export run

end

import .TestStuff

a = TestStuff.run();
