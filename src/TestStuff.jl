module TestStuff

include("StarStuff.jl")
using .StarStuff
# using Gtk
# js"""
# function(){
#     function clearactive(item,index){
        

#     }
#     let thisindex = 1;
#     let tablist = ['particles','bfield','sim','plot'];
#     let thisid = tablist[thisindex];
#     let thistab = document.getElementById(thisid);
#     for (let iname in tablist){
#         let itab = document.getElementById(iname);
#         itab.classList.remove('is-active');
#     }

#     alert('test');


#     thistab.classList.add('is-active');
# }
# """


function run()
    a = Vector([1.,2.,3.])
    b = Vector([4.,5.,6.])

    c = ArrayPartition(a,b)
    d = Coords(c)
    println(d)

    # win = GtkWindow("My First Gtk.jl Program", 400, 200)

    # b = GtkButton("Click Me")
    # push!(win,b)
    
    # function on_button_clicked(w)
    #   println("The button has been clicked")
    # end
    # signal_connect(on_button_clicked, b, "clicked")
    
    # Gtk.showall(win)
    # appdir = joinpath(@__DIR__,"..","app")
    # bulma  = joinpath(appdir,"bulma.css")
    # uijs   = joinpath(appdir,"ui.js")
    # WNode  = WebIO.Node
    # w = Window(async=false)
    # load!(w,bulma)
    # load!(w,uijs)
    # tabitems = WNode(
    #     :ul,
    #     WNode(:li,
    #         WNode(:a,"Particles"),
    #         id="particles",
    #         # className="is-active",
    #         events=Dict(
    #             # "click" => js"tabs.setActive(1)"
    #         )
    #     ),
    #     WNode(:li,
    #         WNode(:a,"Magnetic Field"),
    #         id="bfield",
    #         className="is-active",
    #         events=Dict(
    #             # "click" => js"tabs.setActive(2)"
    #         )

    #     ),
    #     WNode(:li,
    #         WNode(:a,"Simulate"),
    #         id="sim",
    #         events=Dict(
    #             # "click" => js"tabs.setActive(3)"
    #         )

    #     ),
    #     WNode(:li,
    #         WNode(:a,"Plot"),
    #         id="plot",
    #         events=Dict(
    #             # "click" => js"tabs.setActive(4)"
    #         )

    #     )

    # )
    # tabs = WNode(
    #     :div,
    #     tabitems,
    #     attributes=Dict(:class=>"tabs is-boxed",:id=>"star-tabs")
    # )


    # body!(w,tabs)

    # ********** ParticleTree widget test **********
    # pt = pt_manualdt_wdg()

    # # printbutton = Interact.button("Print")
    # # Interact.@on println("Number of particles = ",length(getproperty(&pt,:nodes)))
    # Interact.@on println("*** New Particle ***\n",getindex(&pt,2))
    # # Interact.@map (&printbutton; println(pt[]))

    # w = Window()
    # ui = dom"div"(
    #     pt
    # )
    # body!(w,ui)

    # ********** Creation of empty ParticleTree and appending nodes **********
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
