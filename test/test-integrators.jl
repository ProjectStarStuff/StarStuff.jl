#%%
using DrWatson
@quickactivate "crprop"

using Plots
using StatsPlots
using Unitful
using UnitfulAstro
using UnitfulAtomic
using LinearAlgebra
using DifferentialEquations
using Random
#%%
include(srcdir("Bfield.jl"))
include(srcdir("Proptools.jl"))
#%%
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
Differential equation for position
    * Use SI units for now
"""
function dxdt!(dx, p, x, params, t)
    dx .= Proptools.vFromP(p,params[1])
end
#%%
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
Differential equation for momentum
    * Use SI units for now
"""
function dpdt!(dp, p, x, params, t)
    b = Bfield.model[Symbol(params[3])](x...)
    dp .= params[2]*cross(Proptools.vFromP(p,params[1]),b)
end
#%%
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
usistrip(val) = ustrip(upreferred(val))
#%%
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo

"""
    getArrayPartition(i,sol)

Returns the ith partition from an ArrayPartition
"""
getArrayPartition(i,sol) = sol.x[i]
#%%
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
    getPosition(sol; units = u"m")

Returns the position coordinates of a solved trajectory
"""
function getPosition(sol; units = u"m")
    x = getArrayPartition.(2,sol)
    x = hcat(x...)
    x = x*u"m"
    x = uconvert.(units,x)
    return ustrip(x)
end #getx
#%%
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo

me = 511.0u"keV/c^2"
mp = 938.0u"MeV/c^2"
qp = 1u"q"
qe = -1u"q"

e0 = [1.0, 10.0, 100.0, 1000.0]*u"GeV"
#%%
# TODO: Change pe0 (and others) to an array of vectors
pe0 = Proptools.momentum.(e0,me)
pe0 = [vec([0.0*unit(pe0[i]) pe0[i] 0.0*unit(pe0[i])]) for i=1:length(pe0)]
#pp0 = Proptools.momentum.(e0,mp)
#%%
rge = Proptools.gyroradius.(norm.(pe0),qp,6.0e-10u"T")
ve0 = Proptools.vFromP.(pe0, me)
#vp0 = Proptools.vFromP.(pp0, mp)
x0 = [vec([rge[i] 0.0*unit(rge[i]) 0.0*unit(rge[i])]) for i=1:length(rge)]
#%%
nx0 = [usistrip.(x0[i]) for i=1:length(x0)]
npe0 = [usistrip.(pe0[i]) for i=1:length(pe0)]
#%%

# γe = usistrip(e0/(me*1u"c^2"))
# γp = usistrip(e0/(mp*1u"c^2"))

# βe = β(γe)
# βp = β(γp)
#%%

# create dictionary for input parameters
rkDict = Dict(
    :params => [(upreferred(me).val, upreferred(qe).val,"galacticConstantZ")],
    :x0 => [nx0[1]],
    :p0 => [npe0[1]],
    :integrator => [Tsit5],
    :tspan => [(0.0, 1.0e4)],
    :dt => [nothing],
    :units => "SI",
    :sol => [nothing],
    :evaltime => [nothing]
)

symplecticDict = Dict(
    :params => [(upreferred(me).val, upreferred(qe).val,"galacticConstantZ")],
    :x0 => [nx0[1]],
    :p0 => [npe0[1]],
    :integrator => [Yoshida6, Nystrom4,DPRKN6,McAte5],
    :tspan => [(0.0, 1.0e4)],
    :dt => [10.0],
    :units => "SI",
    :sol => [nothing],
    :evaltime => [nothing]
)

dicts = dict_list(rkDict)
dicts2 = dict_list(symplecticDict)

append!(dicts, dicts2)
#%%
function prop!(d::Dict)
    println("**************** Constructing Problem ******************")

    prob = DynamicalODEProblem(
    dpdt!,
    dxdt!,
    d[:p0],
    d[:x0],
    d[:tspan],
    d[:params]
    )

    println(prob)
    println("**************** Beginning Propagation ******************")
    sol = nothing

    if d[:dt] == nothing
        sol = @timed solve(prob,d[:integrator](),maxiters=1e5)#,dt=d[:dt])
    else
        sol = @timed solve(prob,d[:integrator](),dt = d[:dt], maxiters=1e5)
    end

    println("**************** Finished Propagation ******************")
    d[:sol] = sol.value
    d[:evaltime] = sol.time

end

#%%
for (i, d) in enumerate(dicts)
    prop!(d)
end
#%%

for (i, d) in enumerate(dicts)
    out = dicts[i][:sol]
    dt = dicts[i][:dt]
    npts = 10001
    tInterp = collect(range(min(out.t...),max(out.t...),length=npts))
    interp = out.(tInterp)
    #%%
    # get x coordinates in kilometers
    x = getPosition(interp,units=u"km")
    #%%

    xkm = ustrip.(x[1,:])
    ykm = ustrip.(x[2,:])

    iname = string(nameof(d[:integrator]))
    #%%
    plot(xkm,ykm,label=iname,ratio=:equal,linecolor=:viridis,line_z=tInterp)
    #plot(out,vars=(4,5),label="Tsit5",alpha=0.5)

    rg = uconvert(u"km",rge[1])

    rgkm = ustrip(rg)

    θ = collect(range(0.0,2*π,length=100))

    #plot!(rgkm*cos.(θ), rgkm*sin.(θ),color=:red,label="Actual",linewidth=3)

    if d[:dt] == nothing
        plot!(title="1 GeV electron, 6μG B-field, adaptive time step")
    else
        plot!(title="1 GeV electron, 6μG B-field, dt=$dt")
    end

    plot!(colorbar_title="Time (seconds)")

    xlabel!("Distance (km)")
    ylabel!("Distance (km)")

    if d[:dt] == nothing
        savefig(plotsdir("test-integrators","e-=1GeV_b=6uG_$(iname)_dt=adaptive_no-actual.png"))
    else
        savefig(plotsdir("test-integrators","e-=1GeV_b=6uG_$(iname)_dt=$(dt)_no-actual.png"))
    end

end
