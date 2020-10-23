#%%
using DrWatson
@quickactivate "crprop"

using Plots
using StatsPlots
using Unitful
using UnitfulAstro
using UnitfulAtomic
using LinearAlgebra
using RecursiveArrayTools
using DifferentialEquations
using Random
using Printf

include(srcdir("Bfield.jl"))
include(srcdir("Kinematics.jl"))
include(srcdir("RelativisticIntegrator.jl"))
include(srcdir("PropUtils.jl"))

using .Bfield
using .Kinematics
using .RelativisticIntegrator
using .PropUtils

# particle mass and charge definitions
el = cr(0,0,1)
pr = cr(1,0,0)

# energies to test
e0 = [1.0, 10.0, 100.0, 1000.0] * u"GeV"

# Calculate magnitude of momentum given energy and rest mass and assign to the y-direction

# for electrons
pe0 = Kinematics.momentum.(e0, mcr(el))
pe0 = [
    vec([0.0 * unit(pe0[i]), pe0[i], 0.0 * unit(pe0[i])]) for i = 1:length(pe0)
]

seed = 1234
rng = MersenneTwister(seed)

randpm1(rng) = 2.0*rand(rng)-1
randθ(rng) = π*rand(rng)
randϕ(rng) = 2.0*π*rand(rng)

# give a random momentum to particles in galacticDipole field
for i=1:length(e0)
    θ = randθ(rng)
    ϕ = randϕ(rng)
    momentum = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]*norm(pe0[i])
    push!(pe0, momentum)
end


# set initial position to be one gyroradius (given a 6μG field) in the x-direction

# for electrons
rge = Kinematics.gyroradius.(norm.(pe0), abs(qcr(el)), 6.0e-10u"T")
x0 = [
    vec([rge[i], 0.0 * unit(rge[i]), 0.0 * unit(rge[i])]) for i = 1:length(e0)
]

# initial position at (-8.5, 0, 0) kpc for second half
earthPos = uconvert(u"m",-8.5u"kpc")
for i = 1:4
    push!(x0,vec([earthPos,0.0u"m",0.0u"m"]))
end


nx = [PropUtils.usistrip.(x0[i]) for i = 1:length(x0)]
npe = [PropUtils.usistrip.(pe0[i]) for i = 1:length(pe0)]

q = PropUtils.usistrip(qcr(el))
m = PropUtils.usistrip(mcr(el))

efield(x::Array{T}) where {T<:Number} = vec([0.0, 0.0, 0.0])
bfield1 = Bfield.model[:galacticConstantZ]
bfield2 = Bfield.model[:galacticDipole]


# define initial conditions
dicts = []
params1 = Dict(
        :t => 0.0,
        :q => q,
        :m => m,
        :efield => efield,
        :bfield => bfield1
)
params2 = Dict(
        :t => 0.0,
        :q => q,
        :m => m,
        :efield => efield,
        :bfield => bfield2
)

gyroPeriod = 2*π*PropUtils.usistrip.(rge) ./ norm.(Kinematics.vFromP.(npe,mcr(el).val))
for i=1:length(x0)
    newDict = Dict(
        :prob => DiscreteProblem(rintegrator[:boris],ArrayPartition(nx[i],npe[i]),(1.0,gyroPeriod[i]*100),i < 5 ? params1 : params2),
        :gp => gyroPeriod[i],
        :sol => -1,
        :evaltime => -1
    )
    push!(dicts,newDict)
end



# solve problems in dicts
for i=1:length(dicts)
    solution = @timed solve(dicts[i][:prob], FunctionMap(), dt = dicts[i][:gp]/100)
    dicts[i][:sol] = solution.value
    dicts[i][:evaltime] = solution.time
end
dicts[1][:sol].u

# plot everything
for i=1:length(dicts)
    out = dicts[i][:sol]

    ename = string(e0[(i-1)%4+1])
    eval = e0[(i-1)%4+1].val
    eval = @sprintf "%.1E" eval

    dt = out.t[2]-out.t[1]
    dt = @sprintf "%.2E" dt
    if i < 5
        plot(out,vars=(1,2),ratio=:equal)#,linecolor=:viridis,line_z=out.t)
        plot!(title="$(ename) electron, 6μG B-field, dt=$(dt) sec")
        plot!(colorbar_title="Time (seconds)")
        xlabel!("Distance (m)")
        ylabel!("Distance (m)")


        savefig(plotsdir("test-integrators","boris","e-=$(eval)GeV_b=6uG_Boris_dt=$(dt).png"))
    else
        plot(out,vars=(1,2,3),linecolor=:viridis,line_z=out.t)
        plot!(title="$(ename) electron, Galactic Dipole, dt=$(dt) sec")
        #plot!(colorbar_title="Time (seconds)")


        savefig(plotsdir("test-integrators","boris","e-=$(eval)GeV_b=dipole_Boris_dt=$(dt).png"))
    end
end

function bfield3(x::Array{T}) where {T <: Number}
    return Bfield.dipole(x,m = [0.,0.,1.1e28])
end

bfield3(nx[1])


params3 = Dict(
        :t => 0.0,
        :q => q,
        :m => m,
        :efield => efield,
        :bfield => bfield3
)

newDict = Dict(
    :prob => DiscreteProblem(rintegrator[:boris],ArrayPartition(nx[1],npe[5]),(1.0,gyroPeriod[1]*8), params3),
    :gp => gyroPeriod[1],
    :sol => -1,
    :evaltime => -1
)

solution = @timed solve(newDict[:prob], FunctionMap(), dt = newDict[:gp]/1000)
newDict[:sol] = solution.value
newDict[:evaltime] = solution.time

plot(newDict[:sol],vars=(1,2,3), camera=(50,20), ratio=:equal)
plot!(size=(800,600))
anim = Animation()
for i=0:40
    plot!(camera=(abs(20-i)+60,30))
    frame(anim)
end
filename = plotsdir("test-integrators","boris","e-=1.0GeV_mdipole=1.1e28JT-1_Boris_t=100gp_dt=gpE-3_longer-run.gif")

gif(anim, filename, fps=15)
#%%
plot(newDict[:sol],vars=(1,2,3), camera=(90,90), ratio=:equal)
plot!(size=(800,600))

filename = plotsdir("test-integrators","boris","e-=1.0GeV_mdipole=1.1e28JT-1_Boris_t=100gp_dt=gpE-3_top.png")

savefig(filename)
