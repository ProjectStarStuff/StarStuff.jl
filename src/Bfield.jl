################################################################################
# Bfield module
################################################################################
module Bfield
using Unitful
using UnitfulAtomic
using UnitfulAstro
using LinearAlgebra
using Measurements
using IterTools
using Measurements
using LinearAlgebra

function constant(x::Array{T};
    direction::Array = [0.0, 0.0, 6.0e-10]) where {T <: Number}
    out = vec(direction)
end #constx

"""
    dipole(x::Number,y::Number,z::Number; m = [0.,0.,1.1e59])

Returns dipole field given a position and magnetic moment.
Assumes units are in SI.
"""
function dipole(x::Array{T}; m::Array = [0.,0.,1.1e59]) where {T <: Number}

    # println("pos = ", pos)
    rMag = norm(x)
    # println("pos = ", rMag)
    if rMag == 0
        rMag = 1.0e-6
    end
    rHat = x/rMag
    # println("rHat = ", rHat)
    μ0 = 1u"μ0"
    μ0 = μ0.val
    num = μ0*(3.0*rHat*(m⋅rHat)-m)
    den = 4*π*rMag^3
    out = num/den

end #dipole

g = 1e-4u"T"
μg = 1e-6*g

janssonDict = Dict(
    :bi    =>[measurement(0.1  , 1.8 )*μg,
              measurement(3.0  , 0.6 )*μg,
              measurement(-0.9 , 0.8 )*μg,
              measurement(-0.8 , 0.3 )*μg,
              measurement(-2.0 , 0.1 )*μg,
              measurement(-4.2 , 0.5 )*μg,
              measurement(0.   , 1.8 )*μg,
              measurement(2.7  , 1.8 )*μg
             ],
    :bring => measurement(0.1  , 0.1 )*μg,
    :hdisk => measurement(0.4  , 0.03)*u"kpc",
    :wdisk => measurement(0.27 , 0.08)*u"kpc",
    :bn    => measurement(1.4  , 0.1 )*μg,
    :bs    => measurement(-1.1 , 0.1 )*μg,
    :rn    => measurement(9.22 , 0.08)*u"kpc",
    :rs    => measurement(16.7 , 0.0 )*u"kpc",
    :wh    => measurement(0.2  , .12 )*u"kpc",
    :z0    => measurement(5.3  , 1.6 )*u"kpc",
    :bx    => measurement(4.6  , 0.3 )*μg,
    :Θ0x   => measurement(49.  , 1.  )*u"°",
    :rcx   => measurement(4.8  , 0.2 )*u"kpc",
    :rx    => measurement(2.9  , 0.1 )*u"kpc",
    :γ     => measurement(2.92 , 0.14),
    :rnx   => [5.1,6.3,7.1,8.3,9.8,11.4,12.7,15.5]*u"kpc",
    :rinner  => 1.0u"kpc",
    :rdisklims => [3.0,5.0,20.0]*u"kpc", #[molecular ring, start arms, stop arms]
    :θarms => 11.5u"°"
    )

export janssonDict

"""
    spiral(ϕ,rnx,b)

Returns ``r(ϕ) = rnx*exp(b*ϕ) ``
"""
spiral(ϕ,rnx,b) = rnx*exp(b*ϕ)
export spiral

"""
    galspiral(ϕ, params::Dict)

# Returns shifted logarithmic spiral for jansson2012
"""
function galspiral(ϕ, rnx, θarms = janssonDict[:θarms])
    b = tan(θarms)
    return spiral(ϕ-π,rnx,b)
end
export galspiral

logistic(z,h,w) = 1.0/(1.0 + exp(-2.0*(uconvert(Unitful.NoUnits,(abs(z)-h)/w))))
export logistic


function jansson2012(x::Array{T}; params::Dict = janssonDict) where {T <: Quantity}
    @assert typeof(x[1]) <: Unitful.Length "x-coordinate requires units of distance"
    # define μgauss in this way so it plays nicely with various notebooks (Pluto, Jupyter...)
    g = 1e-4u"T"
    μg = 1e-6*g

    # initialize field output as zero
    out = zeros(3)*μg

    # ensure that coordinates are in kpc
    x = map(d->uconvert(u"kpc",d),x)
    ρ = norm(x)

    # return zero field if inside the unfitted center sphere    
    if ρ < params[:rinner]
        return out
    end

    # calculate cylindrical coordinates
    r = sqrt(x[1]^2+x[2]^2)

    # return zero field if outside the galaxy
    if r > params[:rdisklims][end]
        return out
    end

    # define coordinate array with proper units store cylindrical coordinates
    ϕ = x[1].val == 0.0 ? 0.0u"rad" : atan(x[2],x[1])*u"rad"
    z = x[3]

    # generate cylindrical unit vectors
    rhat = [cos(ϕ),sin(ϕ),0.0]
    ϕhat = [-sin(ϕ),cos(ϕ),0.0]
    zhat = [0.0,0.0,1.0]

    # DISK COMPONENT
    if r > params[:rdisklims][2] # arm region
        # determine spiral arm in which point resides

        #   TODO: spiral search could be faster generate spiral transition
        #   radii at given ϕ
        """
            Returns index of spiral in which point resides
        """
        function whichspiral(r,ϕ,rnx = params[:rnx],θarms = params[:θarms])
            for ϕ0 in [-2*π,0.0,2*π,4*π]
                for (i,r0) in enumerate(rnx)
                    rSpiral = galspiral(ϕ+ϕ0,r0,θarms)
                    if r < rSpiral
                        return i
                    end
                end
            end 
            @assert true "There is a bug in whichspiral"
        end

        # used to get field magnitude later
        iSpiral = whichspiral(r,ϕ)
        # disk field direction
        bhat = sin(params[:θarms])*rhat + cos(params[:θarms])*ϕhat
        # the spiral arm disk field scales as 1/r
        bdisk = Measurements.value(params[:bi][iSpiral])*uconvert(Unitful.NoUnits,params[:rdisklims][2]/r)*bhat
        # use logistic function to smoothly transition to halo field 
        bdisk *= 1.0 - logistic(z,Measurements.value(params[:hdisk]),Measurements.value(params[:wdisk]))
        # superpose disk field onto total field
        out += bdisk

    elseif r > params[:rdisklims][1] # molecular ring
        # molecular ring is purely toroidal
        bdisk = Measurements.value(params[:bring])*ϕhat
        # use logistic function transition
        bdisk *= 1.0 - logistic(z,Measurements.value(params[:hdisk]),Measurements.value(params[:wdisk]))
        # superpose disk field onto total field
        out += bdisk        
    end

    # TOROIDAL HALO COMPONENT
    bhmag = exp(-1*uconvert(Unitful.NoUnits,abs(z)/Measurements.value(params[:z0])))
    bhmag *= logistic(z,Measurements.value(params[:hdisk]),Measurements.value(params[:wdisk]))
    if z.val >= 0.0
        bhmag *= Measurements.value(params[:bn])*(1.0 
            - logistic(z,Measurements.value(params[:rn]),Measurements.value(params[:wh])))
    else
        bhmag *= Measurements.value(params[:bs])*(1.0 
            - logistic(z,Measurements.value(params[:rs]),Measurements.value(params[:wh])))
    end
    out += bhmag*ϕhat

    # X FIELD COMPONENT

    # get values from parameters for easier reading
    # TODO: substitute these back in after some readability adjustments
    Θ0x = Measurements.value(params[:Θ0x])
    rx = Measurements.value(params[:rx])
    bx0 = Measurements.value(params[:bx])
    rcx = Measurements.value(params[:rcx])

    rp = r - abs(z)/tan(Θ0x)
    rp = rp.val < 0.0 ? 0.0*u"kpc" : rp
    if rp > rcx
        bx = bx0*exp(uconvert(Unitful.NoUnits,-1*rp/rx))
        bmag = bx*uconvert(Unitful.NoUnits,rp/r)
        bhat = z.val > 0.0 ? [cos(ϕ)*cos(Θ0x),sin(ϕ)*cos(Θ0x),sin(Θ0x)] : [cos(ϕ)*cos(Θ0x),sin(ϕ)*cos(Θ0x),-1*sin(Θ0x)]
        out += bmag*bhat
    else
        rp = r*rcx/(rcx+abs(z)/tan(Θ0x))
        bx = bx0*exp(uconvert(Unitful.NoUnits,-1*rp/rx))
        bmag = bx*uconvert(Unitful.NoUnits,rp/r)^2
        Θx = r == rp ? 0.0u"rad" : atan(abs(z)/(r-rp))
        bhat = z.val > 0.0 ? [cos(ϕ)*cos(Θx),sin(ϕ)*cos(Θx),sin(Θx)] : [cos(ϕ)*cos(Θx),sin(ϕ)*cos(Θx),-1*sin(Θx)]
        out += bmag*bhat
    end

    return out
end #jansson2012

export jansson2012
########################### B-field wrappers ##################################
function galacticDipole(x::Array{T}) where {T <: Number}
    dipole(x)
end #galacticDipole

function galacticConstantZ(x::Array{T}) where {T <: Number}
    constant(x)
end #galacticConstantX

function jansson2012m(x::Array{T}) where {T <: Number}
    xin = uconvert.(u"kpc",x*u"m")
    return ustrip.(u"T",jansson2012(xin)) 
end

function jansson2012kpc(x::Array{T}) where {T <: Number}
    xin = x*u"kpc"
    return ustrip.(u"T",jansson2012(xin))
end
######################## end B-field wrappers ##################################

"""
`Bfield.model`: dictionary containing commonly used magnetic fields.
"""
model = Dict(
    :galacticDipole         => galacticDipole,
    :galacticConstantZ      => galacticConstantZ,
    :jansson2012m           => jansson2012m,
    :jansson2012kpc         => jansson2012kpc
)

end #Bfield

export Bfield
