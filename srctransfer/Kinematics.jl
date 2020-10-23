################################################################################
# Convenience functions for kinematics
################################################################################

#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo

"""
	β(γ::Number)

Returns β given γ.
"""
β(γ::Number) = √(γ^2 - 1) / γ

export β
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
    function momentum(E::Quantity,id::ParticleID)

Returns momentum given total energy and rest mass.
"""
function momentum(E::Quantity, id::ParticleID)
    @assert typeof(E) <: Unitful.Energy "E requires units of energy"
    m = mass(id)
    return uconvert(u"GeV/c",√(E^2 - m^2 * Unitful.c^4) / Unitful.c)
end #momentum

"""
	function momentum(E::Quantity,m::Quantity)

Returns momentum given total energy and rest mass.
"""
function momentum(E::Quantity, m::Quantity)
    @assert typeof(E) <: Unitful.Energy "E requires units of energy"
    @assert typeof(m) <: Unitful.Mass "m requires units of mass"
    return uconvert(u"GeV/c",√(E^2 - m^2 * Unitful.c^4) / Unitful.c)
end #momentum


"""
	function momentum(E::Number,m::Number)

Returns momentum given total energy and rest mass in SI units.
"""
function momentum(E::Number, m::Number)
    return √(E^2 - m^2 * Unitful.c0.val^4) / Unitful.c0.val
end #momentum

export momentum
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo

"""
	γFromβ(β::Number)

Returns γ given β
"""
function γFromβ(β::Number)
    @assert β <= 1.0 "β ≦ 1 required"
    return 1 / sqrt(1 - β^2)
end #γFromβ

#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo

"""
    γFromP(p::Array{T}, m::Quantity) where {T<:Quantity}

Returns gamma given momentum and rest mass.
"""
function γFromP(p::Array{T}, m::Quantity) where {T<:Quantity}
    pnorm = norm(p)
    γFromP(pnorm,m)
end #γFromP

"""
    γFromP(p::Array{T}, m::Number) where {T<:Number}

Returns gamma given momentum and rest mass in SI units.
"""
function γFromP(p::Array{T}, m::Number) where {T<:Number}
    pnorm = norm(p)
    γFromP(pnorm,m)
end #γFromP

"""
γFromP(p::Quantity,m::Quantity)

Returns gamma given momentum and rest mass.
"""
function γFromP(p::Quantity, m::Quantity)
    @assert typeof(p) <: Unitful.Momentum "p requires units of momentum"
    @assert typeof(m) <: Unitful.Mass "m requires units of mass"
    c = Unitful.c
    βγ = upreferred(p / (m * c))
    return upreferred(sqrt(1 + βγ^2))
end #γFromP


"""
	γFromP(p::Number,m::Number)

Returns gamma given momentum and rest mass in SI units.
"""
function γFromP(p::Number, m::Number)
    c = Unitful.c0.val
    return sqrt(1 + (p / (m * c))^2)
end #γFromP

export γFromP
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
	βFromP(p::Quantity,m::Quantity)

Returns β given particle momentum and rest mass.
"""
function βFromP(p::Quantity, m::Quantity)
    @assert typeof(p) <: Unitful.Momentum "p requires units of momentum"
    @assert typeof(m) <: Unitful.Mass "m requires units of mass"
    c = Unitful.c
    η = upreferred(p / (m * c))
    return sqrt(η^2 / (1 + η^2))
end


"""
	βFromP(p::Number,m::Number)

Returns β given particle momentum and rest mass in SI units
"""
function βFromP(p::Number, m::Number)
    c = Unitful.c0.val
    η = p / (m * c)
    return sqrt(η^2 / (1 + η^2))
end

export βFromP
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
    vFromP(p::Array{T}, m::U) where {T <: Quantity, U <: Quantity}

Returns velocity given particle momentum and rest mass in SI units
"""
function vFromP(p::Array{T}, m::U) where {T <: Quantity, U <: Quantity}
    β = βFromP(norm(p),m)
    c = Unitful.c
    phat = p / norm(p)
    return upreferred.(c * β * phat)
end


"""
    vFromP(p::Array{T}, m::T) where {T <: Number}

Returns velocity given particle momentum and rest mass in SI units
"""
function vFromP(p::Array{T}, m::T) where {T <: Number}
    β = βFromP(norm(p),m)
    c = Unitful.c0.val
    phat = p / norm(p)
    return c * β * phat
end

"""
    vFromP(p::Quantity, m::Quantity)

Returns velocity given particle momentum and rest mass in SI units
"""
function vFromP(p::Quantity, m::Quantity)
    β = βFromP(p,m)
    return upreferred.(Unitful.c * β )
end


"""
    vFromP(p::Array{T}, m::T) where {T <: Number}

Returns velocity given particle momentum and rest mass in SI units
"""
function vFromP(p::Number, m::Number)
    β = βFromP(p,m)
    c = Unitful.c0.val
    return c * β
end

export vFromP
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
	gyroradius(p::Quantity,Z::Quantity,B::Quantity)

Returns gyroradius given dimensionally correct momentum, charge, and B field strength
"""
function gyroradius(p::Quantity, Q::Quantity, B::Quantity)
    @assert typeof(p) <: Unitful.Momentum "p must have units of momentum"
    @assert typeof(Q) <: Unitful.Charge "Z must have units of charge"
    @assert typeof(B*1u"m^2") <: Unitful.MagneticFlux "B must have units of magnetic flux density"

    return upreferred(p/(Q*B))
end
# TODO: add functionality to allow for vector p and B

"""
	gyroradius(p::Number, Z::Number, B::Number)

Returns gyroradius (m) given momentum, atomic (charge) number, and B field strength

# Arguments
- `p::Number`: momentum of particle in GeV/c
- `Z::Number`: charge of particle in units of q_e
- `B::Number`: magnetic flux density in gauss
"""
function gyroradius(p::Number, Z::Number, B::Number)
    33356.40951981521 * p / (Z * B)
end

export gyroradius
