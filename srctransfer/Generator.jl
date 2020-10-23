"""
    unitRandom3Cartesian(rng = MersenneTwister())

Returns a unit vector pointing in a random direction in Cartesian coordinates.
"""
function unitRandom3Cartesian(rng = MersenneTwister())
    θ = π*rand(rng)
    ϕ = 2.0*π*rand(rng)
    return Vector([
        cos(ϕ)*sin(θ),
        sin(ϕ)*sin(θ),
        cos(θ)
    ])
end

export unitRandom3Cartesian

"""
    randbtw(a::T,b::T,rng=MersenneTwister()) where {T<:Real}

Generates a random number between a and b.
"""
function randbtw(a::T,b::T,rng=MersenneTwister()) where {T<:Real}
    if a > b
        return randbtw(b,a,rng)
    end
    return a+rand(rng)*(b-a)
end

export randbtw
"""
    dist(n,func,domain,range,rng=MersenneTwister())

Generates n values following the distribution func in a given domain and range
    - SLOW
    - Can be used on any distribution function
"""
function dist(n,func,domain,range,rng=MersenneTwister())
    out = zeros(n)
    for i=1:n
        val = 0.0
        tryagain = true
        while tryagain
            xtest = randbtw(domain...,rng)
            fxtest = abs(func(xtest))
            if abs(randbtw(range...,rng)) < fxtest
                val = xtest
                tryagain = false
            end
        end
        out[i] = val
    end
    return out
end

export dist

"""
Sampler for a power law distribution
"""
struct PowerLaw <: Sampleable{Univariate, Continuous}
    α::Real
    domain::Tuple{T,T} where {T<:Real}
end

export PowerLaw

"""
    rand(rng::AbstractRNG,s::PowerLaw)

Generate one sample from a power law distribution
"""
function Base.rand(rng::AbstractRNG,s::PowerLaw)
    ξ = rand(rng)
    γ = 1.0-s.α
    e0γ= s.domain[1]^(γ)
    efγ = s.domain[2]^(γ)
    out = (efγ-e0γ)*ξ+e0γ
    return out^(1/γ)
end

"""
Sampler for a broken power law distribution
"""
struct BrokenPowerLaw <: Sampleable{Univariate, Continuous}
    α::Array{T} where {T<:Real}
    domain::Array{T} where {T<:Real}
    function BrokenPowerLaw(α::Array{T}, domain::Array{T}) where {T<:Real}
        @assert length(α)+1 == length(domain) "Requires length(α) == length(domain)-1"
        bOrder = true
        for i=1:(length(domain)-1)
            bOrder *= domain[i] < domain[i+1]
        end

        @assert bOrder "domain breaks must be in order"
        new(α, domain)
    end
end

export BrokenPowerLaw

function afterIndex(val,vec)
	for i=1:(length(vec)-1)
		check = vec[i] < val
		check *= val < vec[i+1]
		if check
			return i
		end
	end
	return 0
end

function Base.rand(rng::AbstractRNG,s::BrokenPowerLaw)
	ξ = rand(rng)
	γ = 1 .- s.α


	cx = ones(length(s.α))
	for i = 1:(length(cx)-1)
		β = s.α[i]/s.α[i+1]
		cx[i+1] = cx[i]^β*s.domain[i+1]^(β-1)
	end
	cx = map(/,cx,γ)


	arts = [
		begin
			(s.domain[i+1]^γ[i] - s.domain[i]^γ[i])*cx[i]
		end for i=1:length(γ)
			]

	parts = vcat([0],arts)
	sumParts = cumsum(parts)

	norm = sumParts[end]
	cnorm = ξ*norm

	supIndex = afterIndex(cnorm,sumParts)

	out = cnorm - sumParts[supIndex]
	out /= cx[supIndex]
	out += s.domain[supIndex]^γ[supIndex]
	out = out^(1/γ[supIndex])
	return out
end
