
"""
	clear!(scene::Scene)

Removes everything from scene except axes
"""
function clear!(scene::Scene) 
	if length(scene.plots) > 2
		delete!(scene,scene.plots[end])		
		clear!(scene)
	elseif length(scene.plots) == 2
		delete!(scene,scene.plots[end])		
	end
end

export clear!

"""
	scaleShift(a::T) where {T<:AbstractVector}

Shifts values of array so it is centered on its mean and range is normalized to 1
"""
function scaleshift(a::T) where {T<:AbstractVector}
	aRange = max(a...) - min(a...)
	aMean  = mean(a)
	out = a .- aMean
	if aRange == 0.0
		return out
	else
		return aMean, aRange, out ./ aRange
	end
end

export scaleshift