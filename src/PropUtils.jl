"""
	function mag10(x::Float64)

Returns the base 10 order of magnitude of x
"""
function mag10(x::Float64)
	if x == 0.0
		return 0
	end
	l10 = log10(x)
	out = trunc(l10)
	if l10 - out == 0.0
		return Int(out)
	end
	if l10 >= 0.0
		return Int(out)
	end
	return Int(out) - 1
end

export mag10

#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
    usistrip(val)

Converts units to SI, then strips them
"""
usistrip(val) = ustrip(upreferred(val))

export usistrip
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
    getArrayPartition(i,sol)

Returns the ith partition from an ArrayPartition
"""
getArrayPartition(i, sol) = sol.x[i]

export getArrayPartition
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
"""
    getPosition(sol; units = u"m")

Returns the position coordinates of a solved trajectory
"""
function getPosition(sol; units = u"m")
    x = getArrayPartition.(2, sol)
    x = hcat(x...)
    x = x * u"m"
    x = uconvert.(units, x)
    return ustrip(x)
end

export getPosition
#..oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo...oooOOO0Q0Q0OOOooo
